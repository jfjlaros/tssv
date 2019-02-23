/*
 * Library for semi-global alignment.
 */
#include <emmintrin.h>
#include <stdlib.h>
#include <string.h>
#include <xmmintrin.h>

#include "sgAlignSSE.h"


/**
 * Calculate the minimum of two values.
 *
 * @arg {char} a - A value.
 * @arg {char} b - A value.
 *
 * @return {char} - The minimum of {a} and {b}.
 */
static inline char _min(char a, char b) {
  if (a < b)
    return a;
  return b;
}

/**
 * Reverse a sequence.
 *
 * @arg {char *} seq - Sequence to reverse.
 * @arg {char *} seqr - Pointer to a buffer that will receive the reversed
 *   sequence.
 * @arg {size_t} len - Length of the sequence, excluding the terminating NULL
 *   byte.
 */
void _revSeq(char *seq, char *seqr, size_t len) {
  char *p = seq, *q = seqr + len;

  *(q--) = 0;
  while (q >= seqr)
    *(q--) = *(p++);
}

/**
 * Initialise a matrix for semi-global alignment.
 *
 * @arg {unsigned int} rows - Number of rows in the matrix.
 * @arg {unsigned int} columns - Number of columns in the matrix. MUST NOT be
 *   less than rows.
 * @arg {unsigned char} indelScore - Penalty score for insertions and
 *   deletions.
 *
 * @return {unsigned char *} - The alignment matrix. The actual matrix starts
 *   at the first 16-byte aligned byte and is stored diagonally.
 */
unsigned char *_makeMatrixSSE(
    unsigned int rows, unsigned int columns, unsigned char indelScore) {
  unsigned int width = ((rows + 14) & ~0x0F) + 16,
               height = columns + rows - 1;
  unsigned char *mem = (unsigned char *)malloc(width * height + 16),
                *matrix = (unsigned char *)(
                  ((unsigned long int)mem + 15) & ~(unsigned long int)0x0F),
                *cell,
                score;
  unsigned int i,
               j;

  // Set the first column to 0.
  for (i = 0, cell = matrix; i < columns; i++, cell += width)
    *cell = 0;

  // Set the first row to 0, 1, 2, 3, 4, ... times the indelScore
  for (i = 0, cell = matrix, score = 0; i < rows; i++, cell += width + 1) {
    *cell = score;
    score = (score > 255 - indelScore)? 255 : score + indelScore;
    for (j = 1; j < (width - i); j++)
      *(cell + j) = 255;  // This protects the second row.
  }

  return mem;
}

/**
 * Fill the alignment matrix.
 *
 * @arg {unsigned char *} mem - The alignment matrix. The actual matrix is
 *   assumed to start at the first 16-byte aligned byte.
 * @arg {unsigned int} rows - Number of rows in the matrix.
 * @arg {unsigned int} columns - Number of columns in the matrix.
 * @arg {char *} seq1 - The first sequence to be aligned.
 * @arg {char *} seq2 - The second sequence to be aligned.
 * @arg {unsigned char} indelScore - Penalty score for insertions and
 *   deletions.
 */
void _alignSSE(
    unsigned char *mem, unsigned int rows, unsigned int columns,
    char *seq1, char *seq2, unsigned char indelScore) {
  unsigned int x = 1,
               y = 1,
               width = ((rows + 14) & ~0x0F) + 16,
               end = columns + _min(15, rows - 1),
               limit;
  unsigned char *matrix = (unsigned char *)(
                  ((unsigned long int)mem + 15) & ~(unsigned long int)0x0F),
                *d = matrix,
                *l = matrix + width,
                *i = l + width + 1;

  const __m128i ones = _mm_set1_epi8(1),
    indelScores = _mm_set1_epi8(indelScore),
    range = _mm_set_epi8(15, 14, 13, 12, 11, 10, 9, 8, 7 ,6, 5, 4, 3, 2, 1, 0);
  __m128i mi,
          md = _mm_load_si128((__m128i *)d),
          ml = _mm_load_si128((__m128i *)l),
          mu = _mm_loadu_si128((__m128i *)(l + 1)),
          mx,
          my;

  // Get copy of seq1 and reverse of seq2, making sure
  // that we can read 16 bytes (of garbage) past the end.
  const size_t seq2len = strlen(seq2);
  char *seq1f = (char *)malloc(strlen(seq1) + 16),
       *seq2r = (char *)malloc(seq2len + 16);
  strcpy(seq1f, seq1);
  _revSeq(seq2, seq2r, seq2len);
  mx = _mm_loadu_si128((__m128i *)seq1f);
  my = _mm_loadu_si128((__m128i *)(seq2r + seq2len - 1));

  while (1) {
    mi = _mm_min_epu8(_mm_adds_epu8(_mm_min_epu8(ml, mu), indelScores),
      _mm_adds_epu8(md, _mm_add_epi8(_mm_cmpeq_epi8(mx, my), ones)));

    limit = y - x + 1;
    if (limit >= 16 || y >= (rows - 1))
      _mm_storeu_si128((__m128i *)i, mi);
    else
      // Need to make sure the top row and left column stay valid.
      // This is VERY MUCH SLOWER than the above.
      _mm_maskmoveu_si128(mi, _mm_cmplt_epi8(
        range, _mm_set1_epi8(limit)), (char *)i);

    if (++y < end) {
      // Move down one row.
      l += width;
      i += width;
      md = ml;
      mu = mi;
      ml = _mm_load_si128((__m128i *)l);

      // Move to the next base in seq2r.
      my = _mm_slli_si128(my, 1);
      if (limit - 1 < columns)
        my = _mm_insert_epi16(my, *(short *)(seq2r + seq2len - limit - 1), 0);
    }
    else if ((x += 16) < rows) {
      // Move right 16 columns.
      y = x;
      end += _min(16, rows - x);
      d += 16 * width + 16;
      l = d + width;
      i = l + width + 1;
      md = _mm_load_si128((__m128i *)d);
      ml = _mm_load_si128((__m128i *)l);
      mu = _mm_loadu_si128((__m128i *)(l + 1));
      mx = _mm_loadu_si128((__m128i *)(seq1f + x - 1));
      my = _mm_loadu_si128((__m128i *)(seq2r + seq2len - 1));
    }
    else break;
  }

  free(seq1f);
  free(seq2r);
}

/**
 * find the minimum distance, ignoring a trailing gap in the sequence
 * associated with the number of rows in an alignment matrix. if the minimum
 * distance is found, also return the row number.
 *
 * it is assumed that the number of rows is larger than the number of columns.
 *
 * @arg {unsigned char *} mem - an {rows} * {columns} matrix. the actual matrix
 *   is assumed to start at the first 16-byte aligned byte.
 * @arg {unsigned int} rows - number of rows in the matrix.
 * @arg {unsigned int} columns - number of columns in the matrix. must not be
 *   less than rows.
 *
 * @return {alignment} - the minimum distance and its row number.
 */
alignment _findMinSSE(
    unsigned char *mem, unsigned int rows, unsigned int columns) {
  unsigned int width = ((rows + 14) & ~0x0f) + 16,
               i;
  unsigned char *matrix = (unsigned char *)(
                  ((unsigned long int)mem + 15) & ~(unsigned long int)0x0f),
                *cell;
  alignment a;

  a.position = 0;
  a.distance = rows - 1;
  cell = matrix + (rows - 1) * width + rows - 1;
  for (i = 0; i < columns; i++, cell += width)
    if (*cell < a.distance) {
      a.distance = *cell;
      a.position = i;
    }

  return a;
}

/**
 * do a semi-global alignment of {seq2} to {seq1}.
 *
 * @arg {char *seq1} - the sequence to be aligned to.
 * @arg {char *seq2} - the sequence to be aligned.
 * @arg {unsigned char} indelScore - penalty score for insertions and
 *   deletions (max 255).
 *
 * @return {alignment} - the minimum distance and its row number.
 */
alignment alignSSE(char *seq1, char *seq2, unsigned char indelScore) {
  alignment a;
  unsigned char *matrix;
  unsigned int rows = strlen(seq1) + 1,
               columns = strlen(seq2) + 1;

  // the alignment is optimised for tall matrices.
  matrix = _makeMatrixSSE(columns, rows, indelScore);
  _alignSSE(matrix, columns, rows, seq2, seq1, indelScore);
  a = _findMinSSE(matrix, columns, rows);
  free(matrix);
  return a;
}
