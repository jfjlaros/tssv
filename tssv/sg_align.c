/*
Library with functions for a semi-global alignment.


Two implementations are contained within this file:
1. An implementation using SSE2 instructions, which is used on systems that
   support these instructions. The SSE2 implementation uses a matrix of
   unsigned chars instead of ints and stores the matrix diagonally. The
   implementation is optimised for tall matrices.
   This implementation achieves a 3-4x performance improvement compaired to
   TSSV 0.2.5, on data consiting predominantly of reads of 200-300bp.
2. A fallback implementation which is used if SSE2 is not available.
*/

#include <string.h>
#include <stdlib.h>
#include "sg_align.h"

/*
Calculate the minimum of two values.

:arg char a: A value.
:arg char b: A value.

:returns char: The minimum of {a} and {b}.
*/
static __inline char _min(char a, char b) {
  if (a < b)
    return a;
  return b;
}

#if defined(__GNUC__) && defined(__SSE2__)
/******************************************************************************
  SSE2-enabled Implementation
******************************************************************************/
#  include <xmmintrin.h>
#  include <emmintrin.h>

/*
Initialise a matrix for semi-global alignment.

:arg unsigned int rows: Number of rows in the matrix.
:arg unsigned int columns: Number of columns in the matrix. MUST NOT be
  less than rows.
:arg unsigned char indel_score: Penalty score for insertions and deletions.

:returns unsigned char *: The alignment matrix. The actual matrix starts at the
  first 16-byte aligned byte and is stored diagonally.
*/
unsigned char *_make_matrix(
    unsigned int rows, unsigned int columns, unsigned char indel_score) {
  unsigned int width = ((rows+14) & ~0x0F) + 16,
               height = columns + rows - 1;
  unsigned char *mem = malloc(width * height + 16),
                *matrix = (unsigned char*)(((unsigned long int)mem + 15) &
                  ~(unsigned long int)0x0F),
                *cell,
                score;
  unsigned int i,
               j;

  // Set the first column to 0.
  for (i = 0, cell = matrix; i < columns; i++, cell += width)
    *cell = 0;

  // Set the first row to 0, 1, 2, 3, 4, ... times the indel_score
  for (i = 0, cell = matrix, score = 0; i < rows; i++, cell += width + 1) {
    *cell = score;
    score = (score > 255 - indel_score)? 255 : score + indel_score;
    for (j = 1; j < (width - i); j++)
      *(cell + j) = 255;  // This protects the second row.
  }

  return mem;
}//_make_matrix

/*
Reverse a sequence.

:arg char *seq: Sequence to reverse.
:arg char *seqr: Pointer to a buffer that will receive the reversed sequence.
:arg size_t len: Length of the sequence, excluding the terminating NUL byte.
*/
void revseq(char *seq, char *seqr, size_t len){
  char *p = seq, *q = seqr + len;
  *(q--) = 0;
  while (q >= seqr)
    *(q--) = *(p++);
}

/*
Fill the alignment matrix.

:arg unsigned char *mem: The alignment matrix. The actual matrix is assumed to
  start at the first 16-byte aligned byte.
:arg unsigned int rows: Number of rows in the matrix.
:arg unsigned int columns: Number of columns in the matrix.
:arg char *seq1: The first sequence to be aligned.
:arg char *seq2: The second sequence to be aligned.
:arg unsigned char indel_score: Penalty score for insertions and deletions.
*/
void _align(
    unsigned char *mem, unsigned int rows, unsigned int columns,
    char *seq1, char *seq2, unsigned char indel_score) {
  unsigned int x = 1,
               y = 1,
               width = ((rows+14) & ~0x0F) + 16,
               end = columns + _min(15, rows - 1),
               limit;
  unsigned char *matrix = (unsigned char*)(((unsigned long int)mem + 15) &
                  ~(unsigned long int)0x0F),
                *d = matrix,
                *l = matrix + width,
                *i = l + width + 1;

  // Get copy of seq1 and reverse of seq2, making sure
  // that we can read 16 bytes (of garbage) past the end.
  const size_t seq2len = strlen(seq2);
  char *seq1f = malloc(strlen(seq1) + 16),
       *seq2r = malloc(seq2len + 16);
  strcpy(seq1f, seq1);
  revseq(seq2, seq2r, seq2len);

  const __m128i ones = _mm_set1_epi8(1),
    indel_scores = _mm_set1_epi8(indel_score),
    range = _mm_set_epi8(15, 14, 13, 12, 11, 10, 9, 8, 7 ,6, 5, 4, 3, 2, 1, 0);
  __m128i mi,
          md = _mm_load_si128((__m128i*)d),
          ml = _mm_load_si128((__m128i*)l),
          mu = _mm_loadu_si128((__m128i*)(l + 1)),
          mx = _mm_loadu_si128((__m128i*)seq1f),
          my = _mm_loadu_si128((__m128i*)(seq2r + seq2len - 1));

  while (1) {
    mi = _mm_min_epu8(_mm_adds_epu8(_mm_min_epu8(ml, mu), indel_scores),
      _mm_adds_epu8(md, _mm_add_epi8(_mm_cmpeq_epi8(mx, my), ones)));

    limit = y - x + 1;
    if (limit >= 16 || y >= (rows - 1))
      _mm_storeu_si128((__m128i*)i, mi);
    else
      // Need to make sure the top row and left column stay valid.
      // This is VERY MUCH SLOWER than the above.
      _mm_maskmoveu_si128(mi, _mm_cmplt_epi8(
        range, _mm_set1_epi8(limit)), (char*)i);

    if (++y < end) {
      // Move down one row.
      l += width;
      i += width;
      md = ml;
      mu = mi;
      ml = _mm_load_si128((__m128i*)l);

      // Move to the next base in seq2r.
      my = _mm_slli_si128(my, 1);
      if (limit - 1 < columns)
        my = _mm_insert_epi16(my, *(short*)(seq2r + seq2len - limit - 1), 0);
    }
    else if ((x += 16) < rows) {
      // Move right 16 columns.
      y = x;
      end += _min(16, rows - x);
      d += 16 * width + 16;
      l = d + width;
      i = l + width + 1;
      md = _mm_load_si128((__m128i*)d);
      ml = _mm_load_si128((__m128i*)l);
      mu = _mm_loadu_si128((__m128i*)(l + 1));
      mx = _mm_loadu_si128((__m128i*)(seq1f + x - 1));
      my = _mm_loadu_si128((__m128i*)(seq2r + seq2len - 1));
    }
    else break;
  }

  free(seq1f);
  free(seq2r);
}//_align

/*
Find the minimum distance, ignoring a trailing gap in the sequence associated
with the number of rows in an alignment matrix. If the minimum distance is
found, also return the row number.

It is assumed that the number of rows is larger than the number of columns.

:arg unsigned char *mem: An {rows} * {columns} matrix. The actual matrix is
  assumed to start at the first 16-byte aligned byte.
:arg unsigned int rows: Number of rows in the matrix.
:arg unsigned int columns: Number of columns in the matrix. MUST NOT be less
  than rows.

:returns alignment: The minimum distance and its row number.
:rtype: alignment
*/
alignment _find_min(
    unsigned char *mem, unsigned int rows, unsigned int columns) {
  unsigned int width = ((rows+14) & ~0x0F) + 16,
               i;
  unsigned char *matrix = (unsigned char*)(((unsigned long int)mem + 15) &
                  ~(unsigned long int)0x0F),
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
}//_find_min

/*
Do a semi-global alignment of {seq2} to {seq1}.

:arg char *seq1: The sequence to be aligned to.
:arg char *seq2: The sequence to be aligned.
:arg unsigned char indel_score: Penalty score for insertions and deletions (max
  255).

:returns alignment: The minimum distance and its row number.
*/
alignment align(char *seq1, char *seq2, char indel_score) {
  alignment a;
  unsigned char *matrix;
  unsigned int rows = strlen(seq1) + 1,
               columns = strlen(seq2) + 1;

  // The alignment is optimised for tall matrices.
  matrix = _make_matrix(columns, rows, indel_score);
  _align(matrix, columns, rows, seq2, seq1, indel_score);
  a = _find_min(matrix, columns, rows);
  free(matrix);
  return a;
}//align

#else
/******************************************************************************
  Plain Implementation
******************************************************************************/

/*
Initialise a matrix for semi-global alignment.

:arg char *matrix: The alignment matrix.
:arg int rows: Number of rows in the matrix.
:arg int columns: Number of columns in the matrix.
:arg char indel_score: Penalty score for insertions and deletions.
*/
void _init_matrix(char *matrix, int rows, int columns, char indel_score) {
  typedef char array_t[rows][columns];
  array_t *_matrix = (array_t *)matrix;
  int i;

  for (i = 1; i < rows; i++)
    (*_matrix)[i][0] = 0;

  for (i = 0; i < columns; i++)
    (*_matrix)[0][i] = i * indel_score;
}//_init_matrix

/*
Fill the alignment matrix.

:arg char *matrix: The alignment matrix.
:arg int rows: Number of rows in the matrix.
:arg int columns: Number of columns in the matrix.
:arg char *seq1: The sequence to be aligned to.
:arg char *seq2: The sequence to be aligned.
:arg char indel_score: Penalty score for insertions and deletions.
*/
void _align(
    char *matrix, int rows, int columns, char *seq1, char *seq2,
    char indel_score) {
  typedef char array_t[rows][columns];
  array_t *_matrix = (array_t *)matrix;
  int r,
      c;

  for (r = 1; r < rows; r++)
    for (c = 1; c < columns; c++)
      (*_matrix)[r][c] = _min(
        _min((*_matrix)[r - 1][c], (*_matrix)[r][c - 1]) + indel_score,
        (*_matrix)[r - 1][c - 1] + (seq1[r - 1] != seq2[c - 1]));
}//_align

/*
Find the minimum distance, ignoring a trailing gap in the sequence associated
with the number of rows in an alignment matrix. If the minimum distance is
found, also return the row number.

:arg char *matrix: An {rows} * {columns} matrix.
:arg int rows: Number of rows in the matrix.
:arg int columns: Number of columns in the matrix.

:returns alignment: The minimum distance and its row number.
*/
alignment _find_min(char *matrix, int rows, int columns) {
  typedef char array_t[rows][columns];
  array_t *_matrix = (array_t *)matrix;
  alignment a;
  int r;

  a.distance = columns - 1;
  a.position = 0;
  for (r = 1; r < rows; r++)
    if ((*_matrix)[r][columns - 1] < a.distance) {
      a.distance = (*_matrix)[r][columns - 1];
      a.position = r;
    }

  return a;
}//_find_min

/*
Do a semi-global alignment of {seq2} to {seq1}.

:arg char *seq1: The sequence to be aligned to.
:arg char *seq2: The sequence to be aligned.
:arg char indel_score: Penalty score for insertions and deletions.

:returns alignment: The minimum distance and its row number.
*/
alignment align(char *seq1, char *seq2, char indel_score) {
  alignment a;
  int rows = strlen(seq1) + 1,
             columns = strlen(seq2) + 1;
  char *matrix = malloc(rows * columns * sizeof(char));

  _init_matrix(matrix, rows, columns, indel_score);
  _align(matrix, rows, columns, seq1, seq2, indel_score);
  a = _find_min(matrix, rows, columns);
  free(matrix);

  return a;
}//align

#endif
