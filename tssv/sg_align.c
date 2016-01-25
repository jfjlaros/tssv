/*
Library with functions for a semi-global alignment.


Two implementations are contained within this file:
1. An implementation using SSE2 instructions, which is used on systems that
   support these instructions. The SSE2 implementation uses a matrix of unsigned
   chars instead of ints and stores the matrix diagonally. The implementation is
   optimised for tall matrices; it will transpose the matrix if it is wide.
   This implementation achieves a 3-4x performance improvement compaired to
   TSSV 0.2.5, on data consiting predominantly of reads of 200-300bp.
2. A fallback implementation which is used if SSE2 is not available. Compared to
   TSSV 0.2.5, this optimised implementation achieves about a 1.5x performance
   improvement, most of which is attributable to allocating the alignment matrix
   as a single block of memory and using pointers instead of array indices to
   access it.
*/

#include <string.h>
#include <stdlib.h>
#include "sg_align.h"

#if defined(__GNUC__) && defined(__SSE2__)
#  include <xmmintrin.h>
#  include <emmintrin.h>
/******************************************************************************
  SSE2-enabled Implementation
******************************************************************************/


/*
Calculate the minimum of two values.

:arg a: A value.
:type a: char
:arg b: A value.
:type c: char

:returns: The minimum of {a} and {b}.
:rtype: char
*/
inline unsigned char _min(unsigned char a, unsigned char b) {
  if (a < b)
    return a;
  return b;
}//_min


/*
Initialise a matrix for semi-global alignment.

:arg x_size: Size of the x dimension of the matrix.
:type x_size: int
:arg y_size: Size of the y dimension of the matrix. MUST NOT be less
than x_size.
:type y_size: int
:arg transposed: Whether the matrix is transposed.
:type transposed: int

:returns: The alignment matrix. The actual matrix starts at the first
16-byte aligned byte and is stored diagonally.
:rtype: char *
*/
unsigned char *_make_matrix(unsigned int x_size, unsigned int y_size, int transposed) {
  unsigned int width = ((x_size+14) & ~0x0F) + 16,
    height = y_size + x_size - 1;
  unsigned char *mem = malloc(width * height + 16),
    *matrix = (unsigned char*)(((unsigned long int)mem + 15) & ~(unsigned long int)0x0F),
    *cell;
  unsigned int i, j;

  if (transposed) {
      // Set the first column to 0.
      for (i = 0, cell = matrix; i < y_size; i++, cell += width)
        *cell = 0;

      // Set the first row to 0, 1, 2, 3, 4, ...
      for (i = 0, cell = matrix; i < x_size; i++, cell += width + 1) {
        *cell = (i > 255? 255 : i);
        for (j = 1; j < (width - i); j++)
            *(cell + j) = 255;  // This protects the second row.
      }
  }
  else {
      // Set the first column to 0, 1, 2, 3, 4, ...
      for (i = 0, cell = matrix; i < y_size; i++, cell += width)
        *cell = (i > 255? 255 : i);

      // Set the first row to 0.
      for (i = 0, cell = matrix; i < x_size; i++, cell += width + 1) {
        *cell = 0;
        for (j = 1; j < (width - i); j++)
            *(cell + j) = 255;  // This protects the second row.
      }
  }

  return mem;
}//_make_matrix

/*
Reverse a sequence.

:arg seq: Sequence to reverse.
:type seq: char *
:arg seqr: Pointer to a buffer that will receive the reversed sequence.
:type seqr: char *
:arg len: Length of the sequence, excluding the terminating NUL byte.
:type len: size_t
*/
void revseq(char* seq, char* seqr, size_t len){
  char *p = seq, *q = seqr + len;
  *(q--) = 0;
  while (q >= seqr)
    *(q--) = *(p++);
}

/*
Fill the alignment matrix.

:arg mem: The alignment matrix. The actual matrix is assumed to start at
the first 16-byte aligned byte.
:type mem: char *
:arg x_size: Size of the x dimension of the matrix.
:type x_size: int
:arg y_size: Size of the y dimension of the matrix.
:type y_size: int
:arg seq1: The first sequence to be aligned.
:type seq1: char *
:arg seq2: The second sequence to be aligned.
:type seq2: char *
*/
void _align(unsigned char *mem, unsigned int x_size, unsigned int y_size, char *seq1, char *seq2) {
  unsigned int x = 1, y = 1,
    width = ((x_size+14) & ~0x0F) + 16,
    end = y_size + _min(15, x_size - 1),
    limit;
  unsigned char *matrix = (unsigned char*)(((unsigned long int)mem + 15) & ~(unsigned long int)0x0F),
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
    range = _mm_set_epi8(15, 14, 13, 12, 11, 10, 9, 8, 7 ,6, 5, 4, 3, 2, 1, 0);
  __m128i mi,
          md = _mm_load_si128((__m128i*)d),
          ml = _mm_load_si128((__m128i*)l),
          mu = _mm_loadu_si128((__m128i*)(l + 1)),
          mx = _mm_loadu_si128((__m128i*)seq1f),
          my = _mm_loadu_si128((__m128i*)(seq2r + seq2len - 1));

  while (1) {
    mi = _mm_min_epu8(_mm_adds_epu8(_mm_min_epu8(ml, mu), ones),
      _mm_adds_epu8(md, _mm_add_epi8(_mm_cmpeq_epi8(mx, my), ones)));

    limit = y - x + 1;
    if (limit >= 16 || y >= (x_size - 1))
      _mm_storeu_si128((__m128i*)i, mi);
    else
      // Need to make sure the top row and left column stay valid.
      // This is VERY MUCH SLOWER than the above.
      _mm_maskmoveu_si128(mi, _mm_cmplt_epi8(range, _mm_set1_epi8(limit)), (char*)i);

    if (++y < end) {
      // Move down one row.
      l += width;
      i += width;
      md = ml;
      mu = mi;
      ml = _mm_load_si128((__m128i*)l);

      // Move to the next base in seq2r.
      my = _mm_slli_si128(my, 1);
      if (limit - 1 < y_size)
        my = _mm_insert_epi16(my, *(short*)(seq2r + seq2len - limit - 1), 0);
    }
    else if ((x += 16) < x_size) {
      // Move right 16 columns.
      y = x;
      end += _min(16, x_size - x);
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

:arg mem: An {x_size} * {y_size} matrix. The actual matrix is assumed
to start at the first 16-byte aligned byte.
:type mem: char *
:arg x_size: Size of the x dimension of the matrix.
:type x_size: int
:arg y_size: Size of the y dimension of the matrix. MUST NOT be less
than x_size.
:type y_size: int
:arg transposed: Whether the matrix is transposed.
:type transposed: int

:returns: The minimum distance and its row number.
:rtype: alignment
*/
alignment _find_min(unsigned char *mem, unsigned int x_size, unsigned int y_size, int transposed) {
  unsigned int width = ((x_size+14) & ~0x0F) + 16,
    i;
  unsigned char *matrix = (unsigned char*)(((unsigned long int)mem + 15) & ~(unsigned long int)0x0F),
    *cell;
  alignment a;

  a.position = 0;
  if (transposed) {
    a.distance = x_size - 1;
    cell = matrix + (x_size - 1) * width + x_size - 1;
    for (i = 0; i < y_size; i++, cell += width)
      if (*cell < a.distance) {
        a.distance = *cell;
        a.position = i;
      }
  }
  else {
    a.distance = y_size - 1;
    cell = matrix + (y_size - 1) * width;
    width++; // Saves us a +1 in the loop below...
    for (i = 0; i < x_size; i++, cell += width)
      if (*cell < a.distance) {
        a.distance = *cell;
        a.position = i;
      }
  }

  return a;
}//_find_min

/*
Do a semi-global alignment of {seq2} to {seq1}.

:arg seq1: The sequence to be aligned to.
:type seq1: char *
:arg seq2: The sequence to be aligned.
:type seq2: char *

:returns: The minimum distance and its row number.
:rtype: alignment
*/
alignment align(char *seq1, char *seq2) {
  alignment a;
  unsigned char *matrix;
  unsigned int x_size = strlen(seq1) + 1,
    y_size = strlen(seq2) + 1;

  if (x_size > y_size) {
      // The alignment is optimised for tall matrices, transpose it.
      matrix = _make_matrix(y_size, x_size, 1);
      _align(matrix, y_size, x_size, seq2, seq1);
      a = _find_min(matrix, y_size, x_size, 1);
      free(matrix);
  }
  else {
      matrix = _make_matrix(x_size, y_size, 0);
      _align(matrix, x_size, y_size, seq1, seq2);
      a = _find_min(matrix, x_size, y_size, 0);
      free(matrix);
  }
  return a;
}//align

#else
/******************************************************************************
  Plain Implementation
******************************************************************************/

/*
Initialise a matrix for semi-global alignment.

:arg x_size: Size of the x dimension of the matrix.
:type x_size: int
:arg x_size: Size of the y dimension of the matrix.
:type y_size: int

:returns: The alignment matrix.
:rtype: int *
*/
unsigned int *_make_matrix(unsigned int x_size, unsigned int y_size) {
  unsigned int i,
    *matrix = malloc(x_size * y_size * sizeof(int)),
    *cell;

  for (i = 0, cell = matrix; i < y_size; i++, cell++)
    *cell = i;

  for (i = 0, cell = matrix; i < x_size; i++, cell += y_size)
    *cell = 0;

  return matrix;
}//_make_matrix

/*
Calculate the minimum of two values.

:arg a: A value.
:type a: int
:arg b: A value.
:type c: int

:returns: The minimum of {a} and {b}.
:rtype: int
*/
inline unsigned int _min(unsigned int a, unsigned int b) {
  if (a < b)
    return a;
  return b;
}//_min

/*
Fill the alignment matrix.

:arg matrix: The alignment matrix.
:type matrix: int *
:arg x_size: Size of the x dimension of the matrix.
:type x_size: int
:arg y_size: Size of the y dimension of the matrix.
:type y_size: int
:arg seq1: The sequence to be aligned to.
:type seq1: char *
:arg seq2: The sequence to be aligned.
:type seq2: char *
*/
void _align(unsigned int *matrix, unsigned int x_size, unsigned int y_size, char *seq1, char *seq2) {
  unsigned int x, y,
      *d = matrix,
      *l = d + 1,
      *u = d + y_size,
      *i = l + y_size;

  for (x = 1; x < x_size; x++, d++, l++, u++, i++)
    for (y = 1; y < y_size; y++, d++, l++, u++, i++)
      *i = _min(_min(*l, *u) + 1,
        *d + (unsigned int)(seq1[x - 1] != seq2[y - 1]));
}//_align

/*
Find the minimum distance, ignoring a trailing gap in the sequence associated
with the number of rows in an alignment matrix. If the minimum distance is
found, also return the row number.

:arg matrix: An {x_size} * {y_size} matrix.
:type matrix: int *
:arg x_size: Size of the x dimension of the matrix.
:type x_size: int
:arg x_size: Size of the y dimension of the matrix.
:type y_size: int

:returns: The minimum distance and its row number.
:rtype: alignment
*/
alignment _find_min(unsigned int *matrix, unsigned int x_size, unsigned int y_size) {
  alignment a;
  unsigned int x, *cell;

  a.distance = y_size - 1;
  a.position = 0;
  for (x = 1, cell = matrix + 2 * y_size - 1; x < x_size; x++, cell += y_size)
    if (*cell < a.distance) {
      a.distance = *cell;
      a.position = x;
    }

  return a;
}//_find_min

/*
Do a semi-global alignment of {seq2} to {seq1}.

:arg seq1: The sequence to be aligned to.
:type seq1: char *
:arg seq2: The sequence to be aligned.
:type seq2: char *

:returns: The minimum distance and its row number.
:rtype: alignment
*/
alignment align(char *seq1, char *seq2) {
  alignment a;
  unsigned int *matrix,
      x_size = strlen(seq1) + 1,
      y_size = strlen(seq2) + 1;

  matrix = _make_matrix(x_size, y_size);
  _align(matrix, x_size, y_size, seq1, seq2);
  a = _find_min(matrix, x_size, y_size);
  free(matrix);

  return a;
}//align

#endif