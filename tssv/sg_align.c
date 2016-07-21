/*
Library with functions for a semi-global alignment.
*/

#include <string.h>
#include <stdlib.h>

#include "sg_align.h"

/*
Calculate the minimum of two values.

:arg a: A value.
:type a: int
:arg b: A value.
:type c: int

:returns: The minimum of {a} and {b}.
:rtype: int
*/
int _min(int a, int b) {
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
:arg x_size: Size of the y dimension of the matrix.
:type y_size: int
:arg seq1: The sequence to be aligned to.
:type seq1: char *
:arg seq2: The sequence to be aligned.
:type seq2: char *
*/
void _align(int *matrix, int x_size, int y_size, char *seq1, char *seq2) {
  typedef int array_t[x_size][y_size];
  array_t *_matrix = (array_t *)matrix;
  int x,
      y;

  for (x = 1; x < x_size; x++)
    for (y = 1; y < y_size; y++)
      (*_matrix)[x][y] = _min(
        _min((*_matrix)[x - 1][y] + 1, (*_matrix)[x][y - 1] + 1),
        (*_matrix)[x - 1][y - 1] + (int)(seq1[x - 1] != seq2[y - 1]));
}//_align

/*
Find the minimum distance, ignoring a trailing gap in the sequence associated
with the number of rows in an alignment matrix. If the minimum distance is
found, also return the row number.

It is assumed that the number of rows is larger than the number of columns.

:arg matrix: An {x_size} * {y_size} matrix.
:type matrix: int *
:arg x_size: Size of the x dimension of the matrix.
:type x_size: int
:arg x_size: Size of the y dimension of the matrix.
:type y_size: int

:returns: The minimum distance and its row number.
:rtype: alignment
*/
alignment _find_min(int *matrix, int x_size, int y_size) {
  typedef int array_t[x_size][y_size];
  array_t *_matrix = (array_t *)matrix;
  alignment a;
  int x;

  a.distance = y_size - 1;
  a.position = 0;
  for (x = 1; x < x_size; x++)
    if ((*_matrix)[x][y_size - 1] < a.distance) {
      a.distance = (*_matrix)[x][y_size - 1];
      a.position = x;
    }//if

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
  int x_size = strlen(seq1) + 1,
      y_size = strlen(seq2) + 1,
      *matrix = malloc(x_size * y_size * sizeof(int));

  _align(matrix, x_size, y_size, seq1, seq2);
  a = _find_min(matrix, x_size, y_size);
  free(matrix);

  return a;
}//align
