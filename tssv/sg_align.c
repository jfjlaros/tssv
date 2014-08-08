/*
  Library with functions for a semi-global alignment.
*/

#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "sg_align.h"

/*
  Initialise a matrix for semi-global alignment.
  
  :arg x_size: Size of the x dimension of the matrix.
  :type x_size: int
  :arg x_size: Size of the y dimension of the matrix.
  :type y_size: int
  
  :returns matrix: The alignment matrix.
  :rtype: int **
*/
int **_make_matrix(int x_size, int y_size) {
  int **matrix = malloc(x_size * sizeof(int *)),
      i;

  for (i = 0; i < x_size - 1; i++)
    matrix[i] = malloc(y_size * sizeof(int));
}//_make_matrix

/*
  Free a matrix for semi-global alignment.

  :arg matrix: The alignment matrix.
  :type matrix: int **
  :arg x_size: Size of the x dimension of the matrix.
  :type x_size: int
*/
void _free_matrix(int **matrix, int x_size) {
  int i;

  for (i = 0; i < x_size - 1; i++)
    free(matrix[i]);
  free(matrix);
}//_free_matrix

/*
  Fill the alignment matrix.

  :arg matrix: The alignment matrix.
  :type matrix: int **
  :arg x_size: Size of the x dimension of the matrix.
  :type x_size: int
  :arg x_size: Size of the y dimension of the matrix.
  :type y_size: int
  :arg seq1: The sequence to be aligned to.
  :type seq1: char *
  :arg seq2: The sequence to be aligned.
  :type seq2: char *
*/
void _align(int **matrix, int x_size, int y_size, char *seq1, char *seq2) {
  int x,
      y;

  for (x = 1; x < x_size; x++)
    for (y = 1; y < y_size; y++)
      matrix[x][y] = min(min(matrix[x - 1][y] + 1, matrix[x][y - 1] + 1),
        matrix[x - 1][y - 1] + (int)(seq1[x - 1] != seq2[y - 1]));
}//_align

/*
  Find the minimum distance, ignoring a trailing gap in the sequence
  associated with the number of rows in an alignment matrix. If the minimum
  distance is found, also return the row number.

  It is assumed that the number of rows is larger than the number of columns.

  :arg matrix: An {x_size} * {y_size} matrix.
  :type matrix: int **
  :arg x_size: Size of the x dimension of the matrix.
  :type x_size: int
  :arg x_size: Size of the y dimension of the matrix.
  :type y_size: int

  :returns: The minimum distance and its row number.
  :rtype: int *
*/
int *_find_min(int **matrix, int x_size, int y_size) {
  int *minimum = malloc(2 * sizeof(int)),
      x;

  minimum[0] = y_size - 1;
  minimum[1] = 0;
  for (x = 0; x < x_size; x++)
    if (matrix[x][y_size - 1] < minimum[0]) {
      minimum[0] = matrix[x][y_size - 1];
      minimum[1] = x;
    }//if

    return minimum;
}//_find_min

/*
  Do a semi-global alignment of {seq2} to {seq1}.

  :arg seq1: The sequence to be aligned to.
  :type seq1: char *
  :arg seq2: The sequence to be aligned.
  :type seq2: char *

  :returns: The minimum distance and its row number.
  :rtype: int *
*/
int *align(char *seq1, char *seq2) {
    int **matrix,
        *minimum,
        x_size = strlen(seq1) + 1,
        y_size = strlen(seq2) + 1;

    matrix = _make_matrix(x_size, y_size);
    _align(matrix, x_size, y_size, seq1, seq2);
    minimum = _find_min(matrix, x_size, y_size);
    _free_matrix(matrix, x_size);

    return minimum;
}//align
