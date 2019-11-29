/*
 * Library for semi-global alignment.
 */
#include <string.h>
#include <stdlib.h>

#include "sgAlign.h"


/**
 * Calculate the minimum of two values.
 *
 * @arg {int} a - A value.
 * @arg {int} b - A value.
 *
 * @return {int} - The minimum of {a} and {b}.
 */
static inline int _min(int a, int b) {
  if (a < b)
    return a;
  return b;
}

/**
 * Initialise a matrix for semi-global alignment.
 *
 * @arg {int *} matrix - The alignment matrix.
 * @arg {int} rows - Number of rows in the matrix.
 * @arg {int} columns - Number of columns in the matrix.
 * @arg {int} indelScore - Penalty score for insertions and
 *   deletions.
 */
void _initMatrix(int *matrix, int rows, int columns, int indelScore) {
  int (*_matrix)[columns] = (int (*)[columns])matrix;
  int i;

  for (i = 1; i < rows; i++)
    _matrix[i][0] = 0;

  for (i = 0; i < columns; i++)
    _matrix[0][i] = i * indelScore;
}

/**
 * Fill the alignment matrix.
 *
 * @arg {int *} matrix - The alignment matrix.
 * @arg {int} rows - Number of rows in the matrix.
 * @arg {int} columns - Number of columns in the matrix.
 * @arg {char *} seq1 - The sequence to be aligned to.
 * @arg {char *} seq2 - The sequence to be aligned.
 * @arg {int} indelScore - Penalty score for insertions and
 *   deletions.
 */
void _align(
    int *matrix, int rows, int columns, char *seq1, char *seq2,
    int indelScore) {
  int (*_matrix)[columns] = (int (*)[columns])matrix;
  int r,
               c;

  for (r = 1; r < rows; r++)
    for (c = 1; c < columns; c++)
      _matrix[r][c] = _min(
        _min(_matrix[r - 1][c], _matrix[r][c - 1]) + indelScore,
        _matrix[r - 1][c - 1] + (seq1[r - 1] != seq2[c - 1]));
}

/**
 * Find the minimum distance, ignoring a trailing gap in the sequence
 * associated with the number of rows in an alignment matrix. If the minimum
 * distance is found, also return the row number.
 *
 * @arg {int *} matrix - An {rows} * {columns} matrix.
 * @arg {int} rows - Number of rows in the matrix.
 * @arg {int} columns - Number of columns in the matrix.
 *
 * @return {alignment} - The minimum distance and its row number.
 */
alignment _findMin(int *matrix, int rows, int columns) {
  int (*_matrix)[columns] = (int (*)[columns])matrix;
  alignment a;
  int r;

  a.distance = columns - 1;
  a.position = 0;
  for (r = 1; r < rows; r++)
    if (_matrix[r][columns - 1] < a.distance) {
      a.distance = _matrix[r][columns - 1];
      a.position = r;
    }

  return a;
}

/**
 * Do a semi-global alignment of {seq2} to {seq1}.
 *
 * @arg {char *} seq1 - The sequence to be aligned to.
 * @arg {char *} seq2 - The sequence to be aligned.
 * @arg {int} indelScore - Penalty score for insertions and
 *   deletions.
 *
 * @return {alignment} - The minimum distance and its row number.
 */
alignment align(char *seq1, char *seq2, int indelScore) {
  alignment a;
  int rows = strlen(seq1) + 1,
               columns = strlen(seq2) + 1;
  int *matrix = (int *)malloc(
    rows * columns * sizeof(int));

  _initMatrix(matrix, rows, columns, indelScore);
  _align(matrix, rows, columns, seq1, seq2, indelScore);
  a = _findMin(matrix, rows, columns);
  free(matrix);

  return a;
}
