/*
 * Library for semi-global alignment.
 */
#include <string.h>
#include <stdlib.h>

#include "sgAlign.h"


/*! Calculate the minimum of two values.
 *
 * \param [in] a A value.
 * \param [in] b A value.
 *
 * \return The minimum of `a` and `b`.
 */
static inline int min_(int const a, int const b) {
  if (a < b) {
    return a;
  }
  return b;
}

/*! Initialise a matrix for semi-global alignment.
 *
 * \param [in,out] matrix The alignment matrix.
 * \param [in] rows Number of rows in the matrix.
 * \param [in] columns Number of columns in the matrix.
 * \param [in] indelScore Penalty score for insertions and deletions.
 */
void initMatrix_(
    int *const matrix, size_t const rows, size_t const columns,
    int const indelScore) {
  int (*const matrix_)[columns] = (int (*const)[columns])matrix;

  for (size_t i = 1; i < rows; i++) {
    matrix_[i][0] = 0;
  }
  for (size_t i = 0; i < columns; i++) {
    matrix_[0][i] = i * indelScore;
  }
}

/*! Fill the alignment matrix.
 *
 * \param [in, out] matrix The alignment matrix.
 * \param [in] rows Number of rows in the matrix.
 * \param [in] columns Number of columns in the matrix.
 * \param [in] seq1 The sequence to be aligned to.
 * \param [in] seq2 The sequence to be aligned.
 * \param [in] indelScore Penalty score for insertions and deletions.
 */
void align_(
    int *const matrix, size_t const rows, size_t const columns,
    char const *const seq1, char const *const seq2, int const indelScore) {
  int (*const matrix_)[columns] = (int (*const)[columns])matrix;

  for (size_t r = 1; r < rows; r++) {
    for (size_t c = 1; c < columns; c++) {
      matrix_[r][c] = min_(
        min_(matrix_[r - 1][c], matrix_[r][c - 1]) + indelScore,
        matrix_[r - 1][c - 1] + (seq1[r - 1] != seq2[c - 1]));
    }
  }
}

/*! Find the minimum distance, ignoring a trailing gap in the sequence
 * associated with the number of rows in an alignment matrix. If the minimum
 * distance is found, also return the row number.
 *
 * \param [in] matrix A `rows` * `columns` matrix.
 * \param [in] rows Number of rows in the matrix.
 * \param [in] columns Number of columns in the matrix.
 *
 * \return The minimum distance and its row number.
 */
Alignment findMin_(
    int const *const matrix, size_t const rows, size_t const columns) {
  int const (*const matrix_)[columns] = (int const (*const)[columns])matrix;
  Alignment a = {columns - 1, 0};

  for (size_t r = 1; r < rows; r++) {
    if (matrix_[r][columns - 1] < a.distance) {
      a.distance = matrix_[r][columns - 1];
      a.position = r;
    }
  }
  return a;
}

Alignment align(
    char const *const seq1, char const *const seq2, int const indelScore) {
  Alignment a;
  size_t rows = strlen(seq1) + 1;
  size_t columns = strlen(seq2) + 1;
  int *matrix = (int *)malloc(rows * columns * sizeof(int));

  initMatrix_(matrix, rows, columns, indelScore);
  align_(matrix, rows, columns, seq1, seq2, indelScore);
  a = findMin_(matrix, rows, columns);
  free(matrix);

  return a;
}
