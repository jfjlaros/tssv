#pragma once

typedef struct {
  size_t distance;
  size_t position;
} Alignment;


/*! Do a semi-global alignment of `seq2` to `seq1`.
 *
 * \param [in] seq1 The sequence to be aligned to.
 * \param [in] seq2 The sequence to be aligned.
 * \param [in] indelScore Penalty score for insertions and deletions.
 *
 * \return The minimum distance and its row number.
 */
Alignment align(
  char const *const seq1, char const *const seq2, int const indelScore);
