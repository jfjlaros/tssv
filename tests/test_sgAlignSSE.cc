#include <catch.hpp>

#include "../tssv/sgAlignSSE.h"

/*
 * Private function prototypes.
 */
void _revSeq(char *, char *, size_t);
unsigned char *_makeMatrixSSE(unsigned int, unsigned int, unsigned char);
void _alignSSE(
  unsigned char *, unsigned int, unsigned int, char *, char *, unsigned char);
alignment _findMinSSE(unsigned char *, unsigned int, unsigned int);


TEST_CASE("Reverse a sequence", "[helper]") {
  char *seq1 = (char *)"ACGT",
       *seq2 = (char *)malloc(5 * sizeof(char));

  _revSeq(seq1, seq2, 4);
  REQUIRE(seq2[0] == 'T');
  REQUIRE(seq2[1] == 'G');
  REQUIRE(seq2[2] == 'C');
  REQUIRE(seq2[3] == 'A');
}

TEST_CASE("Create matrix", "[matrix]") {
  unsigned char (*matrix)[10] = (unsigned char (*)[10])_makeMatrixSSE(
    10, 10, 1);

  REQUIRE(matrix[0][0] == 0);
  // TODO: Test more values.
  free(matrix);

  // TODO: Add tests for indelSize > 1.
}
