#include <catch.hpp>

#include "../tssv/sgAlign.h"

/*
 * Private function prototypes.
 */
void _initMatrix(unsigned char *, unsigned int, unsigned int, unsigned char);
void _align(
  unsigned char *, unsigned int, unsigned int, char *, char *, unsigned char);
alignment _findMin(unsigned char *, unsigned int, unsigned int);


TEST_CASE("Internal alignment function", "[matrix]") {
  unsigned char *matrix = (unsigned char *)malloc(1000); //
}

TEST_CASE("Initialise matrix", "[matrix]") {
  unsigned char *matrix = (unsigned char *)malloc(100 * sizeof(unsigned char)),
                (*_matrix)[10] = (unsigned char (*)[10])matrix;

  _initMatrix(matrix, 10, 10, 1);
  REQUIRE(_matrix[0][0] == 0);
  REQUIRE(_matrix[1][0] == 0);
  REQUIRE(_matrix[9][0] == 0);
  REQUIRE(_matrix[0][1] == 1);
  REQUIRE(_matrix[0][9] == 9);

  _initMatrix(matrix, 10, 10, 3);
  REQUIRE(_matrix[0][0] == 0);
  REQUIRE(_matrix[1][0] == 0);
  REQUIRE(_matrix[9][0] == 0);
  REQUIRE(_matrix[0][1] == 3);
  REQUIRE(_matrix[0][9] == 27);
}

TEST_CASE("Alignment", "[matrix]") {
  alignment a = align(
      (char *)"GCCAACTGTTACCAAGGTCCCTCCCATGCATGCTGCTCCCTACAGAGGCATGTGCACAGT",
      (char *)"CTGTTTCCAAGG",
      1);

  REQUIRE(a.distance == 1);
}

