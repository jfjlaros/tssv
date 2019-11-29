#include <catch.hpp>

#include "../tssv/sgAlign.h"

/*
 * Private function prototypes.
 */
void _initMatrix(int *, int, int, int);
void _align(int *, int, int, char *, char *, int);
alignment _findMin(int *, int, int);


TEST_CASE("Internal alignment function", "[matrix]") {
  int *matrix = (int *)malloc(1000); //
}

TEST_CASE("Initialise matrix", "[matrix]") {
  int *matrix = (int *)malloc(100 * sizeof(int)),
      (*_matrix)[10] = (int (*)[10])matrix;

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

