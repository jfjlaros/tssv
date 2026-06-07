#include <catch.hpp>

#include "../tssv/sgAlign.h"

/*
 * Private function prototypes.
 */
void initMatrix_(int *const, size_t const, size_t const, int const);
void align_(
  int *const, size_t const, size_t const,
  char const *const, char const *const, int const);
Alignment findMin_(int const *const, size_t const, size_t const);


TEST_CASE("Internal alignment function", "[matrix]") {
  int *matrix = (int *)malloc(1000); //
}

TEST_CASE("Initialise matrix", "[matrix]") {
  int *const matrix = (int *const)malloc(100 * sizeof(int));
  int (*const matrix_)[10] = (int (*const)[10])matrix;

  initMatrix_(matrix, 10, 10, 1);
  REQUIRE(matrix_[0][0] == 0);
  REQUIRE(matrix_[1][0] == 0);
  REQUIRE(matrix_[9][0] == 0);
  REQUIRE(matrix_[0][1] == 1);
  REQUIRE(matrix_[0][9] == 9);

  initMatrix_(matrix, 10, 10, 3);
  REQUIRE(matrix_[0][0] == 0);
  REQUIRE(matrix_[1][0] == 0);
  REQUIRE(matrix_[9][0] == 0);
  REQUIRE(matrix_[0][1] == 3);
  REQUIRE(matrix_[0][9] == 27);
}

TEST_CASE("Alignment", "[matrix]") {
  Alignment a = align(
      "GCCAACTGTTACCAAGGTCCCTCCCATGCATGCTGCTCCCTACAGAGGCATGTGCACAGT",
      "CTGTTTCCAAGG", 1);

  REQUIRE(a.distance == 1);
}

