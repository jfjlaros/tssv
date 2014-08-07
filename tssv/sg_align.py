#!/usr/bin/env python

"""
Library with functions for a semi-global alignment.
"""

def _make_matrix(x_size, y_size):
    """
    Initialise a matrix for semi-global alignment.

    :arg x_size: Size of the x dimension of the matrix.
    :type x_size: int
    :arg x_size: Size of the y dimension of the matrix.
    :type y_size: int

    :returns matrix: The alignment matrix.
    :rtype: list[][]
    """
    return [range(y_size)] + [[0] * y_size for i in range(x_size - 1)]
#_make_matrix

def _align(matrix, x_size, y_size, seq1, seq2):
    """
    Fill the alignment matrix.

    :arg x_size: Size of the x dimension of the matrix.
    :type x_size: int
    :arg x_size: Size of the y dimension of the matrix.
    :type y_size: int
    :arg seq1: The sequence to be aligned to.
    :type seq1: string
    :arg seq2: The sequence to be aligned.
    :type seq2: string
    """
    for x in range(1, x_size):
        for y in range(1, y_size):
            matrix[x][y] = min(matrix[x -1][y] + 1, matrix[x][y - 1] + 1,
                matrix[x - 1][y - 1] + int(seq1[x - 1] != seq2[y - 1]))
#_align

def _find_min(matrix, x_size, y_size):
    """
    Find the minimum distance, ignoring a trailing gap in the sequence
    associated with the number of rows in an alignment matrix. If the minimum
    distance is found, also return the row number.

    It is assumed that the number of rows is larger than the number of columns.

    :arg matrix: An {x_size} * {y_size} matrix.
    :arg x_size: Size of the x dimension of the matrix.
    :type x_size: int
    :arg x_size: Size of the y dimension of the matrix.
    :type y_size: int

    :returns: The minimum distance and its row number.
    :rtype minimum: tuple(int, int)
    """
    minimum = y_size - 1
    x_min = 0

    for x in range(x_size):
        if (matrix[x][y_size - 1] < minimum):
            minimum = matrix[x][y_size -1]
            x_min = x
        #if
    #for

    return minimum, x_min
#_find_min

def align(seq1, seq2):
    """
    Do a semi-global alignment of {seq2} to {seq1}.

    :arg seq1: The sequence to be aligned to.
    :type seq1: string
    :arg seq2: The sequence to be aligned.
    :type seq2: string

    :returns: The minimum distance and its row number.
    :rtype minimum: tuple(int, int)
    """
    x_size = len(seq1) + 1
    y_size = len(seq2) + 1

    matrix = _make_matrix(x_size, y_size)
    _align(matrix, x_size, y_size, seq1, seq2)

    return _find_min(matrix, x_size, y_size)
#align
