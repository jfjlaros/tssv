#!/usr/bin/env python

"""
Library with functions for a semi-global alignment.
"""

def __makeMatrix(xSize, ySize):
    """
    Initialise a matrix for semi-global alignment.

    @arg xSize: Size of the x dimension of the matrix.
    @type xSize: int
    @arg xSize: Size of the y dimension of the matrix.
    @type ySize: int

    @returns matrix: The alignment matrix.
    @rtype: list[][]
    """
    return [range(ySize)] + [[0] * ySize for i in range(xSize - 1)]
#__makeMatrix

def __align(matrix, xSize, ySize, seq1, seq2):
    """
    Fill the alignment matrix.

    @arg xSize: Size of the x dimension of the matrix.
    @type xSize: int
    @arg xSize: Size of the y dimension of the matrix.
    @type ySize: int
    @arg seq1: The sequence to be aligned to.
    @type seq1: string
    @arg seq2: The sequence to be aligned.
    @type seq2: string
    """
    for x in range(1, xSize):
        for y in range(1, ySize):
            matrix[x][y] = min(matrix[x -1][y] + 1, matrix[x][y - 1] + 1,
                matrix[x - 1][y - 1] + int(seq1[x - 1] != seq2[y - 1]))
#__align

def __findMin(matrix, xSize, ySize):
    """
    Find the minimum distance, ignoring a trailing gap in the sequence
    associated with the number of rows in an alignment matrix. If the minimum
    distance is found, also return the row number.

    It is assumed that the number of rows is larger than the number of columns.

    @arg matrix: An {xSize} * {ySize} matrix.
    @arg xSize: Size of the x dimension of the matrix.
    @type xSize: int
    @arg xSize: Size of the y dimension of the matrix.
    @type ySize: int

    @returns: The minimum distance and its row number.
    @rtype minimum: tuple(int, int)
    """
    minimum = ySize - 1
    xMin = 0

    for x in range(xSize):
        if (matrix[x][ySize - 1] < minimum):
            minimum = matrix[x][ySize -1]
            xMin = x
        #if
    #for

    return minimum, xMin
#__findMin

def align(seq1, seq2):
    """
    Do a semi-global alignment of {seq2} to {seq1}.

    @arg seq1: The sequence to be aligned to.
    @type seq1: string
    @arg seq2: The sequence to be aligned.
    @type seq2: string

    @returns: The minimum distance and its row number.
    @rtype minimum: tuple(int, int)
    """
    xSize = len(seq1) + 1
    ySize = len(seq2) + 1

    matrix = __makeMatrix(xSize, ySize)
    __align(matrix, xSize, ySize, seq1, seq2)

    return __findMin(matrix, xSize, ySize)
#align
