#include <Python.h>

#include "sgAlign.h"
#include "sgAlignSSE.h"


/**
 * Converter for alignment struct.
 */
PyObject *pyAlignment(int distance, int position) {
  return Py_BuildValue(
    "{s: i, s: i}", "distance", distance, "position", position);
}

/**
 * Wrapper for align function.
 */
PyObject *pyAlign(PyObject *self, PyObject *args) {
  char *seq1,
       *seq2;
  int indel_score;
  alignment a;

  if (!PyArg_ParseTuple(args, "ssi", &seq1, &seq2, &indel_score)) {
    return NULL;
  }

  a = align(seq1, seq2, indel_score);

  return pyAlignment(a.distance, a.position);
}

/**
 * Wrapper for alignSSE function.
 */
PyObject *pyAlignSSE(PyObject *self, PyObject *args) {
  char *seq1,
       *seq2,
       indel_score;
  alignment a;

  if (!PyArg_ParseTuple(args, "ssb", &seq1, &seq2, &indel_score)) {
    return NULL;
  }

  a = alignSSE(seq1, seq2, indel_score);

  return pyAlignment(a.distance, a.position);
}

/*
 * Module methods.
 */
PyMethodDef pySgAlignMethods[] = {
  {
    "align", pyAlign, METH_VARARGS,
    "Do a semi-global alignment of {seq2} to {seq1}.\n\n"
    "  :arg str seq1: The sequence to be aligned to.\n"
    "  :arg str seq2: The sequence to be aligned.\n"
    "  :arg int indel_score: Penalty score for insertions and deletions.\n\n"
    "  :returns dict alignment: The minimum distance and its row number.\n"},
  {
    "align_sse", pyAlignSSE, METH_VARARGS,
    "SSE implementation of {align}."},
  {NULL, NULL, 0, NULL}
};

/*
 * Module definition.
 */
struct PyModuleDef sgAlignModule = {
  PyModuleDef_HEAD_INIT,
  "sg_align",
  "Library for semi-global alignment.",
  -1,
  pySgAlignMethods
};

/**
 * Module init function.
 */
PyMODINIT_FUNC PyInit_sg_align(void) {
  return PyModule_Create(&sgAlignModule);
}
