#include <Python.h>

#include "sgAlign.h"


static PyObject *pyAlignment(int distance, int position) {
  return Py_BuildValue(
    "{s: i, s: i}", "distance", distance, "position", position);
}

static PyObject *pyAlign(PyObject *self, PyObject *args) {
  char *seq1,
       *seq2,
       indel_score;
  alignment a;

  if (!PyArg_ParseTuple(args, "ssb", &seq1, &seq2, &indel_score)) {
    return NULL;
  }

  a = align(seq1, seq2, indel_score);

  return pyAlignment(a.distance, a.position);
}

static PyMethodDef pySgAlignMethods[] = {
  {
    "align", pyAlign, METH_VARARGS,
    "Bla"},
  {NULL, NULL, 0, NULL}
};

static struct PyModuleDef sgAlignModule = {
  PyModuleDef_HEAD_INIT,
  "sg_align",
  "Bla.",
  -1,
  pySgAlignMethods
};

PyMODINIT_FUNC PyInit_sg_align(void) {
  PyObject* module = PyModule_Create(&sgAlignModule);

  return module;
}
