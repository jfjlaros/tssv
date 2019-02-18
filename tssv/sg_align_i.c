#include <Python.h>

#include "sg_align.h"


static PyObject *alignment_i(int distance, int position) {
  return Py_BuildValue(
    "{s: i, s: i}", "distance", distance, "position", position);
}


static PyObject *align_i(PyObject *self, PyObject *args) {
  char *seq1,
       *seq2,
       indel_score;
  alignment a;

  if (!PyArg_ParseTuple(args, "ssb", &seq1, &seq2, &indel_score)) {
    return NULL;
  }

  a = align(seq1, seq2, indel_score);

  return alignment_i(a.distance, a.position);
}

static PyMethodDef sg_align_i_methods[] = {
  {
    "align_i", align_i, METH_VARARGS,
    "Bla"},
  {NULL, NULL, 0, NULL}
};


static struct PyModuleDef sg_align_i_module = {
  PyModuleDef_HEAD_INIT,
  "sg_align_i",
  "Bla.",
  -1,
  sg_align_i_methods
};


PyMODINIT_FUNC PyInit_sg_align_i(void) {
  PyObject* module = PyModule_Create(&sg_align_i_module);

  return module;
}
