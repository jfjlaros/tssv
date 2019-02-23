from setuptools import setup
from distutils.core import Extension


setup(
    ext_modules=[Extension(
        'tssv.sg_align',
        ['tssv/sgAlignWrapper.c', 'tssv/sgAlign.c', 'tssv/sgAlignSSE.c'],
        extra_compile_args=['-O3'])]
)
