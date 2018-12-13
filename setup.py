from setuptools import setup

from distutils.core import Extension


setup(
    ext_modules=[Extension(
        'tssv._sg_align',
        ['tssv/sg_align.c', 'tssv/sg_align.i'],
        extra_compile_args=['-O3'])]
)
