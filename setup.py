from setuptools import setup
from distutils.core import Extension


setup(
    ext_modules=[Extension(
        'tssv.sg_align_i',
        ['tssv/sg_align.c', 'tssv/sg_align_i.c'],
        extra_compile_args=['-O3'])]
)
