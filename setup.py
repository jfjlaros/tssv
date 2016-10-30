# Patch for swig.
from distutils.command.build import build
from setuptools.command.install import install

class CustomBuild(build):
    def run(self):
        self.run_command('build_ext')
        build.run(self)


class CustomInstall(install):
    def run(self):
        self.run_command('build_ext')
        self.do_egg_install()

import sys
from setuptools import setup
from distutils.core import Extension

if sys.version_info < (2, 6):
    raise Exception('TSSV requires Python 2.6 or higher.')

# Todo: How does this play with pip freeze requirement files?
requires = ['biopython', 'fastools', 'future', 'requests']
tests_requires = ['fake-open', 'pytest', 'tox']

# Python 2.6 does not include the argparse module.
try:
    import argparse
except ImportError:
    requires.append('argparse')

import tssv as distmeta

setup(
    cmdclass={'build': CustomBuild, 'install': CustomInstall},
    name='tssv',
    ext_modules=[Extension('tssv/_sg_align',
        ['tssv/sg_align.c', 'tssv/sg_align.i'],
        extra_compile_args=['-O3'])],
    version=distmeta.__version__,
    description='Targeted characterisation of short structural variation.',
    long_description=distmeta.__doc__,
    author=distmeta.__author__,
    author_email=distmeta.__contact__,
    url=distmeta.__homepage__,
    license='MIT License',
    platforms=['any'],
    packages=['tssv'],
    install_requires=requires,
    tests_require=tests_requires,
    package_data={'tssv': ['sg_align.h']},
    entry_points = {
        'console_scripts': [
            'tssv = tssv.tssv:main',
            'tssvl = tssv.tssv_lite:main',
            'tannotate = tssv.annotate:main',
        ]
    },
    classifiers = [
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
    ],
    keywords='bioinformatics'
)
