import sys
from setuptools import setup

if sys.version_info < (2, 6):
    raise Exception("TSSV requires Python 2.6 or higher.")

# Todo: How does this play with pip freeze requirement files?
requires = ["biopython", "suds"]

# Python 2.6 does not include the argparse module.
try:
    import argparse
except ImportError:
    requires.append("argparse")

import tssv as distmeta

setup(
    name="tssv",
    version=distmeta.__version__,
    description="Targeted characterisation of short structural variation.",
    long_description=distmeta.__doc__,
    author=distmeta.__author__,
    author_email=distmeta.__contact__,
    url=distmeta.__homepage__,
    license="MIT License",
    platforms=["any"],
    packages=["tssv"],
    install_requires=requires,
    entry_points = {
        "console_scripts": [
            "tssv = tssv.tssv:main",
            "tannotate = tssv.annotate:main",
        ]
    },
    classifiers = [
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering",
    ],
    keywords="bioinformatics"
)
