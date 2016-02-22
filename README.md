# TSSV documentation and information
TSSV is a program that does targeted characterisation of short structural
variation. It can be used for STR analysis, or any other type of targeted
analysis. It characterises any variation between a set of user-defined markers.

Development of TSSV is hosted
[on our GitLab server](https://git.lumc.nl/j.f.j.laros/tssv).


## Easy installation
TSSV depends on `swig` and `gcc` for the compilation of one of the core
libraries. Also a package that provides `Python.h` should be installed. For
Debian based systems, the following command will install these packages:

    apt-get install libpython-dev swig gcc

If you have installed the `pip` package installer, you can easily install
TSSV by typing:

    pip install tssv

## Installing from source
To get the source code, use:

    git clone https://git.lumc.nl/j.f.j.laros/tssv.git

When installing from source, make sure you have the following packages
installed:
* A Python interpreter (versions 2.6 or higher).
* The `python-argparse`, `python-suds` and `python-biopython` libraries should
  be installed.

## Releases
Version                                                                       | Notes
---                                                                           | ---
[0.3.0](https://pypi.python.org/packages/source/t/tssv/tssv-0.3.0.tar.gz)     | Another big performance boost, and introducing TSSV-Lite (`tssvl`).
[0.2.5](https://pypi.python.org/packages/source/t/tssv/tssv-0.2.5.tar.gz)     | Introduced [paired-end read support](doc/paired-end.md).
[0.2.4](https://pypi.python.org/packages/source/t/tssv/tssv-0.2.4.tar.gz)     | Fixed memory leak.
[0.2.3](https://pypi.python.org/packages/source/t/tssv/tssv-0.2.3.tar.gz)     | Big [performance boost](doc/benchmark.md).
[0.2.2](https://pypi.python.org/packages/source/t/tssv/tssv-0.2.2.tar.gz)     | Applied minimum threshold to summary tables.
[0.2.1](https://pypi.python.org/packages/source/t/tssv/tssv-0.2.1.tar.gz)     | Fixed minimum threshold bug.
[0.2.0](https://pypi.python.org/packages/source/t/tssv/tssv-0.2.0.tar.gz)     | First official release.
[0.2.dev](https://pypi.python.org/packages/source/t/tssv/tssv-0.2.dev.tar.gz) | Beta test.
[0.1.dev](https://pypi.python.org/packages/source/t/tssv/tssv-0.1.dev.tar.gz) | Beta test.

Please note that the [latest](https://pypi.python.org/pypi/tssv) version will
always be available on the [Python Package Index](https://pypi.python.org/).

See the [changelog](CHANGELOG.md) for a record of changes made to TSSV per
release version.

## Usage
The `tssv` program does targeted characterisation of short structural
variation. See the help (`-h`) of this program for a full description of the
parameters.

## Feedback
Please submit feature requests and error reports using the
[TSSV bug tracking system](https://git.lumc.nl/j.f.j.laros/tssv/issues).

Clicking "All" will also show solved bugs and implemented feature requests.

Please click "New Issue" in order to submit bugs or feature requests.

## References
Seyed Yahya Anvar, Kristiaan J. van der Gaag, Jaap W. F. van der Heijden,
Marcel H. A. M. Veltrop, Rolf H. A. M. Vossen, Rick H. de Leeuw, Cor Breukel,
Henk P. J. Buermans, J. Sjef Verbeek, Peter de Knijff, Johan T. den Dunnen, and
Jeroen F. J. Laros, **TSSV: a tool for characterization of complex allelic
variants in pure and mixed genomes.** Bioinformatics, first published online
February 13, 2014. doi:10.1093/bioinformatics/btu068

[abstract](http://bioinformatics.oxfordjournals.org/content/early/2014/02/24/bioinformatics.btu068.abstract),
[full text](http://bioinformatics.oxfordjournals.org/content/early/2014/02/24/bioinformatics.btu068.full.pdf+html).
