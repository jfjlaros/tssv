Installation
============

TSSV depends on ``gcc`` for the compilation of one of the core libraries. Also
a package that provides ``Python.h`` should be installed. For Debian based
systems, the following command will install these packages:

::

    apt-get install libpython3-dev gcc

On Windows and macOS, you will be prompted to install a C compiler while
installing TSSV if you do not have one.

The software is distributed via PyPI_, it can be installed with ``pip``:

::

    pip install tssv

If there are multiple versions of Python installed, it might be better to use
the command ``pip3`` instead of ``pip``.


From source
-----------

The source is hosted on GitHub_, to install the latest development version, use
the following commands.

::

    git clone https://github.com/jfjlaros/tssv.git
    cd tssv
    pip install .


.. _swig: http://swig.org/
.. _PyPI: https://pypi.org/project/tssv
.. _GitHub: https://github.com/jfjlaros/tssv.git
