# TSSV changelog

17 November, 2016:
* Support for Windows and MacOS.

30 October, 2016:
* Released version 0.4.1.

30 October, 2016:
* Python 3 support.
* Full testing with tox.
* Implemented more robust SWIG fix.
* Fixed compatibility issue with MSVC compiler.

23 August, 2016:
* Released version 0.4.0.

23 August, 2016:
* Added unit tests.
* Switched to continuous 2d-array in original implementation.
* Enabled more compiler optimisations.
* Removed non-transposed code (never used).
* General refactoring.

9 June, 2016:
* Added -n/--indel-score argument.
* Fixed bug that caused TSSV to install the compiled _sg_align module into the
  system site-packages or dist-packages directory instead of the tssv package
  directory.

22 February, 2016:
* Released version 0.3.0.

25 January, 2016:
* Introduced TSSV-Light (`tssvl`), which is a stripped-down version of TSSV
  without the STR regular expressions and the distinction between known and new
  alleles, and with more condensed output.

24 December, 2015:
* Increased speed by at least 3.5x* on systems supporting SSE2 instructions,
  and by at least 1.5x on other systems by overhauling the alignment algorithm.
  (*Tested on real-world data consisting of 830,000 reads from 24 markers.
  94% of reads were 200-300bp in length, library sequences were 18bp.)

20 February, 2015:
* Repeated units like X(n)X(m) are now collapsed to X(n+m).

18 February, 2015:
* Fixed bug in the detection of the left flank.

17 February, 2015:
* Introduced [paired-end read support](doc/paired-end.md).
* Added -q option to read data from FASTQ files.

19 August, 2014:
* Fixed memory leak.

8 August, 2014:
* Added minimum thresholds to master allele tables and summary.
* Released version 0.2.2.

7 August 2014:
* Fixed bug in allele statistics table generation, when a minimum allele count
  was set, the program crashed.
* Released version 0.2.1.

15 August 2013:
* Testing complete.
* Released first non-beta version (0.2.0).

31 May 2013:
* Minor bug fixes.
* Released version 0.2.dev.

29 April 2013:
* Added annotation functionality.
* Changed the method of identifying a matched marker pair. We do not look at
  the combined alignment score anymore, instead we require both markers to be
  aligned.
* Added human readable descriptions for known alleles.
* Released version 0.1.dev.

12 January 2013:
* Start of log.
* Made first installable package.
