# TSSV changelog

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
