Usage
=====

The ``tssv`` program does targeted characterisation of short structural
variation. See the help (``-h``) of this program for a full description of the
parameters.


Paired-end read support
-----------------------

With paired-end data in which the insert size was less than twice the read
length, the two reads of each pair may overlap partially or entirely. In this
setting, it is possible to use a read merging tool to combine the two reads of
each pair, possibly giving rise to a longer combined read. The combined read
may contain a known allele that does not fit completely in a single read.

Because repetitive sequences like STRs may perfectly align in many ways, such
tools will be unable to detect whether the combined read actually contains all
repeats that were originally present in the DNA molecule. If the molecule was
too long, the two reads may not overlap at all; however, with a highly
repetitive sequence, they may still appear to overlap significantly. It is only
certain that all repeats have been captured if all of them fitted on a single
read - but the flanking sequences may not be present in their entirety for such
long alleles. Merging the two reads together may still yield the complete
flanking sequence, allowing the allele to be detected by TSSV.

To make sure that any potentially truncated sequences are rejected, the
combined reads can be provided to TSSV with upper case letters for the
overlapped part (i.e., the middle, where both reads in the pair overlap) and
lower case letters for the non-overlapping parts (i.e., those parts only
present in one of the two reads). If TSSV aligns a flanking sequence to a
completely lower case part of the combined read, this means that the flanking
sequence was only present in one of the two reads and had fallen off the end of
the other. In this case, it is impossible to tell whether the sequence between
the flanking sequences contains all repeats that were originally present in the
DNA molecule and therefore it is not meaningful to perform a regular expression
match against it.

One particular paired-end read merging tool that supports this kind of output
is a fork of FLASH 1.2.11 that was specifically made for this purpose, which is
hosted on GitHub_. This fork introduces an optional command-line argument
``--lowercase-overhang`` that, when specified, enables output compatible with
the paired-end read support of TSSV.

Cheat sheet: Interpretation of letter case in TSSV input and output
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Letter case in the library file:

-  All sequences in the library file should be in upper case completely.
   Regular expressions containing lower case letters will never match. Lower
   case letters in flanking sequences will always count as a mismatch in
   alignments.

Interpretation of letter case in input sequences:

-  Completely upper case reads are processed normally.
-  TSSV will be unable to detect any flanking sequences in lower case
   sequences.
-  In mixed case sequences, the regions to which the flanking sequences are
   aligned should contain at least one upper case letter. When a flanking
   sequence aligns to a completely lower case region, the match is rejected.
   When a pair of flanking sequences is found in the read (and neither match is
   rejected), the corresponding regular expression is matched case
   insensitively to the sequence between the flanks.

Letter case in output files:

-  Sequences in output FASTA files have the same case as the corresponding
   input sequences.
-  Sequences in output CSV files will be completely in upper case.

Implementation notes
^^^^^^^^^^^^^^^^^^^^

TSSV will do all alignments and pattern matching against an explicitly upper
case copy of each read. The library should therefore contain only upper case
sequences. When a flanking sequence (as provided in the library file) has been
detected in a read, TSSV will check whether the section of the read that
corresponds to this flanking sequence contains any upper case letters. If this
section of the read was completely lower case, the match is rejected.


References
----------

Seyed Yahya Anvar, Kristiaan J. van der Gaag, Jaap W. F. van der Heijden,
Marcel H. A. M. Veltrop, Rolf H. A. M. Vossen, Rick H. de Leeuw, Cor Breukel,
Henk P. J. Buermans, J. Sjef Verbeek, Peter de Knijff, Johan T. den Dunnen, and
Jeroen F. J. Laros, **TSSV: a tool for characterization of complex allelic
variants in pure and mixed genomes.** Bioinformatics, first published online
February 13, 2014. doi:10.1093/bioinformatics/btu068

- Abstract_.
- `Full text`_.


.. _GitHub: https://github.com/Jerrythafast/FLASH-lowercase-overhang
.. _Abstract: http://bioinformatics.oxfordjournals.org/content/early/2014/02/24/bioinformatics.btu068.abstract
.. _`full text`: http://bioinformatics.oxfordjournals.org/content/early/2014/02/24/bioinformatics.btu068.full.pdf+html
