#!/usr/bin/env python
"""
TSSV-Lite: Targeted characterisation of short structural variation.


The library file consists of three tab-separated columns:
- name of the marker pair
- marker 1
- marker 2

Without the -d option, TSSV-Lite will write a table containing a count
of each unique sequence detected to standard output, and a table of
per-marker statistics to the standard error stream. This output can be
redirected to a file using the -o and -r options, respectively.

If the -d option is used, a folder will be created instead, and the
output is written to two files named sequences.csv and statistics.csv in
this folder. Additionally, a subfolder is created for each marker
containing separate tables of the detected sequences and split FASTA
files.

Copyright (c) 2016 Leiden University Medical Center <humgen@lumc.nl>
Copyright (c) 2016 Jeroen F.J. Laros <j.f.j.laros@lumc.nl>
Copyright (c) 2016 Jerry Hoogenboom <j.hoogenboom@nfi.minvenj.nl>
Copyright (c) 2012 Jaap W.F. van der Heijden

Licensed under the MIT license, see the LICENSE file.
"""

import argparse
import math
import sys
import os

from Bio import Seq, SeqIO
from fastools import fastools

from . import ProtectedFileType, version
from .align_pair import align_pair

usage = __doc__.split("\n\n\n")

def parse_library(library_handle, threshold):
    """
    Read a library of flanking sequences from an open file handle.

    :arg library_handle: Open readable handle to a library file.
    :type library_handle: stream
    :arg threshold: Number of allowed mismatches per nucleotide.
    :type threshold: float

    :returns: Dict containing library data.
    :rtype: dict{str: tuple(tuple(str, str), tuple(int, int))}
    """
    return {data[0]: (
        (data[1].upper(), Seq.reverse_complement(data[2].upper())), (
            int(math.ceil(len(data[1]) * threshold)),
            int(math.ceil(len(data[2]) * threshold))))
                for data in (
                    line.split() for line in library_handle.readlines())}

def prepare_output_dir(dir, markers, file_format):
    """
    Create output directories and return a nested dict of output files.

    :arg dir: Output directory name.
    :type dir: str
    :arg markers: Iterable providing the marker names.
    :type markers: iterable
    :arg file_format: Output file format ('fasta' or 'fastq').
    :type file_format: str

    :returns: Nested dict with open writable streams to output files.
    :rtype: dict
    """
    # Create output directories.
    os.mkdir(dir)
    for marker in markers:
        os.mkdir(os.path.join(dir, marker))

    # Open output files.
    return {
        "sequences": open(os.path.join(dir, "sequences.csv"), "w"),
        "statistics": open(os.path.join(dir, "statistics.csv"), "w"),
        "unknown": open(os.path.join(dir, "unknown.f" + file_format[-1]), "w"),
        "markers": {
            marker: {
                "sequences": open(os.path.join(dir, marker,
                    "sequences.csv"), "w"),
                "paired": open(os.path.join(dir, marker,
                    "paired.f" + file_format[-1]), "w"),
                "noend": open(os.path.join(dir, marker,
                    "noend.f" + file_format[-1]), "w"),
                "nostart": open(os.path.join(dir, marker,
                    "nostart.f" + file_format[-1]), "w"),
            } for marker in markers
        }
    }


def find_pair(seq, pair, thresholds, indel_score=1):
    """
    Align a pair of markers to the forward reference sequence.  The
    reverse complement is used to align the second element of the pair
    (which is also reverse complemented).

    :arg seq: Sequence to align to, and its reverse complement.
    :type seq: list[str, str]
    :arg pair: A pair (forward, reverse) of markers to align.
    :type pair: list[str, str]
    :arg thresholds: Maximum number of allowed mismatches.
    :type thresholds: list[int, int]
    :arg indel_score: Penalty score for insertions and deletions per
    nucleotide (optional, default: 1)
    :type indel_score: int

    :returns: Bools indicating whether the two flanks are found, and the
    sequence between the two hits (or None if there was no paired hit).
    :rtype: tuple(list[bool, bool], str)
    """
    seq_up = map(str.upper, seq)
    alignments = align_pair(seq_up[0], seq_up[1], pair, indel_score)
    matches = [False, False]

    if alignments[0][0] <= thresholds[0]:
        # Left marker was found.
        cutout = seq[0][max(0, alignments[0][1]-len(pair[0])):alignments[0][1]]
        if cutout.lower() != cutout:
            matches[0] = True
    if alignments[1][0] <= thresholds[1]:
        # Right marker was found.
        cutout = seq[0][alignments[1][1] : alignments[1][1]+len(pair[1])]
        if cutout.lower() != cutout:
            matches[1] = True

    if matches[0] and matches[1] and alignments[0][1] < alignments[1][1]:
        # Matched pair.
        return matches, seq[0][alignments[0][1] : alignments[1][1]]
    return matches, None


def process_file(input_handle, file_format, library, outfiles=None,
                 indel_score=1):
    """
    Process a FASTA or FASTQ file, attempt to link every read to a
    marker with the given library of flanking sequences.

    If given, the outfiles argument should be a nested dict with the
    following keys, each containing an open writable handle to a FASTA
    or FASTQ file:
    outfiles['unknown'],
    outfiles['markers'][<marker name>]['paired'],
    outfiles['markers'][<marker name>]['noend'],
    outfiles['markers'][<marker name>]['nostart']

    :arg input_handle: Open readable stream with FASTA or FASTQ data.
    :type input_handle: stream
    :arg file_format: The name of the file format ('fasta' or 'fastq').
    :type file_format: str
    :arg library: Dict containing library data.
    :type library: dict{str: tuple(tuple(str, str), tuple(int, int))}
    :arg outfiles: Dict with open writable streams to which FASTA or
    FASTQ data will be written (optional).
    :type outfiles: dict
    :arg indel_score: Penalty score for insertions and deletions per
    nucleotide (optional, default: 1)
    :type indel_score: int

    :returns: The total number of reads, the number of unrecognised
    reads, a dict with counts of flanking sequence hits per marker, and
    a dict with sequences per marker.
    :rtype: tuple(int, int, dict{str: dict{str: int}},
    dict{str: dict{str: list[int, int]}})
    """
    sequences = {marker: {} for marker in library}
    counters = {marker: {key: 0 for key in
            ("fPaired", "rPaired", "fLeft", "rLeft", "fRight", "rRight")}
        for marker in library}
    total_reads = 0
    unrecognised = 0
    for record in SeqIO.parse(input_handle, file_format):
        ref = (str(record.seq), Seq.reverse_complement(str(record.seq)))
        total_reads += 1
        recognised = False

        for marker in library:

            # Search in the forward strand.
            matches, seq = find_pair(ref, library[marker][0],
                library[marker][1], indel_score)
            if matches[0]:
                recognised = True
                counters[marker]["fLeft"] += 1
            if matches[1]:
                recognised = True
                counters[marker]["fRight"] += 1
            if seq is not None:
                counters[marker]["fPaired"] += 1
                if seq not in sequences[marker]:
                    sequences[marker][seq] = [1, 0]
                else:
                    sequences[marker][seq][0] += 1
                if outfiles:
                    SeqIO.write([record],
                        outfiles["markers"][marker]["paired"], file_format)
            elif outfiles:
                if matches[0]:
                    SeqIO.write([record],
                        outfiles["markers"][marker]["noend"], file_format)
                if matches[1]:
                    SeqIO.write([record],
                        outfiles["markers"][marker]["nostart"], file_format)

            # Search in the reverse strand.
            matches, seq = find_pair(ref[-1::-1], library[marker][0],
                library[marker][1], indel_score)
            if matches[0]:
                recognised = True
                counters[marker]["rLeft"] += 1
            if matches[1]:
                recognised = True
                counters[marker]["rRight"] += 1
            if seq is not None:
                counters[marker]["rPaired"] += 1
                if seq not in sequences[marker]:
                    sequences[marker][seq] = [0, 1]
                else:
                    sequences[marker][seq][1] += 1
                if outfiles:
                    SeqIO.write([record],
                        outfiles["markers"][marker]["paired"], file_format)
            elif outfiles:
                if matches[0]:
                    SeqIO.write([record],
                        outfiles["markers"][marker]["noend"], file_format)
                if matches[1]:
                    SeqIO.write([record],
                        outfiles["markers"][marker]["nostart"], file_format)

        if not recognised:
            unrecognised += 1
            if outfiles:
                SeqIO.write([record], outfiles["unknown"], file_format)

    # Count number of unique sequences per marker.
    for marker in library:
        counters[marker]["unique_seqs"] = len(sequences[marker])

    # Return counters and sequences.
    return total_reads, unrecognised, counters, sequences


def make_sequence_tables(sequences, minimum=1):
    """
    Return a dict of tables of the counts of all sequences for all
    markers, with a specified minimum number of reads per sequence.

    :arg sequences: Nested dict of read counts for all sequences of all
    markers.
    :type sequences: dict{str: dict{str: list[int, int]}}
    :arg minimum: The minimum number of reads per sequence to include in
    the table.
    :type minimum: int

    :returns: A tuple with two elements: a tuple of column names and a
    dict with a table of sequences per marker.
    :rtype: tuple(tuple(str), dict{str: list[list]})
    """
    header = ("marker", "sequence", "total", "forward", "reverse")
    tables = {}
    for marker in sequences:
        tables[marker] = []
        for total, sequence in sorted(
                ((sum(sequences[marker][s]), s) for s in sequences[marker]),
                reverse=True):
            if total < minimum:
                break
            tables[marker].append(
                [marker, sequence, total] + sequences[marker][sequence])
    return header, tables


def make_sequence_table(sequences, minimum=1, outfiles=None):
    """
    Return a tab-separated table of the counts of all sequences for all
    markers, with a specified minimum number of reads per sequence.

    :arg sequences: Nested dict of read counts for all sequences of all
    markers.
    :type sequences: dict{str: dict{str: list[int, int]}}
    :arg minimum: The minimum number of reads per sequence to include in
    the table.
    :type minimum: int
    :arg outfiles: Dict with the marker names as keys and open writable
    streams as values, for separate tables for each marker (optional).
    :type outfiles: dict{str: stream}

    :returns: A tab-separated table of the input data.
    :rtype: str
    """
    header, tables = make_sequence_tables(sequences, minimum)
    header = "\t".join(header)
    for marker in sorted(tables):
        tables[marker] = "\n".join(
            "\t".join(map(str, line)) for line in tables[marker])
        if outfiles:
            outfiles[marker].write(
                "\n".join((header, tables[marker])))
    return "\n".join([header] + [tables[marker] for marker in sorted(tables)])


def make_statistics_table(counters):
    return "\n".join([
        "\t".join(("marker", "unique_seqs", "tPaired", "fPaired", "rPaired",
            "tLeft", "fLeft", "rLeft", "tRight", "fRight", "rRight"))] + [
        "\t".join(map(str, [
            marker,
            counters[marker]["unique_seqs"],
            counters[marker]["fPaired"] + counters[marker]["rPaired"],
            counters[marker]["fPaired"],
            counters[marker]["rPaired"],
            counters[marker]["fLeft"] + counters[marker]["rLeft"],
            counters[marker]["fLeft"],
            counters[marker]["rLeft"],
            counters[marker]["fRight"] + counters[marker]["rRight"],
            counters[marker]["fRight"],
            counters[marker]["rRight"]])) for marker in sorted(counters)])


def tssv_lite(input_handle, file_format, library, minimum, dir, tee,
              output_handle, report_handle, indel_score):
    """
    """
    if dir:
        outfiles = prepare_output_dir(dir, library, file_format)
    else:
        outfiles = None

    # Process the input file and obtain counters and sequences.
    total_reads, unrecognised, counters, sequences = process_file(
        input_handle, file_format, library, outfiles, indel_score)

    # Make table of sequences.
    table_of_sequences = make_sequence_table(sequences, minimum,
        None if outfiles is None else
        {m: outfiles["markers"][m]["sequences"] for m in outfiles["markers"]})

    # Make table of statistics.
    statistics = "\n".join((
        make_statistics_table(counters),
        "",  # Empty line.
        "total reads\t%i" % total_reads,
        "unrecognised reads\t%i" % unrecognised))

    if not outfiles or tee or output_handle != sys.stdout:
        output_handle.write(table_of_sequences + "\n")
    if not outfiles or tee or report_handle != sys.stderr:
        report_handle.write(statistics + "\n")
    if outfiles:
        outfiles["sequences"].write(table_of_sequences)
        outfiles["statistics"].write(statistics)


def main():
    """
    Main entry point.
    """
    parser = argparse.ArgumentParser(description=usage[0], epilog=usage[1],
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('input_handle', metavar='INPUT',
        type=argparse.FileType('r'), help='input FASTA/FASTQ file')
    parser.add_argument('library_handle', metavar='LIBRARY',
        type=argparse.FileType('r'), help='library of flanking sequences')
    parser.add_argument('-f', '--format', metavar='FORMAT',
        choices=('auto', 'fasta', 'fastq'), default='auto',
        help='file type of the input data: one of %(choices)s '
             '(default=%(default)s)')
    parser.add_argument('-m', '--mismatches', type=float, default=0.08,
        help='mismatches per nucleotide (default=%(default)s)')
    parser.add_argument('-n', '--indel-score', type=int, default=1,
        help='insertions and deletions are penalised this number of times '
             'more heavily than mismatches (default=%(default)s)')
    parser.add_argument('-d', '--dir', type=str,
        help='output directory for verbose output')
    parser.add_argument('-t', '--tee', action='store_true',
        help='if specified together with -d/--dir, output on stdout/stderr '
             'is not suppressed')
    parser.add_argument('-a', '--minimum', type=int, default=0,
        help='minimum count per sequence (default=%(default)s)')
    parser.add_argument('-o', '--output', type=ProtectedFileType('w'),
        default=sys.stdout,
        help='write output to this file (default: write to standard output)')
    parser.add_argument('-r', '--report', type=ProtectedFileType('w'),
        default=sys.stderr,
        help='write a statistical report to this file (default: write to '
             'standard error)')
    parser.add_argument('-v', '--version', action='version',
        version=version(parser.prog))

    args = parser.parse_args()

    try:
        tssv_lite(
            args.input_handle,
            fastools.guess_file_format(args.input_handle)
                if args.format == 'auto' else args.format,
            parse_library(args.library_handle, args.mismatches),
            args.minimum,
            args.dir,
            args.tee,
            args.output,
            args.report,
            args.indel_score)
    except ValueError, error:
        parser.error(error)


if __name__ == '__main__':
    main()
