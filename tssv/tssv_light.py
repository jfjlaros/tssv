#!/usr/bin/env python


import argparse
import math

from Bio import Seq, SeqIO
from fastools import fastools

from . import version
from .align_pair import align_pair


def tssv_light(input_handle, library_handle, mismatches):
    """
    """
    library = []
    for line in library_handle.readlines():
        data = line.split()

        library.append((data[0], (data[1].upper(), data[2].upper()), (
            int(math.ceil(len(data[1]) * mismatches)),
            int(math.ceil(len(data[2]) * mismatches)))))

    file_format = fastools.guess_file_format(input_handle)
    files = dict(map(
        lambda x: (x, open('{}.f{}'.format(x[0], file_format[-1]), 'w')), 
        library))

    for record in SeqIO.parse(input_handle, file_format):
        ref = [
            str(record.seq).upper(),
            Seq.reverse_complement(str(record.seq).upper())]

        for marker in library:
            alignments = (
                align_pair(ref[0], ref[1], marker[1])
                align_pair(ref[1], ref[0], marker[1]))


def main():
    """
    Main entry point.
    """
    usage = ['', '']
    parser = argparse.ArgumentParser(description=usage[0], epilog=usage[1],
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('input_handle', metavar='INPUT',
        type=argparse.FileType('r'), help='input FASTA/FASTQ file')
    parser.add_argument('library_handle', metavar='LIBRARY',
        type=argparse.FileType('r'), help='library of flanking sequences')
    parser.add_argument('-m', dest='mismatches', type=float, default=0.08,
        help='mismatches per nucleotide (default=%(default)s)')
    parser.add_argument('-v', action='version', version=version(parser.prog))

    args = parser.parse_args()

    try:
        tssv_light(args.input_handle, args.library_handle, args.mismatches)
    except ValueError, error:
        parser.error(error)


if __name__ == '__main__':
    main()
