from argparse import ArgumentParser, FileType, RawDescriptionHelpFormatter
from sys import stdout
from xopen import xopen

from . import usage, version
from .tssv import tssv


def main():
    """Main entry point."""
    parser = ArgumentParser(
        description=usage[0], epilog=usage[1],
        formatter_class=RawDescriptionHelpFormatter)

    parser.add_argument(
        'input_handle', metavar='INPUT',
        help='a FASTA/FASTQ file')
    parser.add_argument(
        'library_handle', metavar='LIBRARY', type=FileType('r'),
        help='library of flanking sequences')
    parser.add_argument(
        '-m', dest='threshold', type=float, default=0.08,
        help='mismatches per nucleotide (default=%(default)s)')
    parser.add_argument(
        '-M', dest='mismatches', type=int,
        help='fixed number of mismatches, overrides -m (default=%(default)s)')
    parser.add_argument(
        '-n', dest='indel_score', type=int, default=1,
        help='insertions and deletions are penalised this number of times '
             'more heavily than mismatches (default=%(default)s)')
    parser.add_argument(
        '-r', dest='report_handle', type=FileType('w'), default=stdout,
        help='name of the report file')
    parser.add_argument('-d', dest='path', type=str, help='output directory')
    parser.add_argument(
        '-a', dest='minimum', type=int, default=0,
        help='minimum count per allele (default=%(default)s)')
    parser.add_argument(
        '-s', dest='method_sse', action='store_true',
        help='if specified, use SSE2 alignment implementation')
    parser.add_argument('-v', action='version', version=version(parser.prog))

    args = parser.parse_args()

    # Have a little look in the input file to determine the file format.
    with xopen(args.input_handle, 'r') as fin:
        if next(fin).startswith('>'):
            args.file_format = 'fasta'
        else:
            args.file_format = 'fastq'

    # Now that we we know the file format, we can open the file again and
    # have access to the full file content.
    args.input_handle = xopen(args.input_handle)

    try:
        tssv(**{k: v for k, v in vars(args).items()
             if k not in ('func', 'subcommand')})
    except OSError as error:
        parser.error(error)


if __name__ == '__main__':
    main()
