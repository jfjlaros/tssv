"""Convert a csv files containing alleles and counts to HGVS descriptions of
alleles and single variants. Also report statistics about the variant types.


The input file is typically one of the output files of tssv.
"""
from argparse import ArgumentParser, FileType, RawDescriptionHelpFormatter
from collections import defaultdict
from sys import stdout

from requests import get as req_get


def write_table(data, title, report_handle, minimum):
    """Write a table to a file.

    :arg dict data: Dictionary containing counts per type.
    :arg str title: Name of the first column.
    :arg stream report_handle: Open writeable handle to the report file.
    :arg int minimum: Minimum count.
    """
    report_handle.write('{}\ttotal\tforward\treverse\n'.format(title))

    for i in sorted(data, key=lambda x: x[0], reverse=True):
        if data[i][0] < minimum:
            return
        report_handle.write('{}\t{}\t{}\t{}\n'.format(i, *data[i]))


def annotate(alleles_handle, reference, report_handle, minimum):
    """Convert a csv files containing alleles and counts to HGVS descriptions
    of alleles, single variants and variant types.

    :arg stream alleles_handle: Open handle to the alleles file.
    :arg str reference: The reference sequence.
    :arg stream report_handle: Open writeable handle to the report file.
    :arg int minimum: Minimum count.
    """
    alleles = defaultdict(lambda: [0, 0, 0])
    raw_vars = defaultdict(lambda: [0, 0, 0])
    classification = defaultdict(lambda: [0, 0, 0])

    data = list(map(
        lambda x: x.strip('\n').split('\t'), alleles_handle.readlines()[1:]))
    for i in data:
        allele_description = req_get(
            'https://mutalyzer.nl/json/descriptionExtract?' +
            'reference={}&observed={}'.format(reference, i[0])).json()
        encountered = list(map(int, (i[1:])))

        alleles[allele_description['description']] = list(map(
            sum, zip(alleles[allele_description['description']], encountered)))
        for j in allele_description['allele']:
            raw_vars[j['description']] = list(map(
                sum, zip(raw_vars[j['description']], encountered)))
            classification[j['type']] = list(map(
                sum, zip(classification[j['type']], encountered)))

    write_table(alleles, 'allele', report_handle, minimum)
    report_handle.write('\n')
    write_table(raw_vars, 'variant', report_handle, minimum)
    report_handle.write('\n')
    write_table(classification, 'class', report_handle, minimum)


def main():
    """Main entry point."""
    usage = __doc__.split('\n\n\n')
    parser = ArgumentParser(
        description=usage[0], epilog=usage[1],
        formatter_class=RawDescriptionHelpFormatter)

    parser.add_argument(
        'alleles', metavar='alleles', type=FileType('r'),
        help='the alleles file')
    parser.add_argument(
        'reference', metavar='reference', type=str,
        help='sequence of the reference allele')
    parser.add_argument(
        '-r', dest='report', type=FileType('w'), default=stdout,
        help='name of the report file')
    parser.add_argument(
        '-a', dest='minimum', type=int, default=0,
        help='minimum count (default=%(default)s)')

    args = parser.parse_args()

    try:
        annotate(args.alleles, args.reference, args.report, args.minimum)
    except OSError as error:
        parser.error(error)


if __name__ == '__main__':
    main()
