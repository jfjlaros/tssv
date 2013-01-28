#!/usr/bin/env python

"""
Convert a csv files containing alleles and counts to HGVS descriptions of
alleles and single variants. Also report statistics about the variant types.


The input file is typically one of the output files of tssv.
"""

import sys
import numpy
import argparse
import collections
from suds import client

URL = "https://mutalyzer.nl/services/?wsdl"

def writeTable(data, title, reportHandle, minimum):
    """
    Write a table to a file.

    @arg data: Dictionary containing counts per type.
    @type data: dict
    @arg title: Name of the first column.
    @type title: str
    @arg reportHandle: Open writeable handle to the report file.
    @type reportHandle: stream
    @arg minimum: Minimum count.
    @type minimum: int
    """
    reportHandle.write("%s\ttotal\tforward\treverse\n" % title)

    for i in sorted(data, key=lambda x: x[0], reverse=True):
        if data[i][0] < minimum:
            return

        reportHandle.write("%s\t%i\t%i\t%i\n" % tuple([i] + list(data[i])))
    #for
#writeTable

def annotate(allelesHandle, reference, reportHandle, minimum):
    """
    Convert a csv files containing alleles and counts to HGVS descriptions of
    alleles, single variants and variant types.

    @arg allelesHandle: Open handle to the alleles file.
    @type allelesHandle: stream
    @arg reference: The reference sequence.
    @type reference: str
    @arg reportHandle: Open writeable handle to the report file.
    @type reportHandle: stream
    @arg minimum: Minimum count.
    @type minimum: int
    """
    alleles = collections.defaultdict(lambda: numpy.array([0, 0, 0]))
    rawVars = collections.defaultdict(lambda: numpy.array([0, 0, 0]))
    classification = collections.defaultdict(lambda: numpy.array([0, 0, 0]))
    service = client.Client(URL).service

    data = map(lambda x: x.strip('\n').split('\t'),
        allelesHandle.readlines()[1:])
    for i in data:
        alleleDescription = service.descriptionExtract(reference=reference,
            observed=i[0])
        encountered = map(int, (i[1:]))

        alleles[alleleDescription.description] += encountered
        for j in alleleDescription.allele:
            for k in j[1]:
                rawVars[k.hgvs] += encountered
                classification[k.type] += encountered
            #for
    #for

    writeTable(alleles, "allele", reportHandle, minimum)
    reportHandle.write("\n")
    writeTable(rawVars, "variant", reportHandle, minimum)
    reportHandle.write("\n")
    writeTable(classification, "class", reportHandle, minimum)
#annotate

def main():
    """
    Main entry point.
    """
    usage = __doc__.split("\n\n\n")
    parser = argparse.ArgumentParser(description=usage[0], epilog=usage[1],
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("alleles", metavar="alleles",
        type=argparse.FileType("r"), help="the alleles file")
    parser.add_argument("reference", metavar="reference", type=str,
        help="sequence of the reference allele")
    parser.add_argument("-r", dest="report", type=argparse.FileType("w"),
        default=sys.stdout, help="name of the report file")
    parser.add_argument("-a", dest="minimum", type=int, default=0,
        help="minimum count (default=%(default)s)")

    args = parser.parse_args()

    try:
        annotate(args.alleles, args.reference, args.report, args.minimum)
    except OSError, error:
        parser.error(error)
#main

if __name__ == "__main__":
    main()
