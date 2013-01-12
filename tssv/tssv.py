#!/usr/bin/env python

"""
Targeted characterisation of short structural variation.


The library file consists of four tab-separated columns:
- name of the marker pair
- marker 1
- marker 2
- space-separated description of the expected pattern

If the -d option is used, a folder will be created containing the library
table and per marker a subfolder containing the new alleles and split FASTA
files.
"""

import os
import re
import sys
import math
import argparse
import collections
from Bio import Seq, SeqIO

import sg_align

fileNames = {
    "markers" : "markers.csv"
}
""" Names of the global report files. """

markerFileNames = {
    "known"        : "known.fa",
    "new"          : "new.fa",
    "unknown"      : "unknown.fa",
    "knownalleles" : "knownalleles.csv",
    "newalleles"   : "newalleles.csv"
}
""" Names of the marker specific report files. """

def parseLibrary(libHandle, threshold):
    """
    Parse the library file and put the data in a nested dicrionary containing
    per marker the two forward flanking sequences, the two reverse flanking
    sequences and a regular expression pattern object.

    @arg libHandle: Open readable handle to a library file.
    @type libHandle: stream
    @arg threshold: Number of allowed mismatches per nucleotide.
    @type threshold: int

    @returns: Nested dictionary containing library data.
    @rtype: dict
    """
    library = {}
    data = map(lambda x: x.strip().split('\t'), libHandle.readlines())

    for i in data:
        pat = i[3].split()
        pattern = "^%s$" % ''.join(map(lambda x: ("(%s){%s,%s}" % (pat[x],
            pat[x + 1], pat[x + 2])), range(0, len(pat), 3)))

        library[i[0]] = {
            "flanks": [i[1], Seq.reverse_complement(i[2])],
            "counts": [0, 0, 0, 0],
            "pairMatch": [0, 0],
            "thresholds": [int(math.ceil(len(i[1]) * threshold)),
                int(math.ceil(len(i[2]) * threshold))],
            "threshold": int(math.ceil((len(i[1]) + len(i[2])) * threshold)),
            "regExp": re.compile(pattern),
            "new": collections.defaultdict(lambda: [0, 0]),
            "known": collections.defaultdict(lambda: [0, 0])
        }
    #for

    return library
#parseLibrary

def openFiles(path, markers):
    """
    Make a directory structure and return a nested dictionary containing open
    writable handles to the files in the newly created directory.

    @arg path: Name of the output folder.
    @type path: str
    @arg markers: Name of the subfolders.
    @type markers: list[str]

    @returns: Nested dictionary containing writable file handles.
    @rtype: dict
    """
    os.mkdir(path)
    files = dict(map(lambda x: (x, open("%s/%s" % (path, fileNames[x]), "w")),
        fileNames))
    for i in markers:
        markerPath = "%s/%s" % (path, i)

        os.mkdir(markerPath)
        files[i] = dict(map(lambda x: (x, open("%s/%s" % (markerPath,
            markerFileNames[x]), "w")), markerFileNames))
    #for

    return files
#openFiles

def alignPair(reference, referenceRc, pair):
    """
    Align a pair of markers to the forward reference sequence. The reverse
    complement is used to align the second element of the pair (which is also
    reverse complemented).

    @arg reference: Reference sequence to align to.
    @type reference: str
    @arg referenceRc: Reverse complement of the reference sequence.
    @type referenceRc: str
    @arg pair: A pair (forward, reverse) of markers to align.
    @type pair: list[str, str]

    @returns: A tuple (score, position) of the best alignment.
    @rtype: tuple(int, int)
    """
    rightAl = sg_align.align(referenceRc, pair[1])

    return sg_align.align(reference, pair[0]), (rightAl[0], len(reference) -
        rightAl[1])
#alignPair

def alignFwRev(reference, referenceRc, pair):
    """
    Align a pair of markers to both the forward and the reverse complement of
    a sequence. Return the scores and positions for all four alignments.

    @arg reference: Reference sequence to align to.
    @type reference: str
    @arg referenceRc: Reverse complement of the reference sequence.
    @type referenceRc: str
    @arg pair: A pair (forward, reverse) of markers to align.
    @type pair: list[str, str]

    @returns: Two tuples (score, position) of the alignments.
    @rtype: tuple(tuple(int, int), tuple(int, int))
    """
    return alignPair(reference, referenceRc, pair), alignPair(referenceRc,
        reference, pair)
#alignFwRev

def libTable(library, handle):
    """
    Write the library statistics to a file.

    @arg library: Nested dictionary containing library data.
    @type library: dict
    @arg handle: Open writable handle to the output file.
    @type handle: stream
    """
    handle.write("name\tfPaired\trPaired\tfLeft\trLeft\tfRight\trRight\n")

    for i in library:
        handle.write("%s\t%i\t%i\t%i\t%i\t%i\t%i\n" % tuple(
            [i] + library[i]["pairMatch"] + library[i]["counts"]))
#libTable

def alleleTable(newAl, handle, minimum):
    """
    Write the allele statistics to a file.

    @arg newAl: Dictionary with count data of new alleles.
    @type newAl: dict
    @arg handle: Open writable handle to the output file.
    @type handle: stream
    @arg minimum: Minimum count per allele.
    @type minimum: int
    """
    handle.write("allele\ttotal\tforward\treverse\n")

    for i in sorted(newAl, key=lambda x: sum(newAl[x]), reverse=True):
        if sum(newAl[i]) < minimum:
            return

        handle.write("%s\t%i\t%i\t%i\n" % tuple([i] +
            [sum(newAl[i])] + newAl[i]))
    #for
#alleleTable

def makeReport(total, library, handle, minimum):
    """
    Make an overview of the results.

    @arg total: Total number of reads in the FASTA file.
    @type total: int
    @arg library: Nested dictionary containing library data.
    @type library: dict
    @arg handle: Open writable handle to the report file.
    @type handle: stream
    @arg minimum: Minimum count per allele.
    @type minimum: int
    """
    handle.write("total reads: %i\n" % total)
    handle.write("matched pairs: %i\n" % 
        sum(map(lambda x: sum(library[x]["pairMatch"]), library)))
    handle.write("new alleles: %i\n\n" % 
        sum(map(lambda x: len(library[x]["new"]), library)))

    libTable(library, handle)

    for i in library:
        handle.write("\nknown alleles for marker %s:\n" % i)
        alleleTable(library[i]["known"], handle, minimum)
        handle.write("\nnew alleles for marker %s:\n" % i)
        alleleTable(library[i]["new"], handle, minimum)
    #for
#makeReport

def tssv(fastaHandle, libHandle, reportHandle, path, threshold, minimum):
    """
    Do the short structural variation analysis.

    @arg fastaHandle: Open readable handle to a FASTA file.
    @type fastaHandle: stream
    @arg libHandle: Open readable handle to a library file.
    @type libHandle: stream
    @arg reportHandle: Open writable handle to the report file.
    @type reportHandle: stream
    @arg path: Name of the output folder.
    @type path: str
    @arg threshold: Number of allowed mismatches per nucleotide.
    @type threshold: int
    @arg minimum: Minimum count per allele.
    @type minimum: int
    """
    total = 0
    library = parseLibrary(libHandle, threshold)

    if path:
        files = openFiles(path, library)

    for record in SeqIO.parse(fastaHandle, "fasta"):
        ref = [str(record.seq), Seq.reverse_complement(str(record.seq))]
        total += 1

        for i in library:
            alignments = alignFwRev(ref[0], ref[1], library[i]["flanks"])
            scores = map(lambda x: x[0][0] + x[1][0], alignments)
            classification = "unknown"

            if min(scores) <= library[i]["threshold"]: # A matching pair.
                hit = 0
                if scores[0] > scores[1]:              # Match reverse.
                    hit = 1

                library[i]["pairMatch"][hit] += 1
                pat = ref[hit][alignments[hit][0][1]:alignments[hit][1][1]]

                classification = "new"
                if library[i]["regExp"].match(pat):
                    classification = "known"

                library[i][classification][pat][hit] += 1
            #if

            if alignments[0][0][0] <= library[i]["thresholds"][0]:
                library[i]["counts"][0] += 1
            if alignments[0][1][0] <= library[i]["thresholds"][1]:
                library[i]["counts"][2] += 1
            if alignments[1][0][0] <= library[i]["thresholds"][0]:
                library[i]["counts"][1] += 1
            if alignments[1][1][0] <= library[i]["thresholds"][1]:
                library[i]["counts"][3] += 1

            if path:
                SeqIO.write(record, files[i][classification], "fasta")
        #for
    #for

    if path:
        libTable(library, files["markers"])

        for i in library:
            alleleTable(library[i]["known"], files[i]["knownalleles"], minimum)
            alleleTable(library[i]["new"], files[i]["newalleles"], minimum)
        #for
    #if

    makeReport(total, library, reportHandle, minimum)
#tssv

def main():
    """
    Main entry point.
    """
    usage = __doc__.split("\n\n\n")
    parser = argparse.ArgumentParser(description=usage[0], epilog=usage[1],
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("fasta", metavar="FASTA", type=argparse.FileType("r"),
        help="a FASTA file")
    parser.add_argument("lib", metavar="LIBRARY", type=argparse.FileType("r"),
        help="library of flanking sequences")
    parser.add_argument("-m", dest="mismatches", type=float, default=0.08,
        help="mismatches per nucleotide (default=%(default)s)")
    parser.add_argument("-r", dest="report", type=argparse.FileType("w"),
        default=sys.stdout, help="name of the report file")
    parser.add_argument("-d", dest="path", type=str, help="output directory")
    parser.add_argument("-a", dest="minimum", type=int, default=0,
        help="minimum count per allele (default=%(default)s)")

    args = parser.parse_args()

    try:
        tssv(args.fasta, args.lib, args.report, args.path, args.mismatches,
            args.minimum)
    except OSError, error:
        parser.error(error)
#main

if __name__ == "__main__":
    main()
