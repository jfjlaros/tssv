#!/usr/bin/env python

"""
Header.


Footer.
"""

import argparse
import os
import re
from Bio import Seq, SeqIO

import sg_align

fileNames = {
    "report"     : "report.txt",
    "nomarker"   : "NoFoundMarker.txt",
    "diffstruct" : "DifferentStructures.txt"
}
""" Names of the global report files. """

markerFileNames = {
    "noend"           : "NoEnd.txt",
    "nobegin"         : "NoBeginning.txt",
    "newallele"       : "NewAllele.txt",
    "goodallele"      : "GoodAllele.txt",
    "difforientation" : "DifferentOrientation.txt"
}
""" Names of the marker specific report files. """

def parseLibrary(libHandle):
    """
    Parse the library file and put the data in a nested dicrionary containing
    per marker the two forward flanking sequences, the two reverse flanking
    sequences and a regular expression pattern object.

    @arg libHandle: Open readable handle to a library file.
    @type libHandle: stream

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
            "regExp": re.compile(pattern)
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
    @arg reference:
    @type reference: str
    @arg referenceRc:
    @type referenceRc: str
    @arg pair:
    @type pair: list[str, str]
    """
    rightAl = sg_align.align(referenceRc, pair[1])

    return sg_align.align(reference, pair[0]), (rightAl[0], len(reference) -
        rightAl[1])
#alignPair

def alignFwRev(reference, referenceRc, pair):
    """
    @arg reference:
    @type reference: str
    @arg referenceRc:
    @type referenceRc: str
    @arg pair:
    @type pair: list[str, str]

    @returns:
    @rtype: tuple(tuple(int, int), tuple(int, int))
    """
    return alignPair(reference, referenceRc, pair), alignPair(referenceRc,
        reference, pair)
#alignFwRev

def amplicon(fastaHandle, libHandle, path, threshold):
    """
    Do the amplicon analysis.

    @arg fastaHandle: Open readable handle to a FASTA file.
    @type fastaHandle: stream
    @arg libHandle: Open readable handle to a library file.
    @type libHandle: stream
    @arg path: Name of the output folder.
    @type path: str
    @arg threshold:
    @type threshold: int
    """
    library = parseLibrary(libHandle)

    for record in SeqIO.parse(fastaHandle, "fasta"):
        revComp = Seq.reverse_complement(record.seq)

        for i in library:
            #print "%s\n%s %s\n" % (record.seq,
            #    library[i]["flanks"][0], library[i]["flanks"][1])
            alignments = alignFwRev(record.seq, revComp, library[i]["flanks"])
            scores = map(lambda x: x[0][0] + x[1][0], alignments)

            if min(scores) <= threshold:   # A matching pair.
                if scores[0] <= scores[1]: # Match forward.
                    pat = record.seq[alignments[0][0][1]:alignments[0][1][1]]
                else:                      # Match reverse.
                    pat = revComp[alignments[1][0][1]:alignments[1][1][1]]

                if library[i]["regExp"].match(str(pat)):
                    print "matching known allele"
            #if
            else:
                pass
            #else
        #for
    #for

    #openFiles(path, library)
#amplicon

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
    parser.add_argument("-m", dest="mismatches", type=int, default=2,
        help="allowed mismatches per 25nt (default=%(default)s)")
    parser.add_argument("-d", dest="path", type=str, default="./out",
        help="output directory (default=%(default)s)")

    args = parser.parse_args()

    try:
        amplicon(args.fasta, args.lib, args.path)
    except OSError, error:
        parser.error(error)
#main

if __name__ == "__main__":
    main()
