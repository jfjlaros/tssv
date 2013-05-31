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
    "unknown": "unknown.fa",
    "markers": "markers.csv",
    "known"  : "knownalleles.csv",
    "new"    : "newalleles.csv",
    "nostart": "nostart.csv",
    "noend"  : "noend.csv",
    "summary": "summary.csv"
}
""" Names of the global report files. """

markerFileNames = {
    "known"       : "known.fa",
    "new"         : "new.fa",
    "noend"       : "noend.fa",
    "nostart"     : "nostart.fa",
    "knownalleles": "knownalleles.csv",
    "newalleles"  : "newalleles.csv"
}
""" Names of the marker specific report files. """

headers = {
  "markers"   : "name\tfPaired\trPaired\tfLeft\trLeft\tfRight\trRight\n",
  "allele"    : "allele\ttotal\tforward\treverse\n",
  "nostartend": "name\tforward\treverse\ttotal\n",
  "overview"  : "name\tforward\treverse\ttotal\tallele\n"
}
""" Headers for various tables."""

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

def writeTable(table, header, handle):
    """
    General function for saving tables.

    @arg table: Table content.
    @type table: list
    @arg header: Table header.
    @type header: str
    @arg handle: Open writable handle to the output file.
    @type handle: stream
    """
    if header:
        handle.write(header)

    if table:
        for i in table:
            handle.write("%s\n" % '\t'.join(map(str, i)))
#writeTable

def rewrite(regExp, pattern):
    """
    Make a pattern that matches a regular expression more human readable.

    @arg regExp: A compiled regular expression object.
    @type regExp: object
    @arg pattern: A pattern that matches {regExp}.
    @type pattern: str

    @returns: A human readable version of {pattern}.
    @rtype: str
    """
    newPattern = ""
    match = regExp.match(pattern)

    regs = [((0, 0), None)] + filter(lambda y: y != ((-1, -1), None),
        map(lambda x: (match.regs[x], match.group(x)),
        range(1, len(match.regs))))

    for i in range(len(regs) - 1):
        newPattern += "%s(%i)" % (regs[i + 1][1], (regs[i + 1][0][1] -
            regs[i][0][1]) / (regs[i + 1][0][1] - regs[i + 1][0][0]))

    return newPattern
#rewrite

def alleleTable(newAl, minimum):
    """
    Make an allele statistics table.

    @arg newAl: Dictionary with count data of new alleles.
    @type newAl: dict
    @arg minimum: Minimum count per allele.
    @type minimum: int

    @returns: Allele statistics table.
    @rtype: list
    """
    l = []

    for i in sorted(newAl, key=lambda x: sum(newAl[x]), reverse=True):
        if sum(newAl[i]) < minimum:
            return

        l.append([i] + [sum(newAl[i])] + newAl[i])
    #for

    return l
#alleleTable

def makeTables(total, unrecognised, library, minimum):
    """
    Make overview tables of the results.

    @arg total: Total number of reads in the FASTA file.
    @type total: int
    @arg unrecognised: Number of unrecognised reads in.
    @type unrecognised: int
    @arg library: Nested dictionary containing library data.
    @type library: dict
    @arg minimum: Minimum count per allele.
    @type minimum: int

    @returns: A nested dictionary containing overview tables.
    @rtype: dict
    """
    known = []
    new = []
    noStart = []
    noEnd = []

    tables = {
        "library":  map(lambda x: [x] + library[x]["pairMatch"] +
            library[x]["counts"], library),
        "allele" : collections.defaultdict(dict)
    }

    for i in library:
        for j in library[i]["known"]:
            fr = library[i]["known"][j]
            known.append([i] + fr + [sum(fr), j])
        #for
        for j in library[i]["new"]:
            fr = library[i]["new"][j]
            new.append([i] + fr + [sum(fr), j])
        #for

        noStart.append([i,
            library[i]["counts"][2] - library[i]["pairMatch"][0],
            library[i]["counts"][3] - library[i]["pairMatch"][1]])
        noEnd.append([i,
            library[i]["counts"][0] - library[i]["pairMatch"][0],
            library[i]["counts"][1] - library[i]["pairMatch"][1]])

        tables["allele"][i]["known"] = alleleTable(library[i]["known"],
            minimum)
        tables["allele"][i]["new"] = alleleTable(library[i]["new"], minimum)
    #for

    tables["known"] = sorted(known, key=lambda x: (x[0], x[4]))
    tables["new"] = sorted(new, key=lambda x: (x[0], x[3]), reverse=True)
    tables["nostart"] = map(lambda x: x + [sum(x[1:])], sorted(noStart))
    tables["noend"] = map(lambda x: x + [sum(x[1:])], sorted(noEnd))

    tables["summary"] = [
        ["total reads", total],
        ["matched pairs", sum(map(lambda x: sum(library[x]["pairMatch"]),
            library))],
        ["new alleles",  sum(map(lambda x: x[3], tables["new"]))],
        ["new unique alleles", sum(map(lambda x: len(library[x]["new"]),
            library))],
        ["no start", sum(map(lambda x: x[3], tables["nostart"]))],
        ["no end", sum(map(lambda x: x[3], tables["noend"]))],
        ["unrecognised reads", unrecognised]
    ]

    return tables
#makeTables

def makeReport(tables, handle):
    """
    Make an overview of the results.

    @arg tables: A nested dictionary containing overview tables.
    @type tables: dict
    @arg handle: Open writable handle to the report file.
    @type handle: stream
    """
    writeTable(tables["summary"], "", handle)
    handle.write("\n")

    writeTable(tables["library"], headers["markers"], handle)

    for i in tables["allele"]:
        handle.write("\nknown alleles for marker %s:\n" % i)
        writeTable(tables["allele"][i]["known"], headers["allele"], handle)

        meanLength = 0
        sumOfLengths = sum(map(lambda x: len(x[0]) * x[1],
            tables["allele"][i]["new"]))
        numberOfAlleles = sum(map(lambda x: x[1], tables["allele"][i]["new"]))
        if numberOfAlleles:
            meanLength = sumOfLengths / numberOfAlleles

        handle.write("\nnew alleles for marker %s (mean length %i):\n" % (i,
            meanLength))
        writeTable(tables["allele"][i]["new"], headers["allele"], handle)
    #for
#makeReport

def writeFiles(tables, files):
    """
    Write the overview tables to the appropriate files.

    @arg tables: A nested dictionary containing overview tables.
    @type tables: dict
    @arg files: Nested dictionary containing writable file handles.
    @type files: dict
    """
    writeTable(tables["summary"], "", files["summary"])
    writeTable(tables["library"], headers["markers"], files["markers"])
    writeTable(tables["known"], headers["overview"], files["known"])
    writeTable(tables["new"], headers["overview"], files["new"])
    writeTable(tables["nostart"], headers["nostartend"], files["nostart"])
    writeTable(tables["noend"], headers["nostartend"], files["noend"])

    for i in tables["allele"]:
        writeTable(tables["allele"][i]["known"], headers["allele"],
            files[i]["knownalleles"])
        writeTable(tables["allele"][i]["new"], headers["allele"],
            files[i]["newalleles"])
    #for
#writeFiles

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
    unrecognised = 0
    library = parseLibrary(libHandle, threshold)

    if path:
        files = openFiles(path, library)

    for record in SeqIO.parse(fastaHandle, "fasta"):
        ref = [str(record.seq), Seq.reverse_complement(str(record.seq))]
        total += 1
        unknown = True

        for i in library:
            alignments = (alignPair(ref[0], ref[1], library[i]["flanks"]),
                alignPair(ref[1], ref[0], library[i]["flanks"]))
            matches = [False, False, False, False]
            classification = ""

            if alignments[0][0][0] <= library[i]["thresholds"][0]:
                library[i]["counts"][0] += 1
                classification = "noend"
                matches[0] = True
            #if
            if alignments[0][1][0] <= library[i]["thresholds"][1]:
                library[i]["counts"][2] += 1
                classification = "nostart"
                matches[1] = True
            #if
            if alignments[1][0][0] <= library[i]["thresholds"][0]:
                library[i]["counts"][1] += 1
                classification = "noend"
                matches[2] = True
            #if
            if alignments[1][1][0] <= library[i]["thresholds"][1]:
                library[i]["counts"][3] += 1
                classification = "nostart"
                matches[3] = True
            #if

            if (matches[0] and matches[1]) or (matches[2] and matches[3]):
                hit = int(matches[2] and matches[3])

                library[i]["pairMatch"][hit] += 1
                pat = ref[hit][alignments[hit][0][1]:alignments[hit][1][1]]

                classification = "new"
                if library[i]["regExp"].match(pat):
                    classification = "known"

                library[i][classification][pat][hit] += 1
            #if

            if classification:
                unknown = False

                if path:
                    SeqIO.write([record], files[i][classification], "fasta")
            #if
        #for

        if unknown:
            unrecognised += 1

            if path:
                SeqIO.write([record], files["unknown"], "fasta")
        #if
    #for

    tables = makeTables(total, unrecognised, library, minimum)

    # Make the known alleles more human readable.
    for i in tables["allele"]:
        #if tables["allele"][i]["known"]:
        for j in tables["allele"][i]["known"]:
            j[0] = rewrite(library[i]["regExp"], j[0])
    for i in tables["known"]:
        i[4] = rewrite(library[i[0]]["regExp"], i[4])

    if path:
        writeFiles(tables, files)

    makeReport(tables, reportHandle)
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
