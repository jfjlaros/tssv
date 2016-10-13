import argparse
import collections
import math
import os
import re
import sys

from Bio import Seq, SeqIO

from . import ProtectedFileType, doc_split, usage, version
from .align_pair import align_pair


file_names = {
    "unknown": "unknown.fa",
    "markers": "markers.csv",
    "known"  : "knownalleles.csv",
    "new"    : "newalleles.csv",
    "nostart": "nostart.csv",
    "noend"  : "noend.csv",
    "summary": "summary.csv"
}
""" Names of the global report files. """

marker_file_names = {
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


def parse_library(library_handle, threshold):
    """
    Parse the library file and put the data in a nested dictionary containing
    per marker the two forward flanking sequences, the two reverse flanking
    sequences and a regular expression pattern object.

    :arg library_handle: Open readable handle to a library file.
    :type library_handle: stream
    :arg threshold: Number of allowed mismatches per nucleotide.
    :type threshold: float

    :returns: Nested dictionary containing library data.
    :rtype: dict
    """
    library = {}
    data = map(lambda x: x.strip().split('\t'), library_handle.readlines())

    for i in data:
        pat = i[3].split()
        pattern = "^%s$" % ''.join(map(lambda x: ("(%s){%s,%s}" % (pat[x],
            pat[x + 1], pat[x + 2])), range(0, len(pat), 3)))

        library[i[0]] = {
            "flanks": [i[1], Seq.reverse_complement(i[2])],
            "counts": [0, 0, 0, 0],
            "pair_match": [0, 0],
            "thresholds": [int(math.ceil(len(i[1]) * threshold)),
                int(math.ceil(len(i[2]) * threshold))],
            "reg_exp": re.compile(pattern),
            "new": collections.defaultdict(lambda: [0, 0]),
            "known": collections.defaultdict(lambda: [0, 0])
        }

    return library


def open_files(path, markers):
    """
    Make a directory structure and return a nested dictionary containing open
    writable handles to the files in the newly created directory.

    :arg path: Name of the output folder.
    :type path: str
    :arg markers: Name of the subfolders.
    :type markers: list[str]

    :returns: Nested dictionary containing writable file handles.
    :rtype: dict
    """
    os.mkdir(path)
    files = dict(map(lambda x: (x, open("%s/%s" % (path, file_names[x]), "w")),
        file_names))
    for i in markers:
        marker_path = "%s/%s" % (path, i)

        os.mkdir(marker_path)
        files[i] = dict(map(lambda x: (x, open("%s/%s" % (marker_path,
            marker_file_names[x]), "w")), marker_file_names))

    return files


def write_table(table, header, handle):
    """
    General function for saving tables.

    :arg table: Table content.
    :type table: list
    :arg header: Table header.
    :type header: str
    :arg handle: Open writable handle to the output file.
    :type handle: stream
    """
    if header:
        handle.write(header)

    if table:
        for i in table:
            handle.write("%s\n" % '\t'.join(map(str, i)))


def rewrite(regular_expression, pattern):
    """
    Make a pattern that matches a regular expression more human readable.

    :arg regular_expression: A compiled regular expression object.
    :type regular_expression: object
    :arg pattern: A pattern that matches {regular_expression}.
    :type pattern: str

    :returns: A human readable version of {pattern}.
    :rtype: str
    """
    new_pattern = ""
    match = regular_expression.match(pattern)

    regs = reduce(lambda x, y:
            x if y == ((-1, -1), None) else
            x[:-1] + [y] if x[-1][1] == y[1] else
            x + [y],
        map(lambda x: (match.regs[x], match.group(x)),
        range(1, len(match.regs))), [((0, 0), None)])

    for i in range(len(regs) - 1):
        new_pattern += "%s(%i)" % (regs[i + 1][1], (regs[i + 1][0][1] -
            regs[i][0][1]) / (regs[i + 1][0][1] - regs[i + 1][0][0]))

    return new_pattern


def allele_table(new_allele, minimum):
    """
    Make an allele statistics table.

    :arg new_allele: Dictionary with count data of new alleles.
    :type new_allele: dict
    :arg minimum: Minimum count per allele.
    :type minimum: int

    :returns: Allele statistics table.
    :rtype: list
    """
    l = []

    for i in sorted(new_allele, key=lambda x: sum(new_allele[x]),
            reverse=True):
        if sum(new_allele[i]) < minimum:
            break

        l.append([i] + [sum(new_allele[i])] + new_allele[i])

    return l


def summary_table(allele, minimum):
    """
    Filter one of the global allele tables.

    :arg allele: List with count data of alleles.
    :type allele: list
    :arg minimum: Minimum count per allele.
    :type minimum: int

    :returns: Allele statistics table.
    :rtype: list
    """
    return filter(lambda x: x[3] >= minimum, allele)


def make_tables(total, unrecognised, library, minimum):
    """
    Make overview tables of the results.

    :arg total: Total number of reads in the FASTA file.
    :type total: int
    :arg unrecognised: Number of unrecognised reads in.
    :type unrecognised: int
    :arg library: Nested dictionary containing library data.
    :type library: dict
    :arg minimum: Minimum count per allele.
    :type minimum: int

    :returns: A nested dictionary containing overview tables.
    :rtype: dict
    """
    known = []
    new = []
    no_start = []
    no_end = []

    tables = {
        "library":  map(lambda x: [x] + library[x]["pair_match"] +
            library[x]["counts"], library),
        "allele" : collections.defaultdict(dict)
    }

    for i in library:
        for j in library[i]["known"]:
            fr = library[i]["known"][j]
            known.append([i] + fr + [sum(fr), j])
        for j in library[i]["new"]:
            fr = library[i]["new"][j]
            new.append([i] + fr + [sum(fr), j])

        no_start.append([i,
            library[i]["counts"][2] - library[i]["pair_match"][0],
            library[i]["counts"][3] - library[i]["pair_match"][1]])
        no_end.append([i,
            library[i]["counts"][0] - library[i]["pair_match"][0],
            library[i]["counts"][1] - library[i]["pair_match"][1]])

        tables["allele"][i]["known"] = allele_table(library[i]["known"],
            minimum)
        tables["allele"][i]["new"] = allele_table(library[i]["new"], minimum)

    tables["known"] = sorted(summary_table(known, minimum),
        key=lambda x: (x[0], x[4]))
    tables["new"] = sorted(summary_table(new, minimum),
        key=lambda x: (x[0], x[3]), reverse=True)
    tables["nostart"] = map(lambda x: x + [sum(x[1:])], sorted(no_start))
    tables["noend"] = map(lambda x: x + [sum(x[1:])], sorted(no_end))

    tables["summary"] = [
        ["total reads", total],
        ["matched pairs", sum(map(lambda x: sum(library[x]["pair_match"]),
            library))],
        ["new alleles",  sum(map(lambda x: x[3], tables["new"]))],
        ["new unique alleles", sum(map(lambda x:
            len(allele_table(library[x]["new"], minimum)), library))],
        ["no start", sum(map(lambda x: x[3], tables["nostart"]))],
        ["no end", sum(map(lambda x: x[3], tables["noend"]))],
        ["unrecognised reads", unrecognised]
    ]

    return tables


def make_report(tables, handle):
    """
    Make an overview of the results.

    :arg tables: A nested dictionary containing overview tables.
    :type tables: dict
    :arg handle: Open writable handle to the report file.
    :type handle: stream
    """
    write_table(tables["summary"], "", handle)
    handle.write("\n")

    write_table(tables["library"], headers["markers"], handle)

    for i in tables["allele"]:
        handle.write("\nknown alleles for marker %s:\n" % i)
        write_table(tables["allele"][i]["known"], headers["allele"], handle)

        mean_length = 0
        sum_of_lengths = sum(map(lambda x: len(x[0]) * x[1],
            tables["allele"][i]["new"]))
        number_of_alleles = sum(map(lambda x: x[1],
            tables["allele"][i]["new"]))
        if number_of_alleles:
            mean_length = sum_of_lengths / number_of_alleles

        handle.write("\nnew alleles for marker %s (mean length %i):\n" % (i,
            mean_length))
        write_table(tables["allele"][i]["new"], headers["allele"], handle)


def write_files(tables, files):
    """
    Write the overview tables to the appropriate files.

    :arg tables: A nested dictionary containing overview tables.
    :type tables: dict
    :arg files: Nested dictionary containing writable file handles.
    :type files: dict
    """
    write_table(tables["summary"], "", files["summary"])
    write_table(tables["library"], headers["markers"], files["markers"])
    write_table(tables["known"], headers["overview"], files["known"])
    write_table(tables["new"], headers["overview"], files["new"])
    write_table(tables["nostart"], headers["nostartend"], files["nostart"])
    write_table(tables["noend"], headers["nostartend"], files["noend"])

    for i in tables["allele"]:
        write_table(tables["allele"][i]["known"], headers["allele"],
            files[i]["knownalleles"])
        write_table(tables["allele"][i]["new"], headers["allele"],
            files[i]["newalleles"])


def tssv(fasta_handle, library_handle, report_handle, path, threshold,
        minimum, is_fastq, indel_score):
    """
    Do the short structural variation analysis.

    :arg fasta_handle: Open readable handle to a FASTA file.
    :type fasta_handle: stream
    :arg library_handle: Open readable handle to a library file.
    :type library_handle: stream
    :arg report_handle: Open writable handle to the report file.
    :type report_handle: stream
    :arg path: Name of the output folder.
    :type path: str
    :arg threshold: Number of allowed mismatches per nucleotide.
    :type threshold: float
    :arg minimum: Minimum count per allele.
    :type minimum: int
    :arg is_fastq: Read FASTQ file instead of FASTA.
    :type is_fastq: bool
    :arg indel_score: Penalty score for insertions and deletions per
    nucleotide
    :type indel_score: int
    """
    total = 0
    unrecognised = 0
    library = parse_library(library_handle, threshold)

    if path:
        files = open_files(path, library)

    for record in SeqIO.parse(fasta_handle, "fastq" if is_fastq else "fasta"):
        ref = [str(record.seq), Seq.reverse_complement(str(record.seq))]
        ref_up = map(str.upper, ref)
        total += 1
        unknown = True

        for i in library:
            # Align against all-uppercase reference sequence.
            alignments = (
                align_pair(
                    ref_up[0], ref_up[1], library[i]["flanks"], indel_score),
                align_pair(
                    ref_up[1], ref_up[0], library[i]["flanks"], indel_score))
            matches = [False, False, False, False]
            classification = ""

            if alignments[0][0][0] <= library[i]["thresholds"][0]:
                cutout = ref[0][
                    max(0, alignments[0][0][1]-len(library[i]["flanks"][0])):
                    alignments[0][0][1]]
                if cutout.lower() != cutout:
                    library[i]["counts"][0] += 1
                    classification = "noend"
                    matches[0] = True
            if alignments[0][1][0] <= library[i]["thresholds"][1]:
                cutout = ref[0][
                    alignments[0][1][1]:
                    alignments[0][1][1]+len(library[i]["flanks"][1])]
                if cutout.lower() != cutout:
                    library[i]["counts"][2] += 1
                    classification = "nostart"
                    matches[1] = True
            if alignments[1][0][0] <= library[i]["thresholds"][0]:
                cutout = ref[1][
                    max(0, alignments[1][0][1]-len(library[i]["flanks"][0])):
                    alignments[1][0][1]]
                if cutout.lower() != cutout:
                    library[i]["counts"][1] += 1
                    classification = "noend"
                    matches[2] = True
            if alignments[1][1][0] <= library[i]["thresholds"][1]:
                cutout = ref[1][
                    alignments[1][1][1]:
                    alignments[1][1][1]+len(library[i]["flanks"][1])]
                if cutout.lower() != cutout:
                    library[i]["counts"][3] += 1
                    classification = "nostart"
                    matches[3] = True

            if (matches[0] and matches[1]) or (matches[2] and matches[3]):
                hit = int(matches[2] and matches[3])

                library[i]["pair_match"][hit] += 1
                pat = ref_up[hit][alignments[hit][0][1]:alignments[hit][1][1]]

                classification = "new"
                if library[i]["reg_exp"].match(pat):
                    classification = "known"

                library[i][classification][pat][hit] += 1

            if classification:
                unknown = False

                if path:
                    SeqIO.write([record], files[i][classification], "fasta")

        if unknown:
            unrecognised += 1

            if path:
                SeqIO.write([record], files["unknown"], "fasta")

    tables = make_tables(total, unrecognised, library, minimum)

    # Make the known alleles more human readable.
    for i in tables["allele"]:
        for j in tables["allele"][i]["known"]:
            j[0] = rewrite(library[i]["reg_exp"], j[0])
    for i in tables["known"]:
        i[4] = rewrite(library[i[0]]["reg_exp"], i[4])

    if path:
        write_files(tables, files)

    make_report(tables, report_handle)


def main():
    """
    Main entry point.
    """
    parser = argparse.ArgumentParser(description=usage[0], epilog=usage[1],
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("fasta", metavar="FASTA", type=argparse.FileType("r"),
        help="a FASTA file")
    parser.add_argument("lib", metavar="LIBRARY", type=argparse.FileType("r"),
        help="library of flanking sequences")
    parser.add_argument("-q", dest="fastq", action="store_true",
        help="if specified, treat input as FASTQ instead of FASTA")
    parser.add_argument("-m", dest="mismatches", type=float, default=0.08,
        help="mismatches per nucleotide (default=%(default)s)")
    parser.add_argument('-n', '--indel-score', type=int, default=1,
        help='insertions and deletions are penalised this number of times '
             'more heavily than mismatches (default=%(default)s)')
    parser.add_argument("-r", dest="report", type=argparse.FileType("w"),
        default=sys.stdout, help="name of the report file")
    parser.add_argument("-d", dest="path", type=str, help="output directory")
    parser.add_argument("-a", dest="minimum", type=int, default=0,
        help="minimum count per allele (default=%(default)s)")
    parser.add_argument('-v', action="version", version=version(parser.prog))

    args = parser.parse_args()

    try:
        tssv(args.fasta, args.lib, args.report, args.path, args.mismatches,
            args.minimum, args.fastq, args.indel_score)
    except OSError as error:
        parser.error(error)


if __name__ == "__main__":
    main()
