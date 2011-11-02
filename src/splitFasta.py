#!/usr/bin/python

from Bio import SeqIO
import sys

def addToLibrary(library, line) :
    """
    """

    library.append(line.split())
    library[-1].append(open(library[-1][0] + ".txt", "w"))
    library[-1].append(open(library[-1][0] + "_counted.txt", "w"))
    library[-1].append(0)
    library[-1].append({})
#addToLibrary

def writeOutput(entry, description, sequence) :
    """
    """

    entry[3].write(">%s\n" % description)
    entry[3].write("%s\n" % sequence)
    entry[5] += 1
    if entry[6].has_key(sequence) :
        entry[6][sequence][1] += 1
    else :
        entry[6][sequence] = [description, 1]
#writeOutput

# Read the library file.
handle = open(sys.argv[2], "r")
library = []
line = handle.readline()
while line :
    addToLibrary(library, line)
    line = handle.readline()
#while
handle.close()
addToLibrary(library, "Unrecognised _INVALID_ _INVALID_")

# Split the fasta file.
total = 0
recognised = 0
handle = open(sys.argv[1], "r")
for record in SeqIO.parse(handle, "fasta") :
    hit = False
    for entry in library :
        if entry[1] in record.seq :
            flanks = str(record.seq).split(entry[1])
            writeOutput(entry, record.description, flanks[0] + entry[2] + flanks[1])
            recognised += 1
            hit = True
        #if
        if entry[1] == "_INVALID_" and not hit :
            writeOutput(entry, record.description, str(record.seq))
    #for
    total += 1
#for
handle.close()

print "Total: %i" % total
for entry in library :
    print "%s: %i" % (entry[0], entry[5])
    unique = 0
    maximum = 0
    for i in sorted(entry[6], key = lambda item: -entry[6][item][1]) :
        entry[4].write(">COUNT %04d %s\n" % (entry[6][i][1], entry[6][i][0]))
        entry[4].write("%s\n" % i)
        if entry[6][i][1] > maximum :
            maximum = entry[6][i][1]
        unique += 1
    #for
    print "  unique: %i" % unique
    print "  max: %i" % maximum
    entry[3].close()
    entry[4].close()
#for
print "Unrecognised: %i" % (total - recognised)
