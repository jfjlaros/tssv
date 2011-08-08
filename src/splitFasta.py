#!/usr/bin/python

from Bio import SeqIO
import sys

# Read the library file.
handle = open(sys.argv[2], "r")
library = []
line = handle.readline()
while line :
    library.append(line.split())
    library[-1].append(open(library[-1][0] + ".txt", "w"))
    library[-1].append(open(library[-1][0] + "_counted.txt", "w"))
    library[-1].append(0)
    library[-1].append({})
    line = handle.readline()
#while
handle.close()

# Split the fasta file.
total = 0
recognised = 0
handle = open(sys.argv[1], "r")
for record in SeqIO.parse(handle, "fasta") :
    for entry in library :
        if entry[1] in record.seq :
            flanks = str(record.seq).split(entry[1])
            newSeq = flanks[0] + entry[2] + flanks[1]
            entry[3].write(">%s\n" % record.description)
            entry[3].write("%s\n" % newSeq)
            entry[5] += 1
            if entry[6].has_key(newSeq) :
                entry[6][newSeq][1] += 1
            else :
                entry[6][newSeq] = [record.description, 1]
            recognised += 1
        #if
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
