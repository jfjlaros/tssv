import argparse
import os
import re
import operator

def profile(sequence, name, files, markerlib, tellers, mismatches, report, newallel, error_start, error_end, overlapping):
	"""
	This profile routine does gathers information about read's best alignment score and prints that 

	@arg sequence: the read to be analyzed
	@type sequence: string
	@arg name: read identifier
	@type name: string
	@arg files: library of files to write to
	@type files: dict
	@arg markerlib: library of markers and regular expressions
	@type markerlib: dict
	@arg tellers: counters for read specifications
	@type tellers: array
	@arg mismatches: allowed mismatches per 25 nucleotides in alignment
	@type mismatches: int
	@arg report: dictionary to store report output
	@type report: dict
	@arg newallel: dictionary to store new alleles output
	@type newallel: dict
	@arg error_start: dictionary to store read information with no start
	@type error_start: dict
	@arg error_start: dictionary to store read information with no end
	@type error_start: dict
	@arg overlapping: dictionary to store read information with overlapping markers
	@type overlapping: dict

	"""
	score = {}
	for structure, list in markerlib.iteritems():
		if ((len(sequence) > len(list[0])) & (len(sequence) > len(list[1]))):
			align_all("", sequence, structure, list, score)
			revcompl = complement(sequence)
			align_all("RC", revcompl, structure, list, score)
			select_orientation(score, "begin", structure)
			select_orientation(score, "end", structure)
		else:
			tellers[0] += 1	# too short
			return
	ref_start, smallest_start = sort_smallest(score, "begin", "")
	start = re.sub("begin$", '', ref_start)
	ref_end, smallest_end = sort_smallest(score, "end", ref_start)
	end = re.sub("end$", '', ref_end)

	temporal = sequence
	if score[ref_start][2] == "rev":
		temporal = revcompl
	if (smallest_start <= (len(markerlib[start][0])/25.000000*mismatches[0])) and (smallest_end <= (len(markerlib[end][1])/25.000000*mismatches[0])):
		if start == end: 	# same markers in start and end
			#if the same start and end marker is found
			if (score[ref_start][2] == score[ref_end][2]):	# same orientation
				if (len(temporal)-score[ref_start][1]-score[ref_end][1] > 0):
					begin_allel = score[ref_start][1]
					end_allel = len(temporal) - score[ref_end][1]
					al_repeat = temporal[begin_allel:end_allel]
					reg_found, allel = patternMatching(markerlib[start][2], al_repeat)
					if (reg_found == 0):
						count_orientation_allel(start, ref_start, al_repeat, newallel, score);
						tellers[1] += 1	# new allele counter
						files[start][2].write(">start:"+start+" "+str(smallest_start)+" end:"+end+" "+str(smallest_end)+" name:"+name+" orientation: "+score[ref_start][2]+"\n"+temporal+"\n")
					else:
						count_orientation_allel(start, ref_start, allel, report, score);	
						tellers[8] += 1	# counter found alleles					
						files[start][3].write(">start:"+start+" "+str(smallest_start)+" end:"+end+" "+str(smallest_end)+" name:"+name+" orientation: "+score[ref_start][2]+"\n"+temporal+"\n")
				else:
					tellers[9] += 1	# counter for overlapping markers				
					count_orientation(ref_start, overlapping, score)
			else:
				tellers[2] += 1	# different orientations counter
				files[start][5].write(">start:"+start+" "+str(smallest_start)+" end:"+end+" "+str(smallest_end)+" name:"+name+"\n"+temporal+"\n")
		else:
			tellers[3] += 1	# counter for different markers in start and end
			files["All"][2].write(">start:"+start+" "+str(smallest_start)+" end:"+end+" "+str(smallest_end)+" name:"+name+" orientation: "+score[ref_start][2]+"\n"+temporal+"\n")
	elif (smallest_start <= (len(markerlib[start][0])/25.000000*mismatches[0])): #if no end is found
		tellers[4] += 1	#total reads with no end marker counter
		count_orientation (ref_start, error_end, score);
		files[start][0].write(">start:"+start+" "+str(smallest_start)+" end:"+end+" "+str(smallest_end)+" name:"+name+" orientation: "+score[ref_start][2]+"\n"+temporal+"\n")
	elif (smallest_end <= (len(markerlib[end][1])/25.000000*mismatches[0])):
		tellers[5] += 1	#count total reads with no start marker
		count_orientation (ref_end, error_start, score);
		files[end][1].write(">start:"+start+" "+str(smallest_start)+" end:"+end+" "+str(smallest_end)+" name:"+name+" orientation: "+score[ref_start][2]+"\n"+temporal+"\n")
	else:
		tellers[6] += 1	#count where no marker is found 
		files["All"][1].write(">start:"+start+" "+str(smallest_start)+" end:"+end+" "+str(smallest_end)+" name:"+name+" orientation: "+score[ref_start][2]+"\n"+temporal+"\n")

def count_orientation(ref, dictionary, score):
	"""
	This subroutine counts the observed amount of forward and reverse orientations of the markers found in reads

	@arg ref: The marker in the read
	@type ref: string
	@arg dictionary: The dictionary in which the orientations are counted
	@type dictionary: dict
	@arg score: The scoring library with the orientation information
	@type score: dict
	
	@return score: The dictionary with the alignment information
	@rtype: dict
	@return dictionary: The corresponding dictionary of the partly report
	@rtype dictionary: dict
	"""
	if ref in dictionary:
		dictionary[ref][0] += 1
	else:
		dictionary[ref] = [1, 0, 0]
	if score[ref][2] == "for":
		dictionary[ref][1] += 1
	else:
		dictionary[ref][2] += 1

def count_orientation_allel(ref1, ref2,al_repeat, dict, score):
	"""
	This subroutine counts the observed amount of forward and reverse orientations of the markers found in reads where an allele is found
	
	@arg ref1: The marker that is found in the read
	@type ref1: string
	@arg ref2: The name of the start marker
	@type ref2: string
	@arg al_repeat: The allele found in the read
	@type al_repeat: string
	@arg dictionary: The dictionary in which the orientations are counted
	@type dictionary: dict
	@arg score: The scoring library with the orientation information
	@type score: dict

	@return score: The dictionary with the alignment information
	@rtype: dict
	@return dictionary: The corresponding dictionary of the partly report
	@rtype dictionary: dict
	"""

	if ref1 in dict:
		dict2 = dict[ref1]
	else:
		dict[ref1] = {}
		dict2 = {}
	if al_repeat in dict2:
		dict2[al_repeat][0] +=1
	else:
		dict2[al_repeat] = [1, 0, 0]	# store good alleles with orientation and structure
	if score[ref2][2] == "for":
		dict2[al_repeat][1] += 1
	else:
		dict2[al_repeat][2] += 1
	dict[ref1] = dict2
	
def ReadArgs():
	"""
	This routine does the argument parsing

	@arg argument: fasta file, library, allowed mismatches, output directory
	@type argument: strings or integer
	"""
	result = []
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', nargs = 1, help='A fasta file with the reads', metavar='Fasta file name', type=file, required = True, dest='fasta')
	parser.add_argument('-l', nargs = 1, help='The library with the flanking sequences', metavar='Marker library file name', type=file, required = True, dest='library')
	parser.add_argument('-m', nargs = 1, help='allowed mismatches for 25nt (recommended=2)', metavar='mismatches', type=int, default = [2], dest = "mismatches")
	parser.add_argument('-d', nargs = 1, help='Output directory', metavar='directory name', required = True, dest = 'path')
	args = parser.parse_args()
	return args

def patternMatching(regular, al_repeat):
	"""
	This routine does the pattern matching

	@arg regular: The regular expression used for matching
	@type regular: string
	@arg al_repeat: The allele
	@type al_repeat: string
	
	@return reg_found: boolean to determine whether an allele was matched
	@rtype reg_found: bool
	@return allele: The pattern matched allele
	@rtype: string
	"""
	reg_found = 0
	allel = ""
	p = re.compile(regular)
	for m in p.findall(al_repeat):
		reg_found = 1
		if (not isinstance(m, str)):
			for i in range(len(m)):
				count = 0
				while (re.match(m[i], al_repeat)):
					if (m[i] == ''):
						break
					al_repeat = re.sub("^"+m[i], '', al_repeat)
					count += 1
				if (count != 0):				
					allel += m[i] + "x" + str(count) + "\t"
		else:
			count = 0
			while (re.match(m, al_repeat)):
				al_repeat = re.sub("^"+m, '', al_repeat)
				count += 1
			allel += m + "x" + str(count) + "\t"
	allel = re.sub("\t$", '', allel)
	return(reg_found, allel)

def ReadMarkers(FileName):
	"""
	This routine reads the marker file and places it in the total dictionary
	
	@arg FileName: A marker file
	@type FileName: file

	@return total: The readed in marker library
	@rtype total: dict
	"""
	total = {}
	for line in FileName:
		marker = line.split("\t")
		naam = marker[0]
		marker[0] = marker[1]
		marker[1] = marker[2]
		temp = marker[3].split(" ")
		del marker[3]
		RegularExpression = '^'
		for index in range(len(temp)):
			if (index % 3 == 0):
				RegularExpression += "(" + temp[index] + ")"
			elif (index % 3 == 1):
				RegularExpression += "{" + temp[index] + ","
			elif (index % 3 == 2):
				RegularExpression += temp[index].rstrip() + "}"
		marker[2] = RegularExpression + "$"
		total[naam] = marker
	FileName.close()
	return total

def OpenOutputFiles(path, markerlib):
	"""
	This routine opens the files for writing output

	@arg path: The output directory
	@type path: string
	@arg markerlib: The library containing the markers and the regular expressions
	@type markerlib: dict

	@return files: The opened file names
	@rtype files: dict
	"""
	try: 
		os.mkdir(path)
	except:
		#print "dir {0} bestaat al".format(path)
		pass
	files = {}
	file = []
	file.append(open(path + '/report.txt', 'w'))
	file.append(open(path + '/NoFoundMarker.txt', 'w'))
	file.append(open(path + '/DifferentStructures.txt', 'w'))
	key = "All"
	files[key] = file
	for key in markerlib:
		try: 
			os.mkdir(path + "/" + key)
		except:
			#print "dir {0} bestaat al".format(path + "/" + "key")
			pass
		file = []
		file.append(open(path + "/" + key + "/NoEnd.txt", 'w'))
		file.append(open(path + "/" + key + "/NoBeginning.txt", 'w'))
		file.append(open(path + "/" + key + "/NewAllele.txt", 'w'))
		file.append(open(path + "/" + key + "/GoodAllele.txt", 'w'))
		file.append(open(path + "/" + key + "/DifferentOrientation.txt", 'w'))
		files[key] = file
	return files

def ProcessFasta(FileName, files, markerlib, mismatches):
	"""
	This routine reads through the fasta file calling on profile subroutine which handles the rest of the work

	@arg FileName: The fasta file
	@type FileName: file
	@arg files: dictionary that contains The file names to write to
	@type files: dict
	@arg markerlib: The library containing the markers and the regular expressions
	@type markerlib: dict
	@arg mismatches: The amount of allowed mismatches
	@type mismatches: list
	"""
	tellers = [0,0,0,0,0,0,0,0,0,0]	
							# 0 = counter "too short"
							# 1 = new allele counter
							# 2 = different orientations counter
							# 3 = counter for different markers in start and end
							# 4 = total reads with no end marker counter
							# 5 = count total reads with no start marker
							# 6 = count where no marker is found 
							# 7 = reads counter
							# 8 = counter found alleles
	report      = {}
	newallel    = {}
	error_start = {}
	error_end   = {}
	overlapping = {}
	sequence    = None
	name        = None
	for line in FileName: 
		line = line.rstrip('\r\n')
		match = re.match(r"^>", line)
		if (match != None):
			tellers[7] += 1
			name = line.strip('>')
			if (sequence != None):
				profile(sequence, name, files, markerlib, tellers, mismatches, report, newallel, error_start, error_end, overlapping)
				sequence = None
		else:
			if (sequence == None):
				sequence = line
			else:
				sequence += line
	profile(sequence, name, files, markerlib, tellers, mismatches, report, newallel, error_start, error_end, overlapping)
	printOutput(files, tellers, report, newallel, error_start, error_end, overlapping)


def makeMatrix(xSize, ySize):
	"""
	This routine initializes a matrix. The first row are numbers from 0 to matrix size. The rest are zero's

	@arg xSize: The size of the x component of the matrix
	@type xSize: int
	@arg xSize: The size of the y component of ch matrix
	@type ySize: int

	@return matrix: The alignment matrix
	@rtype: two dimensional array
	"""
	matrix = []
	for i in range(xSize):
		matrix.append([])
		for j in range(ySize):
			if (i == 0):
				matrix[i].append(j)
			else: 
				matrix[i].append(0)
	return matrix

def printMatrix(matrix, xSize, ySize, seq1, seq2):
	"""
	This subroutine is for debugging purposes and can print the alignment matrix
	
	@arg matrix: The matrix
	@type matrix: dict
	@arg xSize: The size of the x component of the matrix
	@type xSize: int
	@arg xSize: The size of the y component of ch matrix
	@type ySize: int
	@arg seq1: The sequence to be aligned to
	@type seq1: string
	@arg seq2: The sequence to be aligned
	@type seq2: string
	"""
	txt = "  - "
	for y in range(ySize-1):# Print the second sequence.
		txt += seq2[y] + " "
	print txt
	for x in range(xSize):
		if x == 0: txt = "-" + " "
		else: txt = seq1[x-1] + " "
		for y in range(ySize):
			txt += str(matrix[x][y]) + " "         # Print the contents of the matrix.
		print txt

def align(matrix, xSize, ySize, seq1, seq2):
	"""
	This routine fills the alignment matrix

	@arg xSize: The size of the x component of ch matrix
	@type xSize: int
	@arg xSize: The size of the y component of the matrix
	@type ySize: int
	@arg seq1: The sequence to be aligned to
	@type seq1: string
	@arg seq2: The sequence to be aligned
	@type seq2: string

	@return matrix: The alignment matrix
	@rtype: two dimensional array
	"""
	for x in range(1, xSize):
		for y in range(1, ySize):
			matrix[x][y] = min(matrix[x -1][y] + 1, matrix[x][y - 1] + 1, matrix[x - 1][y - 1] + int(seq1[x - 1] != seq2[y - 1]))
	return matrix

def findMin(matrix, xSize, ySize):
	"""	
	Find the minimum distance, ignoring a trailing gap in the sequence associated
	with the number of rows in an alignment matrix. If the minimum distance is
	found, also return the row number.
	
	It is assumed that the number of rows is larger than the number of columns.
	
	@arg matrix: An xSize * ySize matrix.
	@arg xSize: The size of the x component of the matrix
	@type xSize: int
	@arg xSize: The size of the y component of the matrix
	@type ySize: int

	@return minimum: The minimum distance.
	@rtype minimum: int
	@return xMin: Row number
	@rtype xMin: int
	"""
	minimum = ySize - 1
	xMin = 0
	for x in range(xSize):
		if (matrix[x][ySize - 1] < minimum):
			minimum = matrix[x][ySize -1]
			xMin = x
	return (minimum, xMin)

def al(seq1, seq2):
	"""
	This routine does the calling on the alignment subroutines
	
	@arg seq1: The sequence to be aligned to
	@type seq1: string
	@arg seq2: The sequence to be aligned
	@type seq2: string

	@return minimum: The minimum distance.
	@rtype minimum: int
	@return xMin: Row number
	@rtype xMin: int
	"""
	xSize = len(seq1) + 1
	ySize = len(seq2) + 1
	matrix = makeMatrix(xSize, ySize)
	align(matrix, xSize, ySize, seq1, seq2)
	minimum, xMin = findMin(matrix, xSize, ySize)
	return (minimum, xMin)
		
def align_all(label, sequence, structure, list, score):
	"""
	This routine calls the alignment subroutines and stores the alignment output information in the score dictionary

	@arg label: Used to distinguish between first and second marker alignment
	@arg type: string
	@arg sequence: The to be aligned read
	@type sequence: string
	@arg structure: The marker entry name
	@type structure: string
	@arg list: The list containing flanking sequences and the regular expression
	@type list: list
	@arg score: The dictionary containing the alignment scores
	@type score: dict 
	
	@return score: The dictionary containing the alignment scores
	@rtype score: dict
	"""
	a, b = al(sequence, list[0])
	score[structure+"begin"+label] = [a, b]
	rev_sequence = sequence[::-1]
	rev_ref = list[1][::-1]
	a, b = al(rev_sequence, rev_ref)
	score[structure+"end"+label] = [a, b]

def complement(sequence):
	"""
	This subroutine determines the reverse complement of a sequence
	
	@arg sequence: The sequence of which the reverse complement is generated
	@type sequence: string

	@returns: reverse complement of a DNA string
	@rtype: string
	"""
	seq = list(sequence)
	temp = ''
	for n in seq:
		if   (n == "A"): temp += "T"
		elif (n == "T"): temp += "A"
		elif (n == "G"): temp += "C"
		elif (n == "C"): temp += "G"
	temp = temp[::-1]
	return temp

def select_orientation(score, label, structure):
	"""
	This subroutine determines whether the normal or the Reverse Complement alignment was best
	
	@arg score: The dictionary containing the alignment information
	@type score: dict
	@arg label: can be beginning or end and is used to distinguish between the first or second aligned structure
	@type label: string
	@arg structure: The marker name
	@type structure: string

	@returns: score dictionary
	@rtype: dict
	"""
	a = score[structure+label]
	b = score[structure+label+"RC"]
	if a[0] > b[0]:
		a[0] = b[0]
		a.append("rev")
		a[1] = b[1]
	else:
		a.append("for")
	del score[structure+label+"RC"]


def sort_smallest(score, label, refstart):
	"""
	This subroutine selects the smallest marker out of the alignment. 
	If multiple smallest end marker are found it stores the one which corresponds to the begin marker
	
	@arg score: The scoring dictionary containing the alignment information
	@type score: dict
	@arg label: To distinguish between first and second marker alignment
	@type label: string

	@return ref: The smallest marker name
	@rtype: string
	@return smallest: The smallest score out of all markers
	@rtype smallest: int
	"""
	smallest = -1
	ref = ""
	for reference, list in score.iteritems():
		if smallest != -1:
			if (list[0]  < smallest) & (re.search(label, reference) != None):
				ref = reference
				smallest = list[0]
			if (label == "end") & (list[0] == smallest) & (reference == re.sub("begin$", "end", refstart)):
				ref = reference
				smallest = list[0]
		else:
			if (re.search(label, reference) != None):
				smallest = list[0]
				ref = reference
	return(ref, smallest)

def printScore(score):
	"""
	Prints the score dictionary

	@arg score: The dictionary with the alignment information
	@type score: dict
	"""
	print "\nTable score:"
	for reference, list in score.iteritems():
		print '{0:15s} {1:15s}'.format(reference, list)
	print "\n"

def print_double(report, reportfile):	#this subroutine does the printing for the report and the new hashes.
	"""
	This routine can print a two dimensional dictionary
	
	@arg report: The report to be printed
	@type report: dict
	@arg reportfile: The file to be printed to
	@type reportfile: file
	"""
	for ref, dict in sorted(report.iteritems()):
		for regexp, list in sorted(dict.iteritems(), key = operator.itemgetter(1), reverse = True):
			reportfile.write(ref + "\tFor:" + str(list[1]) + "\tRev:" + str(list[2]) + "\t" + str(list[0]) + "\t"+ regexp + "\n")

def print_single(report, label, reportfile):
	"""
	This subroutine prints a one dimensional dictionary
	
	@arg report: The report to be printed
	@type report: dict
	@arg label: Label to be printed to distinguish which report is printed
	@type label: string
	@arg reportfile: The file to be printed to
	@type reportfile: file
	"""
	for ref, list in sorted(report.iteritems(), key = operator.itemgetter(1), reverse = True):
		ref = re.sub("begin$", '', ref)
		ref = re.sub("end$", '', ref)
		reportfile.write(label + ": " + ref + "\tFor:" + str(list[1]) + "\tRev:" + str(list[2]) + "\t" + str(list[0]) + "\n")

def printOutput(files, tellers, report, newallel, error_start, error_end, overlapping):
	"""
	This routine prints the overall report
	
	@arg files: a dictionary containing the output files
	@type files: dict
	@arg tellers: counters for all the possible features of the reads
	@type tellers: array
	@arg report: The dictionary containing the found alleles that matched with the regular expression
	@type report: dict
	@arg newallel: The dictionary containing the found alleles that didn't match with the regular expression
	@type newallel: dict
	@arg error_start: The dictionary containing the reads were no start marker was found
	@type error_start: dict
	@arg error_end: The dictionary containing the reads were no end marker was found
	@type error_end: dict
	@arg overlapping: The dictionary containing the reads were markers overlap
	@type overlapping: dict

	"""
	file = files["All"][0]
	file.write("total reads:              "+str(tellers[7])+"\n")#now all variables and hashes can be printed
	file.write("beginning but no end:     "+str(tellers[4])+"\n")
	file.write("end but no beginning:     "+str(tellers[5])+"\n")
	file.write("no beginning no end:      "+str(tellers[6])+"\n")
	file.write("different structures:     "+str(tellers[3])+"\n")
	file.write("different orientations:   "+str(tellers[2])+"\n")
	file.write("new alleles:              "+str(tellers[1])+"\n")
	file.write("good alleles:             "+str(tellers[8])+"\n")
	file.write("too short input sequences:"+str(tellers[0])+"\n")
	file.write("overlapping markers:      "+str(tellers[9])+"\n")
	file.write("\n################report\n")
	print_double(report, files["All"][0])
	file.write("\n################error\n")
	print_single(error_end,   "no end"  , files["All"][0])
	print_single(error_start, "no start", files["All"][0])
	print_single(overlapping, "overlapping", files["All"][0])
	file.write("\n################new\n")
	print_double(newallel, files["All"][0])
		
def main():
	"""
	This subroutine calls to open output files and calls to read the markerlib
	"""
	files = {}
	args = ReadArgs()
	markerlib = ReadMarkers(args.library[0])
	files = OpenOutputFiles(args.path[0], markerlib)
	ProcessFasta(args.fasta[0], files, markerlib, args.mismatches)

if __name__ == "__main__":
	"""
	Main entry point
	"""
	main()
