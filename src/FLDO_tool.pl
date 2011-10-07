#!/usr/local/bin/perl -w


=head1 NAME
	FLDO_tool.pl a module that extracts STR sequences from a fasta file and analyses them.

=head1 SYNOPSIS
	perl FLDO_tool.pl -FLDO reads file- -marker library- -allowed mismatches for 25nt (recommended=2)- -run code-

=head1 DESCRIPTION
	NGProf analyses reads in fasta format created by a sequencer. Reads can contain markers that surround a STR.
	STRs are extracted and analysed by a regular expression.

=head2 METHODS
	The first step is
	the identification of the specific STR from the multiplex mixture based on a
	reference library. The reference library consists of a marker name tag, the
	upstream flanking sequence adjacent to the STR with a recommended amount of
	25 basepairs, the downstream flanking sequence and the repeat structure
	written as a regular expression, an example of the library is:
	
	Name \t start_sequence \t end_sequence \t AATA \s 0 \s 15 \s ACTT \s 3 \s 10

	The analysis to be performed with this
	library is automated for forward and reverse reads by also interpreting the sequences 
	as its reverse complement. For each read four semi-global alignments 
	are done with each pair of the flanking sequences from the library. 
	The aligner determines the position of the first occurrence with
	the minimum edit distance of all substrings from
	the read and STR marker sequence. The aligner gives the right most coordinate
	on the read pointing to the start of the STR. The starting position of the
	second marker is found in the same way by reversing both the marker and the
	read. After the minimum edit distance is determined for all markers in the
	library, the minimum number is calculated for up- and downstream flanking
	sequences separately. The maximum number of mismatches, insertions and
	deletions allowed can be set in this analysis and is used to reject marker
	pairs with too many differences. In general we recommended 2 mismatches per
	25 nucleotides. A report per marker is made about the number of reads where
	one or two of the flanking markers is not recognized. After forward and reverse
	marker positions are identified the STR can be extracted from the read. STRs
	are further analyzed using the given regular expression. NGProf outputs the
	repetitive blocks of sequences clarifying the (complex) structure. STR's that
	don't match the regular expression are counted and returned as potential new
	alleles including their sequence and orientation, so they can be used to identify 
	previously unknown repeat structures. To identify SNPs and indels 
	in marker sequences an extra tool is available in the download section.
	This tool uses the now known repeat structure plus a few basepairs of marker 
	sequence to do a search in the fasta file with a regular expression. 
	Sequences in the file are returned including their marker sequence so markers
	can be searched for SNPs and indels.

=head1 LICENSE

This software is available under GNU public license

=head1 AUTHOR

J.W.F. van der Heijden, K.J. van der Graag, J.F.J. Laros.

=cut


#use strict;
use List::Util qw[min max];

if(scalar(@ARGV) != 4){#check for the amount of parameters in the command code
	print "Please use the correct parameters\nUsage: FLDO_tool.pl -FLDO reads file- -marker library- -allowed mismatches for 25nt (recommended=2)- -run code-\n";
	exit();
}

#Example library: Name \t start_sequence \t end_sequence \t AATA \s 0 \s 15 \s ACTT \s 3 \s 10
open (REF, $ARGV[1]) || die "cannot open library\n";#opens the library with markers and repeat patterns
my %hash=();#the hash in which marker sequences are stored
my %regular=();#the hash in which the regular expression is stored. 
while (<REF>){
	if ($_ eq ""){next;}
	$_ =~ s/\n|\r//g;
	my @element=split(/\t/,$_);
	$hash{$element[0]}{ref1}=$element[1];#places start marker in hash under ref1
	$hash{$element[0]}{ref2}=$element[2];#places end marker in hash under ref2
	for (my $i=3;$i<@element;$i++){
		@repeat_structure=split(/ /,$element[$i]);
		my $reg="";
		for (my $j=0;$j<@repeat_structure;$j+=3){
			$reg.="($repeat_structure[$j]){".$repeat_structure[($j+1)].",".$repeat_structure[($j+2)]."}";#builds regular expression for use in pattern matching
		}
		$regular{$element[0]}{$i}=$reg;#stores the regular expression in "regular"
	}
}
close REF;
my $read=0;#counts the total reads
my $sequence="";#current analysed sequence
my @al1="";#storage for edit distances
my @al2="";#storage for edit distances
my %report=();#storage for output
my %error=();
my %new=();
my $correct_allel=0;#variables for encountering possibilities in reads
my $different_structures=0;
my $different_orientation=0;
my $no_end=0;
my $no_start=0;
my $new_allel=0;
my $nothing=0;
my $too_short=0;
my $temporal="";#for switching between forward and reverse
mkdir "$ARGV[3]";
open (IN, $ARGV[0]) || die "cannot open reads file\n";#opens the reads file
open (OUT, ">$ARGV[3]/report.txt") || die "cannot open file for writing report\n";#opens the output file for writing report
open (OUT2, ">$ARGV[3]/NoFoundMarker.txt") || die "cannot open file for writing report\n";#opens the output file for writing report
open (OUT3, ">$ARGV[3]/DifferentStructures.txt") || die "cannot open file for writing report\n";#opens the output file for writing report
foreach $key (keys %hash){
	mkdir "$ARGV[3]/$key";
	open ($key."1", ">$ARGV[3]/$key/NoEnd.txt") || die "cannot open file for writing marker $key";
	open ($key."2", ">$ARGV[3]/$key/NoBeginning.txt") || die "cannot open file for writing marker $key";
	open ($key."3", ">$ARGV[3]/$key/NewAllele.txt") || die "cannot open file for writing marker $key";
	open ($key."4", ">$ARGV[3]/$key/GoodAllele.txt") || die "cannot open file for writing marker $key";
	open ($key."5", ">$ARGV[3]/$key/DifferentOrientation.txt") || die "cannot open file for writing marker $key";
}

#The while loop below reads through the fasta file.
$name="";
while (<IN>){
	$_ =~ s/\n|\r//g;
	if ($_ =~ m/^>/){
		$read++;#count reads
		&profile($sequence,$name);#the profiling subroutine
		$_ =~ s/>//;
		$name=$_;#read id
	} else {
		$sequence.=$_;#build sequence
	}
}
&profile($sequence,$name);#again the subroutine for last read

close IN;
print OUT "total reads: $read\n";#now all variables and hashes can be printed
print OUT "beginning but no end: $no_end\n";
print OUT "end but no beginning: $no_start\n";
print OUT "no beginning no end: $nothing\n";
print OUT "different structures: $different_structures\n";
print OUT "different orientations: $different_orientation\n";
print OUT "new allels: $new_allel\n";
print OUT "good allels: $correct_allel\n";
print OUT "too short input sequences: $too_short\n"; 
print OUT "\n################report\n";
&print_double(\%report);
print OUT "\n################error\n";
&print_single(\%error_end, "no end");
&print_single(\%error_start, "no start");
print OUT "\n################new\n";
&print_double(\%new);
close OUT;

sub profile {
	$sequence=$_[0];
	$name=$_[1];
	if ($sequence ne ""){
		my %score=();
		my $found=0;
		foreach my $structure (keys %hash){
			if (length ($sequence) > length ($hash{$structure}{ref1})){#it is important that the length of the input sequence is longer than the aligned structure	
				#&align_all($sequence, \%hash, $structure, \%score);			
				@al1=al($sequence,$hash{$structure}{ref1});#gives the edit distance and position on sequence of ref1
				$score{$structure}{begin}=$al1[0];
				$score{$structure}{startpos}=$al1[1];
				
				$rev_sequence=reverse($sequence);				
				$rev_ref=reverse($hash{$structure}{ref2});
				@al2=al($rev_sequence,$rev_ref);#gives the edit distance and position from end on sequence of ref2
				$score{$structure}{end}=$al2[0];
				$score{$structure}{endpos}=$al2[1];
				
				#The sequence is converted to the reverse complement
				$revcompl = complement($sequence);

				@al1=al($revcompl,$hash{$structure}{ref1});#gives the edit distance and position on rev_compl sequence of ref1
				$score{$structure}{beginRC}=$al1[0];
				$score{$structure}{startposRC}=$al1[1];
				
				$rev_sequence=reverse($revcompl);				
				$rev_ref=reverse($hash{$structure}{ref2});
				@al2=al($rev_sequence,$rev_ref);#gives the edit distance and position from end on rev_compl sequence of ref2
				$score{$structure}{endRC}=$al2[0];
				$score{$structure}{endposRC}=$al2[1];
				$found=1;
				#the smallest edit distance is selected between forward and reverse-complement sequence
				&select_orientation(\%score,"begin", "startpos", "orientation", $structure);
				&select_orientation(\%score,"end", "endpos", "orientationend", $structure);				
			}
		}
		if ($found==0){$too_short++;next;}
		
		(my $ref_start,my $smallest_start) = &sort_smallest(\%score,"begin");	
		(my $ref_end,my $smallest_end) = &sort_smallest(\%score,"end",$ref_start);
		
		if ($score{$ref_start}{orientation} eq "rev"){$temporal=$revcompl;}else{$temporal=$sequence;}
		#print $ref_end."\n".$score{$ref_end}{orientationend}."\n";
		if ($smallest_start <= (length ($hash{$ref_start}{ref1})/25*$ARGV[2]) && $smallest_end <= (length ($hash{$ref_end}{ref2})/25*$ARGV[2])){#if the edit distance is small enough for start and end
			if ($ref_start eq $ref_end){#and if the same start and end marker is found
				if ($score{$ref_start}{orientation} eq $score{$ref_end}{orientationend}){#in the same orientation
					$al_repeat=substr($temporal,$score{$ref_start}{startpos},(length($temporal)-$score{$ref_start}{endpos}-$score{$ref_start}{startpos}));#substract the repeat
					my $reg_found=0;
					foreach $teller (keys %{$regular{$ref_start}}){#foreach regular expression
						if ($al_repeat =~ m/^$regular{$ref_start}{$teller}$/i){#see if pattern is found
							$reg_found=1;
							$correct_allel++;
							my $allel="";
							my $pos=0;
							foreach $expr (1..$#-) {
								if (!defined ${$expr}){next;}
								$sub_slice=substr($al_repeat,$pos,($+[$expr]-$pos));#determine repeat structure
								$unit=length($sub_slice)/length(${$expr});
								$pos=$+[$expr];
								$allel.=${$expr}."x".$unit."\t";
							}
							chop $allel;
					#&count_orientation_allel(\%report, \%score, $ref_start, $allel);						
							$report{$ref_start}{$allel}{amount}++;#store good allels with orientation and structure
							if ($score{$ref_start}{orientation} eq "for"){
								$report{$ref_start}{$allel}{amountFor}++;
							} else {
								$report{$ref_start}{$allel}{amountRev}++;
							}


							$temp=$ref_start."4";
							print $temp ">start:$ref_start $smallest_start end:$ref_end $smallest_end name:$name orientation:$score{$ref_start}{orientation}\n$temporal\n";
						}
					}
					if ($reg_found==0){
						#print3dArray(\%new);
						#&count_orientation_allel(\%new, \%score, $ref_start, $al_repeat);
						$new{$ref_start}{$al_repeat}{amount}++;#store new allels with orientation
						if ($score{$ref_start}{orientation} eq "for"){
							$new{$ref_start}{$al_repeat}{amountFor}++;
						} else {
							$new{$ref_start}{$al_repeat}{amountRev}++;
						}
						$new_allel++;
						#print3dArray(\%new);
						$temp=$ref_start."3";
						print $temp ">start:$ref_start $smallest_start end:$ref_end $smallest_end name:$name orientation:$score{$ref_start}{orientation}\n$temporal\n";
					}
				} else {
					$different_orientation++;#count when different orientations are found
					$temp=$ref_start."5";
					#print $temp."\n";
					print $temp ">start:$ref_start $smallest_start end:$ref_end $smallest_end name:$name\n$sequence\n";
				}
			} else {
				$different_structures++;#count when different markers are best in start and end
				print OUT3 ">start:$ref_start $smallest_start end:$ref_end $smallest_end name:$name orientation:$score{$ref_start}{orientation}\n$temporal\n";
			}
		} elsif ($smallest_start <= (length ($hash{$ref_start}{ref1})/25*$ARGV[2])){#if no end is found
			$no_end++;#count total reads with no end marker
			&count_orientation (\%error_end, \%score, $ref_start);			
			$temp=$ref_start."1";
			print $temp ">start:$ref_start $smallest_start end:$ref_end $smallest_end name:$name orientation:$score{$ref_start}{orientation}\n$temporal\n";
		} elsif ($smallest_end <= (length ($hash{$ref_end}{ref2})/25*$ARGV[2])){
			$no_start++;#count total reads with no start marker
			&count_orientation (\%error_start, \%score, $ref_end);
			$temp=$ref_end."2";
			print $temp ">start:$ref_start $smallest_start end:$ref_end $smallest_end name:$name orientation:$score{$ref_end}{orientationend}\n$temporal\n";
		} else {
			$nothing++;#count where no marker is found 
			print OUT2 ">start:$ref_start $smallest_start end:$ref_end $smallest_end name:$name\n$temporal\n";
		}
		$sequence="";
	}
}

sub align_all{
	my $seq=shift;
	my %hash=%{(shift)};
	my $structure=shift;
	my %score=%{(shift)};
	#print "$seq, $hash{$structure}{ref1}\n";
	@al1=al($seq,$hash{$structure}{ref1});#gives the edit distance and position on sequence of ref1
	$score{$structure}{begin}=$al1[0];
	$score{$structure}{startpos}=$al1[1];
	my $rev_sequence=reverse($seq);				
	my $rev_ref=reverse($hash{$structure}{ref2});
	@al2=al($rev_sequence,$rev_ref);#gives the edit distance and position from end on sequence of ref2
	$score{$structure}{end}=$al2[0];
	$score{$structure}{endpos}=$al2[1];
}

sub print3dArray {
	my %A = %{(shift)};
        #my $i, $j, $k;

	foreach my $i (%A){
		foreach my $j ($A{$i}){
			foreach my $k ($A{$i}{$j}) {
				printf("%i ", $A{$i}{$j}{$k});
			}
			printf("\n");
		}
		printf("\n");
	}
}

sub count_orientation_allel{
	my %x=%{(shift)};
	my %y=%{(shift)};
	my $ref=shift;
	my $subject=shift;
	$x{$ref}{$subject}{amount}++;#store good allels with orientation and structure
	if ($y{$ref}{orientation} eq "for"){
		$x{$ref}{$subject}{amountFor}++;
	} else {
		$x{$ref}{$subject}{amountRev}++;
	}
}

sub count_orientation{
	my %x=%{(shift)};
	my %y=%{(shift)};
	my $ref=shift;
	$x{$ref}{amount}++;#store new allels with orientation
	if ($y{$ref}{orientation} eq "for"){
		$x{$ref}{amountFor}++;
	} else {
		$x{$ref}{amountRev}++;
	}
}
sub select_orientation{
	my %x=%{(shift)};
	my $label1=shift;
	my $label1RC=$label1."RC";
	my $label2=shift;
	my $label2RC=$label2."RC";
	my $label3=shift;
	
	my $y=shift;
        #printf("%i\n", $x{$y}{$label1RC});
	if ($x{$y}{$label1} > $x{$y}{$label1RC}){
		$x{$y}{$label1} = $x{$y}{$label1RC};
		$x{$y}{$label3}="rev";
		$x{$y}{$label2} = $x{$y}{$label2RC};
	} else {
		$x{$y}{$label3}="for";
	}
}

sub sort_smallest{#this subroutine sorts out the smallest number in a hash and gives the key with it.
	my %x=%{(shift)};
	my $label=shift;
	my $start=shift;
	my $smallest="";
	foreach my $reference (keys %x){
		if ($smallest ne ""){
			if ($x{$reference}{$label}<$smallest){
				$ref=$reference;
				$smallest=$x{$reference}{$label};
			}
			if ($label eq "end" && $x{$reference}{$label}==$smallest && $reference eq $start){
				$ref=$reference;
				$smallest=$x{$reference}{$label};
			}
		} else {
			$smallest=$x{$reference}{$label};
			$ref=$reference;
		}
	}
	return ($ref, $smallest);
}

sub complement {#this subroutine determines the reverse complement of the sequence
	my $x=shift;	
	@seq=split(//,$x);
	my $temp="";
	my $revcompl="";
	for ($k=0;$k<@seq;$k++){
		if ($seq[$k] eq "A"){$temp.= "T";}
		elsif ($seq[$k] eq "T"){$temp.= "A";}
		elsif ($seq[$k] eq "G"){$temp.= "C";}
		elsif ($seq[$k] eq "C"){$temp.= "G";}
	}
	$revcompl = reverse $temp;
	return ($revcompl);
}

sub print_double {#this subroutine does the printing for the report and the new hashes.
	my %hash = %{(shift)};
	foreach my $structure (sort keys %hash){
		foreach my $allel (sort keys %{$hash{$structure}}){
			print OUT "$structure\tFor:";
			if (exists $hash{$structure}{$allel}{amountFor}){
				print OUT $hash{$structure}{$allel}{amountFor}."\t";
			} else {
				print OUT "0\t";
			}
			print OUT "Rev:";
			if (exists $hash{$structure}{$allel}{amountRev}){
				print OUT $hash{$structure}{$allel}{amountRev}."\t";
			} else {
				print OUT "0\t";
			}
			print OUT "$hash{$structure}{$allel}{amount}\t$allel\n";#prints per marker new allels or good allels with the count
		}
	}
}

sub print_single {#this subroutine does the printing for the error hash
	my(%hash) = %{(shift)};
	my $label = shift;
	foreach my $structure (sort keys %hash){
		print OUT $label.": $structure\tFor:";
		if (exists $hash{$structure}{amountFor}){
			print OUT $hash{$structure}{amountFor}."\t";
		} else {
			print OUT "0\t";
		}
		print OUT "Rev:";
		if (exists $hash{$structure}{amountRev}){
			print OUT $hash{$structure}{amountRev}."\t";
		} else {
			print OUT "0\t";
		}
		print OUT "$hash{$structure}{amount}\n";#prints per marker new allels with the count
	}
}

sub makeMatrix {
  my ($xSize, $ySize) = @_;

  for (my $y = 0; $y < $ySize; $y++) {
    $matrix[0][$y] = $y;   # Initialise the first row with increasing numbers.
  }#for
  for (my $x = 1; $x < $xSize; $x++) {
    for (my $y = 0; $y < $ySize; $y++) {
      $matrix[$x][$y] = 0; # Initialise the rest of the elements with zero.
    }#for
  }#for
  return $matrix
}#makeMatrix

# Print an xSize * ySize matrix and the sequences.
#
# Arguments:
#   matrix ; An xSize * ySize matrix.
#   xSize  ; Number of rows in the matrix.
#   ySize  ; Number of columns in the matrix.
#   seq1   ; The first sequence in the pairwise alignment.
#   seq2   ; The second sequence in the pairwise alignment.
#
sub printMatrix {
  my ($matrix, $xSize, $ySize, $seq1, $seq2) = @_;

  print "  - ";
  for (my $y = 0; $y < $ySize; $y++) {     # Print the second sequence.
    print substr($seq2, $y, 1), " ";
  }#for
  print "\n- ";
  for (my $x = 0; $x < $xSize; $x++) {
    if ($x) {
      print substr($seq1, $x - 1, 1), " "; # Print the first sequence.
    }#if
    for (my $y = 0; $y < $ySize; $y++) {
      print $matrix[$x][$y], " ";          # Print the contents of the matrix.
    }#for
    print "\n";
  }#for
}#makeMatrix

# Global pairwise alignment.
#
# Arguments:
#   matrix ; An xSize * ySize matrix.
#   xSize  ; Number of rows in the matrix.
#   ySize  ; Number of columns in the matrix.
#   seq1   ; The first sequence in the pairwise alignment.
#   seq2   ; The second sequence in the pairwise alignment.
#
# Returns:
#   matrix ; A matrix with in position $xSize - 1, $ySize - 1 the value of the
#            global alignment.
#
sub align {
  my ($matrix, $xSize, $ySize, $seq1, $seq2) = @_;

  for (my $x = 1; $x < $xSize; $x++) {
    for (my $y = 1; $y < $ySize; $y++) {
      # Choosing matrix[x - 1][y] or matrix[x][y - 1] introduces a gap, with
      #   penalty 1. 
      # Choosing matrix[x - 1][y - 1] either introduces a mismatch, or a match,
      #   depening on the match of seq1[x - 1] and seq2[y - 1].
      $matrix[$x][$y] = min($matrix[$x - 1][$y] + 1, $matrix[$x][$y - 1] + 1,
                            $matrix[$x - 1][$y - 1] + 
                            int(substr($seq1, $x - 1, 1) ne
                                substr($seq2, $y - 1, 1)));
    }#for
  }#for
  return $matrix
}#align

# Find the minimum distance, ignoring a trailing gap in the sequence associated
#  with the number of rows in an alignment matrix. If the miniumum distance is
#  found, also return the row number.
#
# It is assumed that the number of rows is larger than the number of columns.
#
# Arguments:
#   matrix ; An xSize * ySize matrix.
#   xSize  ; Number of rows in the matrix.
#   ySize  ; Number of columns in the matrix.
#
# Returns:
#   minimum ; The minimum distance.
#   xMin    ; Row number.
#
sub findMin {
  my ($matrix, $xSize, $ySize) = @_;
  my $minimum = $ySize - 1;
  my $xMin = 0;

  # Look in the last column.
  for ($x = 0; $x < $xSize; $x++) {
    if ($matrix[$x][$ySize - 1] < $minimum) {
      $minimum = $matrix[$x][$ySize - 1];     # We found a new minimum.
      $xMin = $x;
    }#if
  }#for
  return $minimum, $xMin;
}#findMin

# Main entry point.
#

# It is important that seq2 is smaller than seq1.
#$seq1 = "AGCGGAGCGCGGGGCCCCCAGAGAGAGAGATAGATAGATCAGATACAGATAGATAGAC";
#$seq2 = "CCCXAGA";
sub al {
	my ($seq1,$seq2)=@_;
	my $xSize = length($seq1) + 1;
	my $ySize = length($seq2) + 1;
	my $matrix = makeMatrix($xSize, $ySize);

	# Notice that we mis-use the global alignment algorithm, because the matrix is
	#   initialised for local alignment.
	align($matrix, $xSize, $ySize, $seq1, $seq2);
	#printMatrix($matrix, $xSize, $ySize, $seq1, $seq2);

	# Now we can find the minimum number of mismatches and the position in seq1.
	my ($minimum, $xMin) = findMin($matrix, $xSize, $ySize);
	return $minimum, $xMin;
}
