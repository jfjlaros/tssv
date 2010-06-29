#!/usr/local/bin/perl -w

#use strict;
use List::Util qw[min max];

#place the FLDO library in the same directory as the perl script, give the fasta reads file as argument 1, outputs rappor.txt
if(scalar(@ARGV) != 4){
	print "Please use the correct parameters\nUsage: FLDO_tool.pl -FLDO reads file- -marker library- -allowed mismatches for 25nt (recommended=2)- -run code-\n";
	exit();
}

open (REF, $ARGV[1]) || die "cannot open FLDO library\n";#opens the library with markers and repeat patterns
my %hash=();
my %regular=();
while (<REF>){
	$_ =~ s/\n|\r//g;
	my @element=split(/\t/,$_);
	$hash{$element[0]}{ref1}=$element[1];#places start marker in hash under ref1
	$hash{$element[0]}{ref2}=$element[2];#places end marker in hash under ref2
	for (my $i=3;$i<@element;$i++){
		@repeat_structure=split(/ /,$element[$i]);
		my $reg="";
		for (my $j=0;$j<@repeat_structure;$j+=3){
			$reg.="($repeat_structure[$j]){".$repeat_structure[($j+1)].",".$repeat_structure[($j+2)]."}";#builds regular expression for use in pattern detection
		}
		$regular{$element[0]}{$i}=$reg;#stores the regular expression in "regular"
	}
}
close REF;
my $read=0;
my $sequence="";
my @al1="";
my @al2="";
my %rapport=();
my %error=();
my %new=();
my $correct_allel=0;
my $different_structures=0;
my $no_end=0;
my $no_start=0;
my $new_allel=0;
open (IN, $ARGV[0]) || die "cannot open reads file\n";#opens the reads file
open (OUT, ">rapport_$ARGV[3].txt") || die "cannot open file for writing rapport\n";#opens the output file for writing rapport
open (OUT2, ">rapport_$ARGV[3]_whole_sequences.fa") || die "cannot open file for writing rapport\n";#opens the output file for whole sequences
while (<IN>){
	$_ =~ s/\n|\r//g;
	if ($_ =~ m/^>/){
		$read++;
		if ($sequence ne ""){
			my %score=();
			my $found=0;
			foreach my $structure (keys %hash){
				if ($structure eq ""){print "|$structure|\n";next;}
				if (length ($sequence) > length ($hash{$structure}{ref1})){#it is important that the length of the imput sequence is longer than the aligned structure
					@al1=al($sequence,$hash{$structure}{ref1});#gives the edit distance and position on sequence of ref1
					$score{$structure}{begin}=$al1[0];
					$score{$structure}{startpos}=$al1[1];
					$rev_sequence=reverse($sequence);
					$rev_ref=reverse($hash{$structure}{ref2});
					@al2=al($rev_sequence,$rev_ref);#gives the edit distance and position from end on sequence of ref2
					$score{$structure}{end}=$al2[0];
					$score{$structure}{endpos}=$al2[1];
					$found=1;
				}
			}
			if ($found==0){next;}
			my $smallest_start="";
			my $smallest_end="";
			my $ref_start="";
			my $ref_end="";
			foreach my $reference (keys %score){#this loop determines the smallest edit distance per read for each marker and saves the marker in ref_start and ref_end
				if ($reference eq ""){print ":$reference:\n";next;}
				if ($smallest_start ne ""){
					if ($score{$reference}{begin}<$smallest_start){
						$ref_start=$reference;
						$smallest_start=$score{$reference}{begin};
					}
				} else {
					$smallest_start=$score{$reference}{begin};
					$ref_start=$reference;
				}
				if ($smallest_end ne ""){
					if ($score{$reference}{end}<$smallest_end){
						$ref_end=$reference;
						$smallest_end=$score{$reference}{end};
					}
				} else {
					$smallest_end=$score{$reference}{end};
					$ref_end=$reference;
				}
			}
			if ($smallest_start <= (length ($hash{$ref_start}{ref1})/25*$ARGV[2]) && $smallest_end <= $ARGV[2]){#if the edit distance is small enough for start and end
				print OUT2 ">$refstart\n$sequence\n";
				if ($ref_start eq $ref_end){#and if the same start and end marker is found
					$al_repeat=substr($sequence,$score{$ref_start}{startpos},(length($sequence)-$score{$ref_start}{endpos}-$score{$ref_start}{startpos}));#substract the repeat
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
							$rapport{$ref_start}{$allel}{amount}++;#store existing allels
						}
					}
					if ($reg_found==0){
						$new{$ref_start}{$al_repeat}{amount}++;#store new allels
						$new_allel++;
					}
				} else {
					$different_structures++;#count when different markers are best in start and end
				}
			} elsif ($smallest_start <= (length ($hash{$ref_start}{ref1})/25*$ARGV[2])){#if no end is found
				$no_end++;#count total reads with no end marker
				$error{$ref_start}{noend}++;#store per marker
			} elsif ($smallest_end <= $ARGV[2]){
				$no_start++;#count total reads with no start marker
				$error{$ref_start}{nostart}++;
			} else {
				$nothing++;#count where no marker is found 
			}
		$sequence="";
		}
	} else {
		$sequence.=$_;#build sequence
	}
}
close IN;
print OUT "aantal reads: $read\n";
print OUT "beginning but no end: $no_end\n";
print OUT "end but no beginning: $no_start\n";
print OUT "no beginning no end: $nothing\n";
print OUT "different structure: $different_structures\n";
print OUT "new allel: $new_allel\n";
print OUT "good allels: $correct_allel\n";

print OUT "\n################rapport\n";
foreach my $structure (sort keys %rapport){
	foreach my $amount (sort keys %{$rapport{$structure}}){
		print OUT "$structure\t$rapport{$structure}{$amount}{amount}\t$amount\n";#prints per marker the pattern matching allels with the count
	}
}
print OUT "\n################error\n";
foreach my $structure (sort keys %error){
	if ($error{$structure}{nostart}){	
		print OUT "no start\t$structure\t$error{$structure}{nostart}\n";#prints per marker the count where no start is found
	}
	if ($error{$structure}{noend}){
		print OUT "no end\t$structure\t$error{$structure}{noend}\n";#prints per marker the count where no end is found
	}
}
print OUT "\n################new\n";
foreach my $structure (sort keys %new){
	foreach my $amount (sort keys %{$new{$structure}}){
		print OUT "$structure\t$new{$structure}{$amount}{amount}\t$amount\n";#prints per marker new allels with the count
	}
}
close OUT;
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
