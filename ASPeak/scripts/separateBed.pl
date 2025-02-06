#!/usr/bin/env perl
#                                       separateBed.pl
#########################################################################################
# 
# It is used for bed separatation.
#
#########################################################################################
# AUTHORS:
#
# Can Cenik, PhD 
# Alper Kucukural, PhD 
# Hakan Ozadam, PhD
#########################################################################################

use Pod::Usage; 

use strict;
use warnings;

# ./script inputFile outputDirectory

if(scalar(@ARGV) and $ARGV[0] eq "--help"){
  	pod2usage( {
		'-verbose' => 2, 
		'-exitval' => 1,
	} );
} 

if (scalar(@ARGV) < 2) 
{ 
	# when no command line options are present
	# print SYNOPSIS
	pod2usage( {
		'-verbose' => 0, 
		'-exitval' => 1,
	} );
}

# die("Run this script as\n./script inputFile outputFolder\n") if(@ARGV < 2);

my $inputFile = $ARGV[0];
my $outputFolder = $ARGV[1];

open(my $IH, $inputFile) or die("Can not open the input file $inputFile\n");

# Holds the file handlers.
# We have a separate file for each chromosome.
my %chromosomes;

while(my $line = <$IH>){
  my @lineArray = split(" ", $line);
  if(@lineArray >= 4)
  { 
    if(exists $chromosomes{$lineArray[0]}){
      my $IFILE =  $chromosomes{$lineArray[0]};
      print $IFILE $line; 
    }
    else{
      my $fileName = $outputFolder."/".$lineArray[0].".bed";
      local *FILE;
      open(FILE, ">".$fileName) or die("Can not open the file".$lineArray[0]."bed\n");
      print FILE $line;
      $chromosomes{$lineArray[0]} = *FILE;
    }
   }
}

foreach my $key (keys %chromosomes){
  close($chromosomes{$key});
}

close($IH);


__END__

=head1 NAME

separateBedGraph.pl

=head1 SYNOPSIS

This program requires two arguments.
The first one is the input bedGraph file.
The second one is the ouptput directory.
Run this script as

./separateBedGraph.pl input_bg_file output_dir

For detailed information type

./separateBedGraph.pl --help


=head1 DESCRIPTION

If the reads for all chromosomes are in one file, this script partitions this big file
into smaller files according to the chromosomes.
So in the end, each file contains the reads of a single chromosome.
The output files are named according  to the chromosomes.
For example chr1.bg holds the read counts of the first chromosome.

=head1 EXAMPLE

ex: ./separateBedGraph all_reads.bg ./separated_bg_dir/

=head1 AUTHORS

 Can Cenik, PhD 

 Alper Kucukural, PhD
 
 Hakan Ozadam

 
=head1 FILE FORMAT

A bedGraph file (both input and output) has four entries in each row separated by whitespace.

Column 1: chromosome
Column 2: start position
Column 3: end position
Column 4: read counts

Ex:

chr12 5678 5679 12
chr12 5679 5680 16
 
=head1 LICENSE AND COPYING

 This program is free software; you can redistribute it and / or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, a copy is available at
 http://www.gnu.org/licenses/licenses.html
