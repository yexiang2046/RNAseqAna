#!/usr/bin/env perl

#########################################################################################
#                                       makeRPKM.pl
#########################################################################################
# 
#  This program mergesc count files to create a RPKM file. If there are more than one
#  RNASeq file, this program get average RPKM values to create a single RPKM file.
#
#########################################################################################
# AUTHORS:
#
# Can Cenik, PhD 
# Alper Kucukural, PhD 
# Hakan Ozadam, PhD
#########################################################################################

BEGIN
{
  use Cwd 'abs_path';
  use File::Basename;
  my($mainScriptName, $scriptDirectory, $suffix) = fileparse(abs_path($0)); # Get the directory of the working script
  my @scriptPathPieces = split(/\//,$scriptDirectory); # Make an array of subdirectories
  pop(@scriptPathPieces); # Since this script is under "/dir1/dir2/.../dirN/scripts" directory, it pops "scripts" from the array
  my $libPath = join("/", @scriptPathPieces); #This is the main directory containing ASPeak files. Hence  it is our library patth 
  push (@INC, $libPath); #Add the library path to default inclusion paths.
}

############## LIBRARIES AND PRAGMAS ################

 use List::Util qw[min max];
 use strict;
 use File::Basename;
 require lib::io;
 require lib::func;
 use Getopt::Long;
 use Pod::Usage; 
 use Class::Struct;
 use POSIX qw/floor/;
 
#################### VARIABLES ######################
 my $input            = "";
 my $output           = "";
 my $help             = "";
 my $print_version    = "";
 my $version          = "1.0.0";
################### PARAMETER PARSING ####################

my $cmd=$0." ".join(" ",@ARGV); ####command line copy

GetOptions( 
	'input=s'        => \$input,
	'output=s'       => \$output,
	'help'           => \$help, 
	'version'        => \$print_version,
) or die("Unrecognized optioins.\nFor help, run this script with -help option.\n");

if($help){
    pod2usage( {
		'-verbose' => 2, 
		'-exitval' => 1,
	} );
}

if($print_version){
  print "Version ".$version."\n";
  exit;
}

pod2usage( {'-verbose' => 0, '-exitval' => 1,} ) if ( ($input eq "") or ($output eq ""));	


################### MAIN PROGRAM ####################
#    Read input files to merge RPKM values for each interval and write it into output file

my @files=split(/:/, $input);

my @RPKMVals=();
for(my $i=1; $i<@files; $i++)
{
    my %tmphash=();

    getRPKMfile(\%tmphash, $files[$i]);
    push(@RPKMVals, \%tmphash);
}

if (!(-s $files[0])) {
    	print "Please run StepPrepInput for $files[0]!!!\n";
	exit;
}
open(OUT, ">$output");
open (IN, $files[0]);
while (my $line=<IN>)
{  
    chomp($line);
    my @arr=split(/[\t\s]+/, $line);
    my $tot=$arr[5];
    #print "$tot+";
    for(my $i=1; $i<@files; $i++)
    {
	my %tmphash=%{$RPKMVals[$i-1]};
	#print $tmphash{$arr[1]};
	$tot+=$tmphash{$arr[1]};
    }
    #print "=$tot\n";
   
    print OUT $arr[0]."\t".$arr[1]."\t".($tot/@files)."\n";
}
close(IN);
close(OUT);

sub getRPKMfile
{
    my ($hash, $file)=@_;

    if (!(-s $file)) {
    	print "Please run StepPrepInput for $file!!!\n";
	exit;
    }
    open(IN, $file);
    while (my $line=<IN>) 
    {
      chomp($line);
      my @arr=split(/[\t\s]+/, $line);
      #print "$arr[1] = $arr[5]\n";
      ${$hash}{$arr[1]}=$arr[5];
    }
    close(IN);
}
__END__


=head1 NAME

makeRPKM.pl

=head1 SYNOPSIS

makeRPKM.pl -i input file <count format> -o output file <RPKM format>

makeRPKM.pl -help

makeRPKM.pl -version

For help, run this script with -help option.

=head1 OPTIONS

=head2 -i  input file <BED format> 

It merges input count file for each chromosme separately.

Count File
    Column 1  Chromosome\\
    Column 2  Interval name\\
    Column 3  Interval length\\
    Column 4  Read Count\\
    Column 5  Sum of the square of the reads\\
    Column 6  RPKM Value\\


=head2 -o output file <BG format>

RPKM FILE format
  Column 1: The name of the chromosome (eg. chr1, chrX) \\
  Column 2: Interval Name \\
  Column 3: RPKM\\


=head2 -help

Display this documentation.

=head2 -version

Display the version

=head1 DESCRIPTION

  This program mergesc count files to create a RPKM file. If there are more than one
  RNASeq file, this program get average RPKM values to create a single RPKM file.


=head1 EXAMPLE

./makeRPKM.pl -i your_input_dir/chr1.intervalcount:your_input_dir2/chr1.intervalcount -o your_output_dir/chr1.csv 

=head1 AUTHORS

 Can Cenik, PhD 

 Alper Kucukural, PhD
 
 Hakan Ozadam
 
 
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
