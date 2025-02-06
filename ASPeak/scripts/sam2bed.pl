#!/usr/bin/env perl
#########################################################################################
#                                      sam2bed.pl
#########################################################################################
#
# Conversion from sam to bed format
#
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
 use POSIX;
 use Class::Struct;
 use File::Basename;
 use Getopt::Long;
 use Pod::Usage;
 use strict;
 require lib::io;
 require lib::func;

#################### CONSTANTS ###################### 
my $VERSION = '2.0.2';

#debuging  purposes only
my $INPUT_DEBUG = 0;

#################### VARIABLES ######################
my $input="";
my $output ="";
my $help="";
my $print_version="";

############## LIBRARIES AND PRAGMAS ################RAMETER PARSING ####################
my $cmd = $0." ".join(" ",@ARGV); ####command line copy 

GetOptions(
        'inout=s'   => \$input, #
        'output=s'  => \$output,
        'help'      => \$help,
        'version'   => \$print_version,
) or die " unrecognized option(s)!! please refer to the help documentation\n\n";

if ($print_version)
{
        print "Pipeline main script, version $VERSION\n\n";
        exit;
}

if ($help) {
        # print entire POD
        pod2usage( {
                '-verbose' => 2,
                '-exitval' => 1,
        } );
}

pod2usage( {'-verbose' => 0, '-exitval' => 1,} ) if ( ($input eq "") or ($output eq ""));

unless ($input)
{
        $input = shift @ARGV or die "No input file specified! \n use --help\n";
}

unless ($output)
{
        $output = shift @ARGV or die "No output file specified! \n use --help\n";
}


  open(my $infile, $input);
  open(my $outfile, ">$output");
  my @cache;

  while(<$infile>){
	next if(/^@/);#Skip header
	chomp;
	my ($name, $flag, $chrom, $pos, $mapq, undef, undef, undef, undef, $read) = split("\t");
	#If chrom is not in the 3rd column or the definition is different please fix it here.
	
	next if $flag & 4;#Unmapped read. 
	
	#Query strand is in the 5th(index == 4) bit... 
	#i.e. 2*4 = 16 bitwise 16 & 16 == 16 else 0
	my $strand = ($flag & 16) ? '-' : '+';

	#SEQ_REGION_NAME, START, END, FEATURE_NAME, SCORE STRAND
	push @cache, join("\t", ($chrom, $pos, ($pos +length($read) -1), $name, $mapq, $strand));

	if(scalar(@cache) == 1000){
	  print $outfile join("\n", @cache)."\n";
	  @cache = ();
	}
  }

  print $outfile join("\n", @cache)."\n";

  close $infile;
  close $outfile;

__END__


=head1 NAME

sam2bed.pl

=head1 SYNOPSIS

sam2bed.pl -i input.sam -o output.sam

sam2bed.pl -help

sam2bed.pl -version  

For help, run this script with -help option.

=head1 OPTIONS

=head2 -i <inputfile>

Input file in sam format.

=head2 -o <library name>
 
Output file in bed format.

=head2 -help

Display this documentation.

=head2 -version

Display the version

=head1 DESCRIPTION

It converts the sam files to bed format.


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
