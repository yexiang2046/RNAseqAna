#!/usr/bin/env perl
#########################################################################################
#                                       calcFDR.pl
#########################################################################################
# 
# False discovery rate calculation for each region.
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
 use Class::Struct;
 use POSIX; 
 use lib "./";
 require lib::io;
 require lib::func;
 use Getopt::Long;
 use Pod::Usage; 
 use File::Basename;

struct( exon => [ identifier => '$', len => '$', count => '$', sumSquare => '$' , RPKM => '$', average => '$' ]);

#################### VARIABLES ######################

my $libfile       = "";
my $controlfile   = ""; 
my $version = "2.0.1";
my $help = "";
my $print_version = "";

################### CONSTANTS  #####################

my $DEBUG       = 0;
my $PEAKCOLNUM  = 12;

################### PARAMETER PARSING ####################

my $cmd = $0." ".join(" ",@ARGV); ####command line copy


GetOptions(
	'libfile=s'           => \$libfile,
	'controlfile=s'       => \$controlfile,
	'help'                => \$help, # request help
	'version'             => \$print_version,
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

pod2usage( {'-verbose' => 0, '-exitval' => 1,} ) if ( ($libfile eq "") or ($controlfile eq ""));	


################### M A I N      P R O G R A M #############################################
#	Get the counts and the lengths of each interval.
#############################################################################################

my %controlhash=();
calcCounts(\%controlhash, "$controlfile");
my %libhash=();
calcCounts(\%libhash, "$libfile");
my %fdrlib=();
calcFDR(\%fdrlib, \%libhash, \%controlhash);
writePeaks(\%fdrlib, "$libfile");
#my %fdrcontrol=();
#calcFDR(\%fdrcontrol, \%controlhash, \%libhash);
#writePeaks(\%fdrcontrol, "$controlfile");

sub writePeaks
{
    my ($fdr, $file) = @_;
    open(IN, $file);
    open(OUT,">$file.fdr");
    print OUT "#Chr\tStrand\tInterval\tInterval_Length\tTag_Count\tMaximum_Height\tSpan_Size\tStart\tEnd\tRPKM\tWeighted_Center\tp-Value\tFDR\n";

    while (my $line=<IN>)
    {
	chomp($line);
        if ($line!~/^#/)
	{
	  my @arr=split(/\t/, $line);
	  my $pval = sprintf("%e", $arr[$PEAKCOLNUM-1]);
	  $pval=~/\-(\d+)$/;
	  my $key = floor($1);
	  my $f=sprintf("%0.2e", ${$fdr}{$key});
	  $f = "<".$f if (${$fdr}{$key}>0);
	  my $txt="";
	  for (my $i=0; $i<$PEAKCOLNUM; $i++)
	  {
	    $txt.=$arr[$i]."\t";
	  }
	  print OUT $txt.$f."\n";
	}
    }
    close(IN);
    close(OUT);
    $file=~/(.*\/)/;
    my $dir=$1;
    
    if (-s "$dir/control*") {
	 `rm $1/control*`;
    }
    
    `mv $file.fdr $file`; 
}
sub calcFDR
{
    my ($fdr, $libhash, $controlhash) = @_;
    my ($v, $s)=(0, 0);
    
    foreach my $key ( sort { $b <=> $a } keys %{$libhash})
    {
	
	if (exists ${$libhash}{$key})
	{
	  $s += ${$libhash}{$key};
	}
	if (exists ${$controlhash}{$key})
	{
	  $v += ${$controlhash}{$key};
	}
	if (($v+$s)>0)
	{
	  ${$fdr}{$key}=$v/($v+$s);
	   print $key.":$v:$s:".${$fdr}{$key}."\n";
	}
    }
    my $minfdr=1;
    foreach my $key (sort {$a <=> $b} keys %{$fdr})
    {
       if (${$fdr}{$key}<$minfdr)
       {
	  $minfdr=${$fdr}{$key};
       }
       ${$fdr}{$key}=$minfdr;
       print $key.":".${$fdr}{$key}."\n";
    }
}

sub calcCounts
{
    my ($hash, $file) = @_;
    open(IN, $file);
    while (my $line=<IN>)
    {
	chomp($line);
        if ($line!~/^#/)
	{
	 my @arr=split(/\t/, $line);
	 my $pval = sprintf("%e", $arr[$PEAKCOLNUM-1]);
	 $pval=~/\-(\d+)$/;
	 my $key = floor($1);
	 ${$hash}{$key}++;
	}
    }
    close(IN);
}

__END__

=head1 NAME

calcFDR.pl.pl

=head1 SYNOPSIS 

 calcFDR.pl.pl -l <library peaks>  -c <control peaks>

 calcFDR.pl.pl -help

 calcFDR.pl.pl -version

 For help, run this script with the -help option.


=head1 OPTIONS


=head2 -l <library peaks>

Library peaks 

=head2 -c <control peaks>

Control peaks

=head2 -help

Displays this documentation.

=head2 -version

Prints the version  number.


=head1 DESCRIPTION

Using the counts of the intervals (coming from the RIPSeq data) and the abundance of the intervals
(in the RPKM file coming fom the RNASeq data), this script estimates the parameters of the Negative binomial distribution for each interval.
It uses the method of moments for estiamting the parameters r and p.

=head1 INPUT

 This program has two inputs.

=head2 count files directory
  
 It contains the counts of intervals obtained from Rip / Clip Seq data
 The contents of the count files are as follows:

=over 5

=item Column 1 : Chromosome 

=item Column 2 : Interval name

=item Column 3 : Interval length

=item Column 4 : Read Count (the sum of counts for each nucleotidde position)

=item Column 5 : Sum of the square of the reads coming from each nucleotide position

=item Column 6 : RPKM value

=back

=head2 RPKMfile

 RPKMfile is obtained from RNAseq data. It provides the abundance information of the intervals.
 The contents of the RPKM file are as follows.
 
=over 5

=item Column 1 : The name of the chromosome (eg. chr1, chrX)

=item Column 2 : Interval Name 

=item Column 3 : RPKM value


=back

=head1 OUTPUT

The output is written to a single file containing the estimated parameters of the Negative Binomial distribution for each interval.
The columns of the output file are as follows.

=over 5

=item Column 1: interval name (for example: NR_024540_utr3_8_0_chr1_8131_r)

=item Column 2: p (the first NB parameter, 0 < p < 1)

=item Column 3: r (the second NB parameter, r > 0)

=item Column 4: RPKM of the interval

=item Column 5: Interval Count

=item Column 6: Interval Length

=back


=head1 EXAMPLE

Suppose that  the interval count files are in the directory /home/user/intervalcounts
the output directory is /home/user/NBParameters
and the RPKM file is /home/user/RPKM/RPKM.txt

Then the parameters of the NB can be estimated by

estimateNBParameters.pl -i /home/user/intervalcounts -o /home/user/NBParameters -rpkm /home/user/RPKM/RPKM.txt

=head1 AUTHORS

 Can Cenik, PhD 

 Alper Kucukural, PhD
 
 Hakan Ozadam, PhD
 
 
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
