#!/usr/bin/env perl
#########################################################################################
#                                       estimateNBParameters.pl
#########################################################################################
# 
# Given interval counts (coming from RNASeq data) and the sorted RPKM values 
# (coming from RipSeq data), this program calculates the RPKM value of each interval.
# For each interval, we compute the average count by dividing the count 
# by the length in the count file. 
# 
#
# For detailed information and license, run this script with -help option
# or see the plain old documentation at the bottom.
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

my $RPKMfile = "";
my $countdir = ""; 
my $countFile = "";
my $parameterFile = "";
my $version = "2.0.1";
my $help = "";
my $print_version = "";

# Size of the frame to average the counts to find lambda
my $listcount = 0;
my $radius = 0; # the length of the radius around each interval in estimating the parameters 

# The cutoff percentage in calculating  the lambda value of an exon.
# For instance if the listcount is 1000 and $cutoffPercentage is 5,
# then for each exon, we consider the exons in the RPKM file 500 above and below
# the exon. Then in calculating the lambda value, we disregard 50 (5 percent of 1000) exons having
# highest average and 50 exons having the lowest average (5 percent of 1000)
my $cutoffPercentage = 5;

my $DEBUG = 0;

################### PARAMETER PARSING ####################

my $cmd = $0." ".join(" ",@ARGV); ####command line copy


GetOptions(
	'rpkm=s'    => \$RPKMfile,
	'i=s'       => \$countdir,
	'radius=i'  => \$radius, 	
	'help'      => \$help, # request help
	'o=s'       => \$parameterFile,
	'version'   => \$print_version,
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

pod2usage( {'-verbose' => 0, '-exitval' => 1,} ) if ( ($RPKMfile eq "") or ($countdir eq "") or ($parameterFile eq ""));	


my($filename, $outdir, $suffix) = fileparse($parameterFile);

my $logDir = $outdir."/LOG"; 
`mkdir -p $outdir; mkdir -p $logDir`;

my $logname = "summary_estimateNBParameters.log";
my $num = 1;   

if (-s "$logDir/$logname.1")
{
   $num = `ls $logDir/$logname.*|awk '{split(\$1, a, "."); print a[3]; }'|sort -n|tail -n 1`;
   $num++;
}

## The default value of the number of intervals in determining the parameters of an interval is 1000
$listcount = 1000 if ($listcount == 0); 
# The default value of radius is 500
$radius = 500 if ($radius == 0);

open(my $LOG, ">$logDir/$logname.$num") or die("Error: Couldn't open the log file $logDir/$logname.$num\n");

lib::io::wlog($LOG, "program started at: ");    
lib::io::wlog($LOG, "The command used is:");
lib::io::wlog($LOG, "$cmd\n");
lib::io::wlog($LOG, "Parameters specified:");
lib::io::wlog($LOG, "count dir: $countdir ");
lib::io::wlog($LOG, "outdir: $outdir");
lib::io::wlog($LOG, "RPKMfile: $RPKMfile");
lib::io::wlog($LOG, "estimation radius: $radius");

################### M A I N      P R O G R A M #############################################
#	Get the counts and the lengths of each interval.
#	These counts are obtained from the files in  the $countdir.
#	Get the RNAseq data fom the RPKM file and sort the intervals
#	by their RPKM values.
#	For each interval, we compute the average count by dividing the count by the length.
#	Then, in order to find the parameters p and r of an interval,   
#	we find the average value of the intervals having "similar" RPKM values.
#	We discard the outliers by the top 5% and bottom 5% of the intervals in the sorted average list. 
#       The parameters of the Negative Binomial distribution are estimated using the method of moments.
#	As an example, suppose that list count is 1000 and we want  to find the
#	lambda of the interval which is in the 2000 th position in the sorted RPKM list. 
#	We sort the intervals by their averages in %$hashIntervalCounts and %hashExonLengths 
#	from the position 1500 to 2500 in the sorted RPKM list.
#	Then we start  from 50 (exclude the smallest 5%) take the counts and the sum of squares
#	till the position 950 (exclude the largest 5%). 
#	We estimate the moments m1 = totalcount / total length and m2 = sum_of_squares / total_length
#	Then we estimate the parameters p and r by
#	p = m1 / (m2 - m1^2)
# 	r = p*m1 / (1-p)
#############################################################################################

open (OUT, ">$parameterFile") or die("Error: Couldn't open the file $parameterFile for writing.\n");

my %hashIntervalCounts = (); # sum of tag counts
my %hashIntervalCountSquares = (); # sum of the square of tag counts
my %hashIntervalLengths = ();
my $c = 0; #size of the above hashes

# We already calculated the number of read counts 
# and lengths of each interval (using findCounts.pl) and saved the results in
# .intervalcounts file. We read them into the hashes %$hashIntervalCounts and %$hashIntervalLengths
lib::io::wlog($LOG, "getCounts");
getIntervalCounts($countdir,"intervalcount", \$c, \%hashIntervalCounts, \%hashIntervalCountSquares , \%hashIntervalLengths);
lib::io::wlog($LOG, "getCounts: $c DONE");

my %hashRPKM = (); # Chrom_gene_exon => exonRPKM (or geneRPKM, depending on $geneexon)

# Its elements are of the form: Chrom_gene_exon
# These exons are sorted according to their RPKMS
my @sortedRPKM = (); # 

#Number of elements in %hashRPKM and @sortedRPKM
my $hashSize = 0;

#DEBUGGING PURPOSES
# foreach my $key (keys %hashIntervalCounts){
#   print "$key -> $hashIntervalCounts{$key} \t $hashIntervalCountSquares{$key} \n";
# }
# #exit

# Get the RPKMS in sorted  form.
lib::io::wlog($LOG, "getRPKM");

getRPKM($RPKMfile, \$hashSize, \%hashRPKM, \@sortedRPKM);
lib::io::wlog($LOG, "getRPKM: $hashSize DONE");

# #Debugging purposses
# foreach my $an (@sortedRPKM){
#   print "$an -> $hashRPKM{$an}\n";
# }
# exit;

my @arr_Intervals = (); # The elements of this array are the structs defiend above

my $rightSize = 0;
my $leftSize = 0;

# Our frame is centered at the current exon with a radius of $listcount / 2 
my $end = $radius + 1;

 $end = scalar(@sortedRPKM) if ($end > scalar(@sortedRPKM) );

# Initialize @arr_Intervals
initIntervalArray(\@arr_Intervals, \@sortedRPKM, \%hashIntervalCounts, \%hashIntervalCountSquares, \%hashIntervalLengths, 0, $end);
$rightSize = $end - 1;

# printExonAray(\@arr_Intervals);

  my $tempExon;
  my $tempIdentifier;
  my $tempLength;
  my $tempCount;
  my $tempAverage;

  #First find the lambda for the first element in the sortedRPKM list
  $tempIdentifier = $sortedRPKM[0];

  print "listcount = $listcount , radius = $radius\n\n" if $DEBUG;
  
  my ($p, $r) = calculateParameters(\@arr_Intervals, $rightSize + $leftSize + 1);
  
  print "For, the first interval, $tempIdentifier,  p = $p, r = $r\n" if $DEBUG;
  print "===============================================\n\n" if $DEBUG;

  print OUT $tempIdentifier."\t".$p."\t".$r."\t".$hashRPKM{$tempIdentifier}."\t".$hashIntervalCounts{$tempIdentifier}."\t".$hashIntervalLengths{$tempIdentifier}."\n"
      if (exists $hashIntervalCounts{$tempIdentifier});

for(my $i = 1; $i < @sortedRPKM; $i++)
{
  my $totalCount = 0;
  my $totalLength = 0;
     # if the right Size is less than radius, then we have no room
    # to enlarge our  frame on the right.
    # Since the position is shifted 1 unit to the right, right size decreases and left size increases.  
    if(($rightSize < $radius) and ($leftSize < $radius) )
    {
      $rightSize--; $leftSize++;
    }
    elsif( ($rightSize == $radius) and ($leftSize < $radius) )
    {
      if( ($i + $rightSize) < scalar(@sortedRPKM) ) 
      {
	$tempIdentifier = $sortedRPKM[$i + $rightSize];
	pushIntervalToArray(\@arr_Intervals, $tempIdentifier, \%hashIntervalCounts, \%hashIntervalCountSquares,\%hashIntervalLengths);
	$leftSize++;
      }
      else
      {
	$rightSize--; $leftSize++; 
      }	
    }
    elsif( ($rightSize < $radius) and ($leftSize == $radius) )
    {
      shift(@arr_Intervals);
      $rightSize--;
    }
    elsif( ($rightSize == $radius) and ($leftSize == $radius) )
    {
      if($i + $rightSize < scalar(@sortedRPKM))
      {
	shift(@arr_Intervals);
	$tempIdentifier = $sortedRPKM[$i + $rightSize];
	pushIntervalToArray (\@arr_Intervals, $tempIdentifier, \%hashIntervalCounts, \%hashIntervalCountSquares, \%hashIntervalLengths);
      }
      else
      {
	shift(@arr_Intervals);
	$rightSize--;
      }
     }

   
  my $n = $sortedRPKM[$i];

  if (exists $hashIntervalCounts{$n})
  {
    print "\nCalculating the parameters for $n:\n" if $DEBUG;
    ($p, $r) = calculateParameters(\@arr_Intervals, $rightSize + $leftSize + 1);
    print "For $sortedRPKM[$i],  p = $p, r = $r\n" if $DEBUG;
    print "===============================================\n\n" if $DEBUG;
    print OUT $n."\t".$p."\t".$r."\t".$hashRPKM{$n}."\t".$hashIntervalCounts{$n}."\t".$hashIntervalLengths{$n}."\n";  
  }
  else
  {
    print "Warning: hashIntervalCounts does not exist for the key ".$sortedRPKM[$i].". Skipping this.\n" if ($DEBUG); 
  }
}

lib::io::wlog($LOG, "PROGRAM DONE"); 
close(OUT);
close($LOG);


#####################################################################################
#######################  F u n c t i o n s ################################
#####################################################################################

# RPKM FILE
#  Column 1: The name of the chromosome (eg. chr1, chrX) \\
#  Column 2: Interval Name \\
#  Column 3: RPKM\\


sub getRPKM
{
  my ($file, $size, $hashRPKM, $arraySortedRPKM) = @_;
  
  open(my $FH, $file) or die("getRPKM(): Error: Couldn't open the file $file\n");
  my $i = 0;
  
  while(<$FH>)
  {
    chomp($_);
    my @fields = split(/\s+/,$_);
    
    # Get only the nonzero RPKM values
    # Also ignore the rows having less than three elements.
    if( (scalar(@fields) >= 3) and ($fields[2] > 0) )
    {
      ${$hashRPKM}{$fields[1]} = $fields[2];
      $i++;
    }
    
  }
  
  close($FH);
  
  @{$arraySortedRPKM} = sort { ${$hashRPKM}{$a} <=> ${$hashRPKM}{$b}} keys %{$hashRPKM};
  
  ${$size} = $i;
}


###############################################################################
########### Function: g e t Interval C o u n t s  ##############################
###############################################################################
#It reads interval counts

#Count File
#    Column 1  Chromosome\\
#    Column 2  Interval name\\
#    Column 3  Interval length\\
#    Column 4  Read Count
#    Column 5  Sum of the square of the reads	

sub getIntervalCounts
{
  my ($countDir, $extension, $c, $hashIntervalCounts, $hashIntervalCountSquares, $hashIntervalLengths) = @_;
  
  my @countFiles = split(/\n/,`ls $countdir/*.$extension 2>/dev/null`);
  print("Warning: getIntervalCounts(): No files are found in the directory $countdir with the extension $extension\n") if (scalar(@countFiles) == 0);
  my $i = 0;
  
  print("The counts files are\n".join("\n",@countFiles)."\n") if $DEBUG;
  
  foreach my $countfile (@countFiles)
  {
    open(my $CH, $countfile) or die("Error: getIntervalCounts(): Couln't open the count file $countfile\n");  
  
    while(<$CH>)
    {
      chomp($_);
      my @fields = split(/\s+/,$_);
    
      if( scalar(@fields) >= 3 )
      {
	${$hashIntervalCounts}{$fields[1]} = $fields[3];
	${$hashIntervalCountSquares}{$fields[1]} = $fields[4];
	${$hashIntervalLengths}{$fields[1]} = $fields[2];
	$i++;
      }
    }
  
    close($CH);
  }
  
  ${$c} = $i;
}


###############################################################################
########### Function: calculate Parameters  ###################################
###############################################################################
# Using the method of moments, aproximates the parameters of the negative binomial distribution.
# The elements of arr_Intervals is the struct defined above containing interval data.
# 
# Its output is the array ($p, $r) which are the parameters of the NB dist.
# It sorts the array by the average count, namely the abundance.
# The outliers are discarded by adjusting start and end values according to the cutoffpercentage

sub calculateParameters
{
  my ($arr_Intervals, $listcount)=@_;
  
  if(scalar(@{$arr_Intervals}) ==0)
  {
    print "Warning: calculateParameters(): The input array \$arr_Intervals is empty. Returning p=0.5, r =1\n";
    return(0.5,1);
  }
  # The first is r the second is p
  my ($p, $r) = (-1, -1);
  
  my @sorted=sort compareIntervalsByAverage @{$arr_Intervals};
  my $start = floor($listcount * ($cutoffPercentage / 100) );
  my $end = floor($listcount * ( (100 - $cutoffPercentage) / 100) );
  # my $start = 0; my $end = $listcount; 
  
  my ($sum, $sum_of_squares, $total_length) = (0,0,0);

  for(my $i = $start; $i <= $end; $i++) 
  {
      if (defined $sorted[$i])
      {
	$sum += $sorted[$i]->count;
	$sum_of_squares += $sorted[$i]->sumSquare;
	$total_length += $sorted[$i]->len;
      }
  } 

   #DEBUG
   print "!!!! Begin  Calc Lambda\n" if $DEBUG;
   print "   array size = $listcount\t start(inclusive) = $start \t end(inclusive) = $end\n" if $DEBUG;
   print "   The sorted Array is:\n" if $DEBUG;
   print "   " if $DEBUG;
   printExonAray(\@sorted) if $DEBUG;
   
   my ($moment1, $moment2) = (1,1);
   
  if ($total_length > 0)
  {  
    $moment1 = $sum / $total_length;
    $moment2 = $sum_of_squares / $total_length;
  }
  else{
    return(0.5 , 1);
  }
    
  $p = $moment1 / ($moment2 - $moment1**2) if( ($moment2 - $moment1**2) > 0 ) ;
  $p = 0.5 if( ($p < 0.001) or ($p > 0.9999 )); #See the note below
  $r = ($p*$moment1) / (1 - $p); 
  return ($p, $r);
  #NOTE
  #The method of moments can produce negative $p or $p >1  in some data sets.
  # In such a case we let p = 0.5 and r be the first moment i.e., the average value
}


###############################################################################
########### Function: i n i t E x o n A r r a y   #############################
###############################################################################

sub initIntervalArray
{
  my ($arr_Intervals, $sortedRPKM, $hashIntervalCounts, $hashIntervalCountSquares, $hashIntervalLengths, $start, $end) = @_;

  my $tempIdentifier;
  
  for(my $i = $start; $i < $end; $i++)
  {
    $tempIdentifier = ${$sortedRPKM}[$i];
    pushIntervalToArray($arr_Intervals, $tempIdentifier, $hashIntervalCounts, $hashIntervalCountSquares, $hashIntervalLengths) if(exists($hashIntervalCounts{$tempIdentifier}) );
  }
}

###############################################################################
########### Function: c o m p a r e E x o n s B y A v e r a g e    ############
###############################################################################
# This is the function that the sort function uses to sort the exons.
sub compareIntervalsByAverage
{
	$a->average <=> $b->average;
}

###############################################################################
########### Function: p u s h I n t e r v a l T o A r r a y    ########################
###############################################################################
sub pushIntervalToArray
{
  my ($arr_Intervals, $identifier, $hashIntervalCounts, $hashIntervalCountSquares, $hashIntervalLengths) = @_;
  
  my $tempExon;

  my ($tempLength, $tempCount, $tempAverage, $tempSumSquare) = (0, 0, 0, 0);
  
   if( exists(${$hashIntervalLengths}{$identifier}) and exists(${$hashIntervalCounts}{$identifier}) and exists(${$hashIntervalCountSquares}{$identifier}) )
   {
      $tempLength = ${$hashIntervalLengths}{$identifier};
      $tempCount = ${$hashIntervalCounts}{$identifier};
      $tempSumSquare = ${$hashIntervalCountSquares}{$identifier};
   }
   else
   {
    print "Push: Warning: identifier $identifier does not exist in the hashes hashIntervalLengths hashIntervalCounts hashIntervalCountSquares\n" if ($DEBUG);
   }
   
   $tempAverage = ($tempCount / $tempLength) if($tempCount > 0 and $tempLength > 0);

   $tempExon = exon->new(identifier => $identifier, len => $tempLength, count => $tempCount, sumSquare =>$tempSumSquare , average => $tempAverage);
   push(@{$arr_Intervals}, $tempExon);
}

#Debug Purposes
sub printExonAray
{
  my ($exons) = @_;
  
  print "The Intervals are:\n";
  print "Exon \t\t Exon_Average \t\t Count \t\t Squares \t\t Length\n";
  
  my $i = 0;
  
  foreach my $currentExon (@{$exons})
  {
    print $i++." : ".$currentExon->identifier."\t\t".$currentExon->average."\t\t".$currentExon->count."\t\t".$currentExon->sumSquare."\t\t".$currentExon->len."\n";
  }
  print "End of Intervals-----------------------------------------------\n\n";
}


__END__

=head1 NAME

estimateNBParameters.pl

=head1 SYNOPSIS 

 estimateNBParameters.pl -rpkm <RPKM file> -i <count Files diretory> -o <outputfile> [-l listcount]
             
 estimateNBParameters.pl -help

 estimateNBParameters.pl -version

 For help, run this script with the -help option.


=head1 OPTIONS

=head2 -rpkm <RPKM file>

This file contains the RPKM values of the intervals.

=head2 -i <count Files diretory>

This directory contains the count files. The files are separated according to the choromosomes.

=head2 -o <output file>

The output file contains the estimated parameters of the Negative Binomial distribution. 

=head2 -l <listcount>

This parameter is optional. 
It is the number of parameters used in estimating the parameters of the intervals,
Its default value is 1000.
More precisely, intervals are sorted by their RPKM values. 
In this sorted list, for each interval, listcount/2 many intervals above and listcount/2 many intervals below
the interval are taken and 5% from the top and 5% from the bottom, having the highest averagge counts, are removed.
The read counts of the remaining intervals are used to estimate the parameters of the current interval.

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
