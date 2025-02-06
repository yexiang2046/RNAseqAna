#!/usr/bin/env perl

#########################################################################################
#                                       prepInput.pl
#########################################################################################
# 
# This program converts bed files to BedGraph using only the center of the reads.
# It also find all counts for the intervals in region files.
# For detailed information and license, run this script with -help option
# or see the documentation at the bottom.
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
 my $strand           = "";
 my $regionfile       = "";
 my $countfile        = "";
 my $help             = "";
 my $print_version    = "";
 my $version          = "1.0.0";
################### PARAMETER PARSING ####################

my $cmd=$0." ".join(" ",@ARGV); ####command line copy

GetOptions( 
	'input=s'        => \$input,
	'output=s'       => \$output,
	'strand'         => \$strand,
	'count=s'        => \$countfile,
	'regionfile=s'   => \$regionfile,
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
#    Read input file calculate the # of reads centered each nucleotide and write it into output file

 
 open(IN, $input);
 my $chrom=basename($input, ".bed");
  
 my %bgvals=();
 my %diffisoforms=(); #Find all different isoforms for different center point. Divide read count to total isoforms for this position
 my $watcri="a";
 while(my $line=<IN>)
 {   
   chomp($line);   
   my @a=split(/[\t\s]+/,$line );
   my $center=floor(($a[1]+$a[2])/2);
   
   if ($strand) {
      $watcri="w";
      if ($a[5] eq "-") {
	$watcri="c";
      }
   }
   $bgvals{$watcri}{$center}++;
   $diffisoforms{$watcri}{$center}{$a[3]}=1;
 }
 close(IN);

 my @watcriarray=();
 my %totalmapped=();
 if ($strand) {
    push(@watcriarray, "watson");
    push(@watcriarray, "crick");
    $totalmapped{"w"}=getMappedCount($countfile."_watson_mapped.count");
    $totalmapped{"c"}=getMappedCount($countfile."_crick_mapped.count");
    lib::io::makeDir("$output/watson");
    lib::io::makeDir("$output/crick");
    lib::io::makeDir("$output/watsoncounts");
    lib::io::makeDir("$output/crickcounts");
 }
 else
 {
    push(@watcriarray, "all");
    $totalmapped{"a"}=getMappedCount($countfile."_all_mapped.count");
    lib::io::makeDir("$output/all");
    lib::io::makeDir("$output/allcounts");
 }
 
 my %counts=();
 my %squarecounts=();
 
 foreach my $watcri (@watcriarray)
 {
   open(OUT, ">$output/$watcri/$chrom.bg");
   my $wca=substr($watcri,0,1);
   if (exists $bgvals{$wca})
   {
     my %vals=%{$bgvals{$wca}};
     my %diffiso=%{$diffisoforms{$wca}};
     foreach my $pos (sort{$a<=>$b} keys %vals)
     {
       my %diff=%{$diffiso{$pos}};
       my $size = scalar keys %diff;
       my $count = floor($vals{$pos}/$size);
       foreach my $region (keys %diff)
       {
	 #print $region."\n";
	 $counts{$wca}{$region}+=$count;
	 $squarecounts{$wca}{$region}+=$count*$count;
       }
       print OUT "$chrom\t$pos\t".($pos+1)."\t".$count."\n";
     }
   }
   close(OUT);
  }
  undef %bgvals;
  undef %diffisoforms;
  my %countlines=();
  foreach my $watcri (@watcriarray)
  {
    
    my $wca=substr($watcri,0,1);
    if (exists $counts{$wca})
    {
      my %c=%{$counts{$wca}};
      my %sc=%{$squarecounts{$wca}};
      
      foreach my $key (keys %c)
      {
        #print $key."\n";
        $countlines{$wca}{$key}=$c{$key}."\t".$sc{$key};
      }
    }
  }
  undef %counts;
  undef %squarecounts;
  my @region=();
  lib::io::getRegionChrom(\@region, $regionfile, $chrom);
  foreach my $watcri (@watcriarray)
  {
    my $wca=substr($watcri,0,1);
    if (exists $countlines{$wca}) {	
    
    my %clines=%{$countlines{$wca}};
    my $outfile="$output/".$watcri."counts/$chrom.intervalcount";
    open(OUT, ">$outfile");
    foreach my $interval (@region)
    {
	my $count=0;
	my $squarecount=0;
	if (exists $clines{$interval->name})
	{  
	  my @arr=(/\t/,$clines{$interval->name});
          $count=$arr[0];
	  $squarecount=$arr[1];
	}
	
	my $len=$interval->end - $interval->start;
	my $rpkm=0;
	if (($len*$totalmapped{$wca})>0)
	{
	 $rpkm=(10**9*$count)/($len*$totalmapped{$wca});
	}
	print OUT $chrom."\t".$interval->name."\t".$len."\t".$count."\t".$squarecount."\t$rpkm\n";
    }
    close(OUT);
    }
  }
  
sub getMappedCount
{
    my ($file)=@_;
    if (!(-s $file)) {
    	print "Please run StepCountMappedReads!!!";
	exit;
    }
    open (IN, $file);
    my $count=<IN>;
    $count=~/^([\d]+)/;
    $count=$1;
    return $count;
}
__END__


=head1 NAME

prepInput.pl

=head1 SYNOPSIS

prepInput.pl -i input file <BED format> -o output file <BG format>

prepInput.pl -help

prepInput.pl -version

For help, run this script with -help option.

=head1 OPTIONS

=head2 -i  input file <BED format> 

It converts input bed file for each chromosme separately.

=head2 -o output file <BG format>

BedGraph file that is calculated only from read centers

=head2 -help

Display this documentation.

=head2 -version

Display the version

=head1 DESCRIPTION

 This program converts bed files to BedGraph using only the center of the reads.
 The bed files have to be separated to the chromosmes. Otherwise it doesn't work right.


=head1 EXAMPLE

./prepInput.pl -p your_input_dir/chr1.bed -o your_output_dir/chr1.bg 

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
