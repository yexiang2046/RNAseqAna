#!/usr/bin/env perl
#########################################################################################
#                                      findPeaks.pl
#########################################################################################
# 
# This program finds the peaks in a given Rip / Clip Seq data.
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
 use List::Util qw[min max sum];
 use Scalar::Util qw[looks_like_number];
 use POSIX qw/floor/;
 use strict;
 use lib "../";
 require lib::io;
 require lib::func;
 use Getopt::Long;
 use Pod::Usage;
 use Math::CDF;
 use Class::Struct;

#################### CONSTANTS ###################### 
use constant MINTAGCOUNT => 3;
use constant MINTFOOTSIZE => 1;

#use constant SMALLINTERVAL => 20; #
## Set these to 1 to see the related parameters and peakfinding steps.
my $DEBUG = 0;
my $DETECTDEBUG = 0;
my $INDEXDEBUG = 0;
my $NORNASEQDEBUG = 0;
#################### VARIABLES ######################

struct( interval => [ identifier => '$', start => '$', end => '$', strand => '$' , count => '$', sumSquare => '$' , RPKM => '$', p => '$', r => '$' ]);
struct( parameter => [ p => '$', r => '$' ]);

struct( 'Peak', { interval => '$', strand => '$', tagcount => '$', max_height => '$', footsize => '$',
 start => '$', end => '$', intervalLength=> '$', RPKM => '$', weighted_center => '$', block=> '$', pVal => '$',
} );

my $refFile=""; # Hold the annotation coming from bed file

# If no lambda file is given, we will obtain the lambda
# from the bedGraph file  
my $NBparameterFile="";

my $BGfile = "";
my $outdir = "";
my $pcutoff = 0.01;
#my $countsdir = "";
my $RPKMfile = "";
my $chrom = "";

my $frameRadiiString = "";
my $help = "";
my $print_version = "";
my $version = "2.0.0"; 
my $noRNASeq = "";
my $gapNumber = 2;

# Given the $pcutoff value and the lambdavalue,
# in our probability distribution, there willl be a critical count value, say N,
# such that the probability of seeing count greater than N is $pcutoff.
# For a fixed $pcutoff, we store these N values in this hash as:  lambda => criticalCount (namely N)
# This way we avoid computing the same critical cutoff count over and over again.
# We don't expect to have many lambda values so keeping them in the memory is feasible
my %cutoffCounts = ();

## These variables are needed for lambda calculation
## They must be initialized by the function initializeFrameSums
my @frameSums;
my @rightSizes;
my @leftSizes;
my $lastLambdaPosition = -1;
#my @frameSums = ();
my @frameSquareSums = ();

################### PARAMETER PARSING ####################

my $cmd = $0." ".join(" ",@ARGV); ####command line copy 

GetOptions( 
	'r=s'   => \$refFile,
	'l=s'    => \$NBparameterFile, 
	'c=s'   => \$chrom,
	'b=s'   => \$BGfile,
	'o=s'   => \$outdir,
	'p=f'   => \$pcutoff,
	'rpkm=s'  => \$RPKMfile,
	'radii=s'   => \$frameRadiiString,
	'help'  => \$help,
	'version' => \$print_version,
	'noRNASeq' => \$noRNASeq,
	'gap=i' => \$gapNumber,
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

# DEBUG PURPOSES
# The statistical test can be checked with this function
# testNBManually(0.01, 0.05);

pod2usage( {'-verbose' => 0, '-exitval' => 1,} ) 
	  if ( ( (not $refFile) or (not $chrom) or (not $BGfile) or (not $outdir) or (not $RPKMfile) or (not $NBparameterFile) ) 
	    and ( (not $refFile) or (not $chrom) or (not $BGfile) or (not $outdir) or (not $noRNASeq) )
	  );

# Create the output directory and write the log for the peak.
# Note that this is a separate log file for this script.
# It creats "summaryPeak.chr*.log*" files where ** are numbers.
# The main program (findPeaks.pl) itself creates "summary.log.*" files. 
# Note that for each run a separate log file is created.
# This is done in the if block below.

my $logDir = $outdir."/LOGS";
`mkdir -p $outdir`;
`mkdir -p $logDir`;
 my $logname = "summaryPeak.".$chrom.".log";
 my $num = 1;

 #Get the existing log files parse their names using awk and get the number in the tail.
 #Then increment the number and create the new log by labeling it with that number at the end of its name,
 if (-s "$logDir/$logname.1")
 {
   $num = `ls $logDir/$logname.*|awk '{split(\$1, a, "."); print a[3]; }'|sort -n|tail -n 1`;
   $num++;
 }

# These are the default frame radii to compute the lambda value
# if no fram radius and no lambda file is given
my @frameRadii = (1000,10000,50000);
@frameRadii = split("," , $frameRadiiString) if($frameRadiiString);

## Initialize the arrays for lambda calculation
for(my $i=0; $i < scalar(@frameRadii); $i++)
{
  ($frameSums[$i], $frameSquareSums[$i], $rightSizes[$i],  $leftSizes[$i])  = (0, 0, 0, 0);
}

open(my $LOG, ">$logDir/$logname.$num") or die("Can not open the logDir: $logDir/$logname.$num\n");

 lib::io::wlog($LOG, "program started at: ");    
 lib::io::wlog($LOG, "The command used is:");
 lib::io::wlog($LOG, "$cmd\n");
 lib::io::wlog($LOG, "Parameters specified:");
 lib::io::wlog($LOG, "BGfile: $BGfile");
 lib::io::wlog($LOG, "outdir: $outdir");
 lib::io::wlog($LOG, "chrom: $chrom");
 lib::io::wlog($LOG, "refFile: $refFile");
 lib::io::wlog($LOG, "lambdafile: $NBparameterFile");
 lib::io::wlog($LOG, "frameRadii: ".join(",".@frameRadii));
 lib::io::wlog($LOG, "noRNASeq: $noRNASeq");
 lib::io::wlog($LOG, "Gaps Tolerated (gapNumber): $gapNumber");

##################################################################
################### Read Input Files ####################
##################################################################

# Read the contents (actually the first and the second columns) of the $lambdafile
# into %hashLambda
 my %hashNBparameters = ();
 my %hashRPKM = ();

 # If a lambda file is not specified, then the lambda values
 # are stored in this array.
 my @localLambdaValues = ();
 my @expandedCounts = ();
 
   # read the contents of the bed file into intervals
  my @intervals = ();
  my $c = 0;
  lib::io::wlog($LOG, "getIntervals: $chrom");
  getIntervals(\@intervals, $refFile, $chrom);
  lib::io::wlog($LOG, "getGenes:".scalar(@intervals)." DONE\n");

  #printIntervals(\@intervals, "all"); # DEBUG CODE!!!
  
  if(not $noRNASeq)
  {
    lib::io::wlog($LOG, "System is working WITH RNASeq data.");
    lib::io::wlog($LOG, "getRPKM: $chrom");
    getRPKM($RPKMfile, \$c, \%hashRPKM);
    lib::io::wlog($LOG, "getRPKM:$c DONE\n");
    lib::io::wlog($LOG, "get Negative Binomial Parameters: $chrom");
    getNBparameters($NBparameterFile, \%hashNBparameters, \$c);
    lib::io::wlog($LOG, "get Negative Binomial Parameters: $c DONE\n");
  }
  else
  {
    lib::io::wlog($LOG, "System is working WITHOUT RNASeq data.");
    lib::io::wlog($LOG, "No RPKM file is used!"); 
    lib::io::wlog($LOG, "No NB Parameter file is used!");
  }
 
#  DEBUG
#  print "Key -> RPKM\n";
#  foreach my $key (keys %hashRPKM){
#   print "$key -> $hashRPKM{$key}\n";
#  }
  
  my @sortedKeys;
  my %hashtags;
  #Read the bedGraph ile into the 4 variables above  
  lib::io::wlog($LOG, "getBGVals: $BGfile");
  getBGVals($BGfile, \%hashtags, \@sortedKeys, \$c);
  lib::io::wlog($LOG, "getBGVals: $c DONE\n");
  
  # initialize the rpkm values of the intervals, if there is RNASeq data
  if(not $noRNASeq)
  { 
    foreach my $int (@intervals)
    {
      $int->RPKM($hashRPKM{$int->identifier}) if( exists($hashRPKM{$int->identifier}) );
      $int->r($hashNBparameters{$int->identifier}->r) if( exists( ${hashNBparameters}{$int->identifier} ) );
      $int->p($hashNBparameters{$int->identifier}->p) if( exists( ${hashNBparameters}{$int->identifier} ) );
    }
  }


#DEBUG
#  print "RIPSEQ DATA:\n";
#  print "Key -> Count\n";
#  my @bgKeys = sort {$a<=> $b} keys %hashtags;
#  foreach my $key (@bgKeys){
#   print "$key -> $hashtags{$key}\n";
#  }
#print("Sorted Array List Is:\n".join("\n", @sortedKeys)."\n"); 

#DEBUG
#printParameters(\%hashNBparameters);

# printIntervals(\@intervals, "all"); # DEBUG CODE!!!
# exit;
  
##################################################################
################### M A I N     P R O G R A M ####################
##################################################################  

print "MINTAGCOUNT = ",MINTAGCOUNT,"     MINFOOTSIZE = ", MINTFOOTSIZE, "   ALLOWED GAPS = $gapNumber\n" if($DEBUG);

my $noRPKMcount = 0;
# DEBUG
# print("The first position is $sortedKeys[0]\n");
# print("FrameSums:".join(" ", @frameSums)."\n");
# print("FrameSquareSums:".join(" ", @frameSquareSums)."\n");
# exit;

my @peaks = ();
my $i = 0;

open (my $OUTPEAK, ">$outdir/$chrom.peaks") or die("Can not open the file $outdir/$chrom.peaks for writing.\n");
print $OUTPEAK join("\t", "#Chr", "Strand", "Interval", "Interval_Length" , "Tag_Count", "Maximum_Height", "Foot_Size", "Start", "End", "RPKM", "Weighted_Center", "p-Value"  ), "\n";
my $oldi=0;
my $oldintervalend=0;
foreach my $interval (@intervals)
{
  if ($oldintervalend>$interval->start) {
    $i=$oldi;
  }
  $oldintervalend=$interval->end;  
  $oldi=$i;
 if($noRNASeq)
 {  

    # More parameters here later
    
    detectPeaksNORNASEQ(\$i, \@peaks, \%hashtags, \@sortedKeys, $interval, $pcutoff, $gapNumber, \@frameRadii);
    writePeaks(\@peaks, $OUTPEAK, $chrom);
    @peaks = (); # We have written the peaks, so flush the array.
    $noRPKMcount++;
  
 }
 else
 {
    my ($start, $end) = ($interval->start,$interval->end); 
    print "interval=",$interval->identifier, ,"[$start , $end)   i = $i,   pos = $sortedKeys[$i] ","\n" if ($INDEXDEBUG);
    if($interval->RPKM > 0)
    {
      # We have RNASeq data. So the NB parameters should come from the parameters file.
      # Note that those values are already written into the interval structs p and r attributes.   
      detectPeaks(\$i, \@peaks, \%hashtags, \@sortedKeys, $interval, $pcutoff, $gapNumber, \@frameRadii);
      $noRPKMcount++;
    }
    else
    {
      # We don't have RNASeq data for this interval.
      # So we estimate the parameters locally.
      print("Couldn't find RNASeq data for ", $interval->identifier,", rpkm = ", $interval->RPKM ," calculating parameters locally\n") if($NORNASEQDEBUG);
      detectPeaksNORNASEQ(\$i, \@peaks, \%hashtags, \@sortedKeys, $interval, $pcutoff, $gapNumber, \@frameRadii);
    }
        
    writePeaks(\@peaks, $OUTPEAK, $chrom);
    @peaks = (); # We have written the peaks, so flush the array.
  } # end of foreach
}

lib::io::wlog($LOG, "There were $noRPKMcount intervals without RNASeq data. Their NB parameters were determined locally from the interval.\n");
lib::io::wlog($LOG, "findPeaks.pl done.");

# foreach my $interval (@intervals)
# {  
#   for(my $pos = $interval->start; $pos < $interval->end; $pos++)
#   {
#     my ($p, $r) = getParameters(\@frameRadii, \@frameSums, \@frameSquareSums, \@rightSizes, \@leftSizes, \%hashtags, $sortedKeys[0],  $sortedKeys[-1], $pos, \$lastLambdaPosition);
#     print("p = $p   r = $r\n");
#   }
# }

close($OUTPEAK);

############################################################################
################ E N D    O F    M A I N     P R O G R A M #################
############################################################################
 
########################################################################
######################## F U N C T I O N S #############################
########################################################################

########################################################################
#######function :   d e t e c t P e a k s N O R N A S E Q   ############
########################################################################

# If the interval is short, the bacground is computed from the interval.
# If the interval is long enough so that the largest frame can fit into it,
# for each nucleotide position, the background is computed using the counts
# in the frame around the current position. For the next position, the frame is shifted
# to the right.

sub detectPeaksNORNASEQ
{
 my ($i, $peaks, $hashtags, $sarr, $interval, $pcutoff, $gapNumber, $frameRadii) = @_;
 
  print "detectPeaksNORNASEQ: interval = ".$interval->identifier."\n" if ($NORNASEQDEBUG);          
  
 if( ($interval->end - $interval->start) < 2*max( @{$frameRadii} ) )
 {
    # So the interval is short. No need to use local frames.
    my @countsBucket = (); # holds the read counts of the interval
    
    #put all the read counts of the interval into this bucket
    for(my $y = $interval->start; $y < $interval->end; $y++)
    {
      push(@countsBucket, ${$hashtags}{$y} ) if(exists(${$hashtags}{$y}));  
    }
    ## Now give this bucket to the function that estimates the parameters and set the NB parameters of the interval.
    my ($p,$r) = estimateParameters(\@countsBucket, $interval->end - $interval->start);
    $interval->p($p);
    $interval->r($r); 
    ## Now that we know the parameters p and r, call the peaks on this interval as usual.
    print "The whole interval is used as the frame. p = $p , r = $r \n" if ($NORNASEQDEBUG) ;
    detectPeaks($i, $peaks, $hashtags, $sarr, $interval, $pcutoff, $gapNumber);
 }
 else
 {
    # Interval is long.
    # Use local frames to estimate the parameter.
    # This time, the parameters p and r are pcalculated for eac position.
    # Around each position, we consider a frame  [ radius many  positions  , current poistion , radius many positions]
    # Totaling 2*radius+1 nucleotides.
    # Using these nucleotides, we estimate the parameters p and r.
    # For the next position, we shift this frame one unit to the right.
     my $index = ${$i};
     my $endPos = $interval->end;
     my $len =  $interval->end - $interval->start; 
     my $RPKM = $interval->RPKM;
     my ($pos, $pos_old ,$report_footsize, $report_tagcount, $report_max_height, $report_start, $report_end)  = (0,0,0,0,0,0,0);
     my $lookingForNewPeak = 1;
     
    while( ($index < scalar(@{$sarr}) ) and (${$sarr}[$index] < $interval->start) )
    {
      $index++;
    }

    $pos = ${$sarr}[$index] ;
 
    my ($peakFound , $pVal)  = (0 , 1);
    # Continuity test variables
    my $lookAheadPosition = $pos;
    my $continuityTest = 0;
    my $testIndex = $index;
    
    my @frameSums = map {0} @{$frameRadii};
    my @frameSquareSums = map {0} @{$frameRadii};
    my @leftPositions; 
    my @rightPositions;
    
    # In the previous if block we made sure that the initial frame can fit into the interval for all radii 
    # as  ($interval->end - $interval->start) < 2*max( @{$frameRadii} )
    initializeFrameSums ($frameRadii, \@frameSums, \@frameSquareSums, $hashtags, $interval->start, \@rightPositions, \@leftPositions);
    print $interval->identifier, " :\n" if($NORNASEQDEBUG);
    print("Initialized sums are : ", join(", ",@frameSums), "\n") if($NORNASEQDEBUG);
    print("Initialized square sums are : ", join(", ",@frameSquareSums), "\n") if($NORNASEQDEBUG);
    
     my $lastSpikePosition = -$gapNumber;  
     my $lastPos = $pos;
     
     while($pos > 0 && $pos < $endPos && $index < scalar(@{$sarr}))
    {
	$pos = ${$sarr}[$index];
	my $tagcount = ${$hashtags}{$pos};
	
	my $lookAheadPosition = $pos;
	my $continuityTest = 0;
	my $testIndex = $index;
	
	for(my $tempPos = $lastPos; $tempPos <= $pos; $tempPos++)
	{
	    shiftFrameSumsToTheRight($frameRadii, \@frameSums, \@frameSquareSums, $hashtags, $endPos , \@rightPositions, \@leftPositions, $tempPos );
	}
	$lastPos = $pos;
	
	my ($p, $r) = getLocalParameters(\@frameSums, \@frameSquareSums, $frameRadii);
	
	print ("NORNASEQ while loop:\n")  if($NORNASEQDEBUG); 
	print ("pos = $pos, index = $index, count = $tagcount, p = $p , r = $r\n") if($NORNASEQDEBUG) ;
	print("Sums are : ", join(", ",@frameSums), "\n") if($NORNASEQDEBUG);
        print("Square sums are : ", join(", ",@frameSquareSums), "\n") if($NORNASEQDEBUG);
    
	 my @tempSums = @frameSums;
	 my @tempSumSquares =  @frameSquareSums;
	 my @tempLeftPositions =  @leftPositions;
	 my @tempRightPositions =  @rightPositions;
	 my $lastTempPos = $lookAheadPosition;
	 my $tempPVal = $pVal;
	 
	CONTTESTLOOPNORNASEQ: while( ( ($lookAheadPosition - $lastSpikePosition)  <=  ($gapNumber + 1)  ) and ($lookAheadPosition < $endPos) and ($testIndex < scalar(@{$sarr}) ) )
	{ 
	  for(my $tempPos = $lastTempPos; $tempPos <= $lookAheadPosition; $tempPos++)
	  {
	    shiftFrameSumsToTheRight($frameRadii, \@tempSums, \@tempSumSquares, $hashtags, $endPos , \@tempRightPositions, \@tempLeftPositions, $tempPos );
	  }
	  
	  $lastTempPos = $lookAheadPosition;
	  my ($pT, $rT) = getLocalParameters(\@tempSums, \@tempSumSquares, $frameRadii);
	  print "CONTTESTLOOPNORNASEQ: p = $pT , r = $rT\n" if ($NORNASEQDEBUG);
	  
	  my $lookAheadpVal = NegativeBinomial(${$hashtags}{$lookAheadPosition}, $rT, $pT, 1);
	  # We are adding one new position to our current peak (for testing) so increment the summ by the current hashtag and increment the span size by 1 to get the pvalue
	  $tempPVal = NegativeBinomial( $report_tagcount + ${$hashtags}{$lookAheadPosition}, $rT * ($report_end - $report_start + 2) , $pT, 1);
	  
	  if(  ( $lookAheadpVal < $pcutoff) 
		and ( $tempPVal < $pcutoff )  
		and ${$hashtags}{$lookAheadPosition} > MINTAGCOUNT )	  
	  {
	    $continuityTest = 1;
	    last CONTTESTLOOPNORNASEQ;
	  }
	  
	  $lookAheadPosition = ${$sarr}[++$testIndex];	 
	}# end of CONTTESTLOOPNORNASEQ
    
    # If this position is the first position of a possible peak then set $continuityTest to 1.
    # This can be the first position of a possible peak if and only if $lastSpikePosition = -$gapnumber and the read is significantly large.
      if( ($lastSpikePosition == -$gapNumber) and (NegativeBinomial($tagcount, $r, $p, 1) < $pcutoff && $tagcount > MINTAGCOUNT))
      {
	$continuityTest = 1;
      }
    
      ### If this read is significantly large, update the last spike position.
      if( ( $tagcount > MINTAGCOUNT ) and ( NegativeBinomial($tagcount, $r, $p, 1) < $pcutoff ) )
      {
	print "possible peak at $pos\n" if $DETECTDEBUG;
	$lastSpikePosition = $pos;	      
      }
    
      if( $lookingForNewPeak )
      {
	#so we are looking for a new peak           
	if($continuityTest)
	{
	  $report_footsize++;
	
	  if($report_footsize >=  MINTFOOTSIZE)
	  {
	    print "Found a peak. Trying to see if it extends.\n" if($DETECTDEBUG);
	    # Ok, we have sufficiently large counts in sufficiently consecutive places
	    # so this must be a peak!
	    $lookingForNewPeak = 0;
	    $report_start = $pos - (MINTFOOTSIZE) + 1;
	  }	
	  $report_tagcount += $tagcount;
          $report_end = $pos;      
          $report_max_height = $tagcount if ($report_max_height < $tagcount); 
          ## $pVal = NegativeBinomial($report_tagcount, $r * ($report_footsize), $p, 1); # To be deleted
	}
	else
	{
	   #We dont have continuity. 
	   #So we prepare for a new peak by resetting these values to 0. 
	   #So reset the values to 0.
	   ($report_footsize, $report_tagcount,$report_max_height, $report_start, $report_end) = (0, 0, 0, 0, 0);
	   $lastSpikePosition = -$gapNumber;
	   $pVal = 1;
        }
    }
    else
    {
      # we already found a peak we will see if it expands or ends  
	# if the tag counts are large, expand it
	#if( NegativeBinomial($tagcount, $r, $p, 1) < $pcutoff && ($pos == ($pos_old + 1) ) && $tagcount > MINTAGCOUNT )
	if($continuityTest)
	{
	  $report_footsize++;
	  $report_tagcount += $tagcount;
	  $report_end = $pos;
	  $report_max_height = $tagcount if ($report_max_height < $tagcount);
	  $lastSpikePosition = $pos;
	  # $pVal *= NegativeBinomial($tagcount, $r, $p, 1); # to be deleted
	}
	else
	{
	  $pVal = NegativeBinomial($report_tagcount, $r * ($report_end - $report_start + 1), $p, 1);
	  print "Saving peak on $report_start - $report_end, pval = $pVal\n" if($DETECTDEBUG);
	  #tagcount is small, so this is the end of the peak
	  savePeak($peaks, $hashtags, $interval, $p, $r, $report_start, $report_end,$report_max_height, $report_tagcount, $report_footsize, $pcutoff, $pVal);     
	  $lookingForNewPeak = 1;
	  ($report_footsize, $report_tagcount, $report_max_height, $report_start, $report_end) = (0 , 0 , 0 , 0 , 0);
	  $lastSpikePosition = -$gapNumber; 
	}      
    }# end of the else of  if( $lookingForNewPeak )
    
    $pos_old = $pos; 
    $index++;
    $pos = ${$sarr}[$index];
   } #end while($pos > 0 && $pos <= $endPos && $index < @{$sarr})

   my ($p, $r) = getLocalParameters(\@frameSums, \@frameSquareSums, $frameRadii);
   
   if ( not $lookingForNewPeak)
   {
    # We have one last peak to save
    $pVal = NegativeBinomial($report_tagcount, $r * ($report_end - $report_start + 1), $p, 1);
    savePeak($peaks, $hashtags, $interval, $p, $r, $report_start, $report_end,$report_max_height, $report_tagcount, $report_footsize, $pcutoff, $pVal);  
    print "Saving peak on $report_start - $report_end\n" if($DETECTDEBUG);
   } 
 ${$i} = $index;
    
 }# main else block if( ($interval->end - $interval->start) < 2*max( @{$frameRadii} ) )
 
} 



###############################################################################
########### Function: d e t e c t P e a k s   #################################
###############################################################################
#This function is callled for each interval and it finds the peaks in each interval.
#It pushes found peaks to the peaks array.

# Our statistical test is implemeted in the function NegativeBinomial. 
# If its return value is less than $pcutoff, we call the read under consideration as significant or significantly large.

sub detectPeaks
{
 my ($i, $peaks, $hashtags, $sarr, $interval, $pcutoff, $gapNumber) = @_;
 my $index = ${$i};
 my $endPos = $interval->end;
 my $len =  $interval->end - $interval->start; 
 my $RPKM = $interval->RPKM;
 my ($pos, $pos_old ,$report_footsize, $report_tagcount, $report_max_height, $report_start, $report_end)  = (0,0,0,0,0,0,0);
 my $lookingForNewPeak = 1;
 #print ("Beginning index = $index, interval = ".$interval->identifier."\n");
 # If current position (namely ${$sarr}[$index]) les than the start position of the interval,
 # Then this read (this position) is on a further interval. So 
  # move the index forward till the position 
  while( ($index < scalar(@{$sarr}) ) and (${$sarr}[$index] < $interval->start) )
  {
    $index++;
  }

  $pos = ${$sarr}[$index] ;
 
 my $peakFound = 0;
 my $pVal = 1;
 
#  print "interval\tr\tp\n" if($DEBUG);
#  print $interval->identifier,"\t",$interval->r, "\t", $interval->p,"\n" if($DEBUG);

  my ($r, $p) = ($interval->r, $interval->p);
  my $lastSpikePosition = -$gapNumber;
  
 while($pos > 0 && $pos < $endPos && $index < scalar(@{$sarr}))
 {
    $pos = ${$sarr}[$index];
    my $tagcount = ${$hashtags}{$pos};
    print ("Beginning of while loop: pos = $pos, index = $index, count = $tagcount, interval = ".$interval->identifier."\n") if($DETECTDEBUG) ;
   
    ### We allow to be gapNumber many gaps between the positions.
    ### We test this here. We call this the continuity test.
    ### Two reads are said to be continuous if they both  pass the NB statistical test and
    ### the diference between their positions indices is less than or equal to the gapnumber 
    ### The last peak candidate is in $lastSpikePosition.
    ### If we haven't seen any peak candidates before or we saved them in the previous peaks
    ### We set $lookAheadPosition to -gapNumber to indicate the beginning of the new peaak.
    ### This way we can make the first position of the peak "continuous" on its own.
    my $lookAheadPosition = $pos;
    my $continuityTest = 0;
    my $testIndex = $index;
    my $tempPVal = $pVal;
    
    CONTTESTLOOP: while( ( ($lookAheadPosition - $lastSpikePosition)  <=  ($gapNumber + 1)  ) and ($lookAheadPosition < $endPos) and ($testIndex < scalar(@{$sarr}) ) )
    { 	
      my $lookAheadpVal = NegativeBinomial(${$hashtags}{$lookAheadPosition}, $r, $p, 1);
      # note we test  what happens if the tagcount  at the lookaheadPosition is added to the current peak
      # so we increment the report tagcount by ${$hashtags}{$lookAheadPosition} and increment the peak span size  by 1 which is 1 + ($report_end - $report_start + 1)
      $tempPVal = NegativeBinomial($report_tagcount + ${$hashtags}{$lookAheadPosition}, $r * ($report_end - $report_start + 2) , $p, 1);
      
      if(  ($lookAheadpVal < $pcutoff) 
	    and ( $tempPVal < $pcutoff )
	    and ${$hashtags}{$lookAheadPosition} > MINTAGCOUNT )
      {
	$continuityTest = 1;
	last CONTTESTLOOP;
      }
      $lookAheadPosition = ${$sarr}[++$testIndex];
    }
    
    # If this position is the first position of a possible peak then set $continuityTest to 1.
    # This can be the first position of a possible peak if and only if $lastSpikePosition = -$gapnumber and the read is significantly large.
    if( ($lastSpikePosition == -$gapNumber) and (NegativeBinomial($tagcount, $r, $p, 1) < $pcutoff && $tagcount > MINTAGCOUNT) )
    {
      $continuityTest = 1;
    }
    
    ### If this read is significantly large, update the last spike position.
    if(  ($tagcount > MINTAGCOUNT) and  ( NegativeBinomial($tagcount, $r, $p, 1) < $pcutoff ) )
    {
      print "possible peak at $pos\n" if $DETECTDEBUG;
      $lastSpikePosition = $pos;	      
    }
    
    if( $lookingForNewPeak )
    {
      #so we are looking for a new peak           
      if($continuityTest)
      {
	$report_footsize++;
	
	if($report_footsize >=  MINTFOOTSIZE)
	{
	  print "Found a peak. Trying to see if it extends.\n" if($DETECTDEBUG);
	  # Ok, we have sufficiently large counts in sufficiently consecutive places
	  # so this must be a peak!
	  $lookingForNewPeak = 0;
	  $report_start = $pos - (MINTFOOTSIZE) + 1;
	}	
	  $report_tagcount += $tagcount;
          $report_end = $pos;      
          $report_max_height = $tagcount if ($report_max_height < $tagcount); 
      }
      else
      {
	   #We dont have 
	   #So we repare for a new peak by resetting these values to 0. 
	   #So reset the values to 0.
	   $report_footsize = 0;
	   $report_tagcount = 0;
	   $report_max_height = 0;
	   $report_start = 0;
	   $report_end = 0;
	   $lastSpikePosition = -$gapNumber;
	   $pVal = 1;
      }
    }
    else
    {
      # we already found a peak
      # we will see if it expands or ends  
	# if the tag counts are large, expand it
	#if( NegativeBinomial($tagcount, $r, $p, 1) < $pcutoff && ($pos == ($pos_old + 1) ) && $tagcount > MINTAGCOUNT )
	if($continuityTest)
	{
	  $report_footsize++;
	  $report_tagcount += $tagcount;
	  $report_end = $pos;
	  $report_max_height = $tagcount if ($report_max_height < $tagcount);
	  $lastSpikePosition = $pos;
	}
	else
	{
	  $pVal = NegativeBinomial($report_tagcount, $r * ($report_end - $report_start + 1), $p, 1);
	  print "Saving peak on $report_start - $report_end, pVal = $pVal\n" if($DETECTDEBUG);
	  #tagcount is small, so this is the end of the peak
	  savePeak($peaks, $hashtags, $interval, $p, $r, $report_start, $report_end,$report_max_height, $report_tagcount, $report_footsize, $pcutoff, $pVal);     
	  $lookingForNewPeak = 1;
	  ($report_footsize, $report_tagcount, $report_max_height, $report_start, $report_end) = (0 , 0 , 0 , 0 , 0);
	  $lastSpikePosition = -$gapNumber; 
	  $pVal = 1; # Reset the pVal for the peaks to come.
	}      
    }
    
  $pos_old = $pos; 
  $index++;
  $pos = ${$sarr}[$index];
 } #end while($pos > 0 && $pos <= $endPos && $index < @{$sarr})

 if ( not $lookingForNewPeak)
 {
   # We have one last peak to save 
   $pVal = NegativeBinomial($report_tagcount, $r * ($report_end - $report_start + 1), $p, 1);
   savePeak($peaks, $hashtags, $interval, $p, $r, $report_start, $report_end,$report_max_height, $report_tagcount, $report_footsize, $pcutoff, $pVal);  
   print "Saving peak on $report_start - $report_end\n" if($DETECTDEBUG);
   $pVal = 1; # Reset the pVal for the peaks to come.
 } 
 ${$i} = $index;
}

###############################################################################
########### Function: s a v e P e a k    ######################################
###############################################################################
sub savePeak
{
my ($peaks, $hashtags, $interval, $p, $r, $report_start, $report_end,$report_max_height, $report_tagcount, $report_footsize, $pcutoff, $pVal) = @_;
    
    #my $pVal = NegativeBinomial($report_tagcount, $r,$p, $report_footsize);

    if($pVal < $pcutoff)
    {
     my $weighted_center = 0;  
     if($interval->strand eq '+')
     {
        my $tot_tag_c = 0;
        my $k = $report_start;     
        while($tot_tag_c <= ($report_tagcount/2))
        {
          $tot_tag_c += ${$hashtags}{$k};
          $weighted_center = $k;
          $k++;
        }            
     }
     else
     {
         my $tot_tag_c = 0;
         my $k = $report_end;     
         while($tot_tag_c <= ($report_tagcount / 2))
         {
           $tot_tag_c += ${$hashtags}{$k};
           $weighted_center = $k;
           $k--;
         }               
     }
     
     # NOTE Just like in bedGraph files, peak->start is inclusive and peak->end is exclusive
     # So we add 1 to report_end
     my $peak = Peak->new( tagcount => $report_tagcount, start => $report_start, end => $report_end + 1, max_height => $report_max_height , weighted_center => $weighted_center, strand => $interval->strand );
     $peak->footsize($report_footsize);
     my $block = ($peak->footsize == 1) ? 1 : 0;
     $pVal=sprintf("%0.2e", $pVal);
     $pVal = 1e-45 if($pVal <= 1e-45);
     $peak->pVal($pVal);
     $peak->RPKM($interval->RPKM);
     $peak->block($block);
     $peak->interval($interval->identifier);
     $peak->intervalLength($interval->end - $interval->start);
     push(@{$peaks}, $peak);
    }
}

################################################################################
#######   getIntervals
################################################################################
sub getIntervals
{
  my($intervalArray, $refFile, $chromosome)= @_;
  # This struct is defined above. IntervalArray holds these structs 
  #struct( interval => [ identifier => '$', start => '$', end => '$', strand => '$' , count => '$', sumSquare => '$' , RPKM => '$' ]);
  open(my $FH, $refFile) or die("Error: getIntervals(): Couldn't open the bed file $refFile\n");
  
  GETINTLOOP: while(<$FH>)
  {
    chomp($_);
    my @contents = split(/\s+/, $_);
    next GETINTLOOP if(scalar(@contents) < 6); # just in case, verifying that bed file has sufficint columns
    if($contents[0] eq $chromosome)
    {
      my $currentInterval = interval->new( identifier => $contents[3], start => $contents[1], end => $contents[2], strand => $contents[5], count => 0, sumSquare => 0, RPKM => -1, r => 0, p=> 0.5 );
      push(@{$intervalArray}, $currentInterval);
    }
  }
  
  close($FH);
}


###############################################################################
##### Function: Negative Binomial     ###################################
###############################################################################

# Returns the probability of observing k or more reads in the NB distribution.

# Note that pnbinom(x,r,p) returns the prob of observing x or less reads.
# Therefore 1 - pnbinom(x,r,p) returns the probability of observing reads STRICTLY greater than x.
# We want to find the robability of observing x or more reads.
# So we subtract 1 from x, that is we use 1 - pnbinom(x-1,r,p)
# This gives us the probability of  observing reads strictly greater than x-1
# In other words the probability of observing x or more reads.
sub NegativeBinomial{
  my($k, $r, $p, $windowsize) = @_;
  
  return 1 if($k <= 0); # pValue of 0 is 1.
  my $tag_c_window = ($k-1) / $windowsize;
  
  return(1 - Math::CDF::pnbinom($tag_c_window, $r, $p) );  
}


###############################################################################
##### Function: getLambda  ##############################
###############################################################################
# returns the lambda value of the given $position
# make calls to this function in ascending order (with respect to position)
#
sub getParameters
{
  my ($frameRadii, $frameSums, $frameSquareSums, $rightSizes, $leftSizes, $hashtags, $indexStartPosition, $indexEndPosition,$position, $lastPosition, $rightEndPointers) = @_;
  
  die("getLambda: Error!\n\$position(= $position) is less than \$lastPosition(=${$lastPosition}).  Aborting\n") if($position < ${$lastPosition});
  
  for(my $i = 1; $i <= $position - ${$lastPosition}; $i++)
  {
    shiftFrameSums($frameRadii, $frameSums,$frameSquareSums, $rightSizes, $leftSizes, $hashtags, $indexStartPosition, $indexEndPosition, ${$lastPosition} + $i, $rightEndPointers);
  }
  
  ${$lastPosition} = $position;
  
  my @averages = ();
  my ($maxAvg, $maxIndex) = (0,0);
  
  for(my $k = 0; $k < scalar(@{$frameSums}); $k++)
  {
    $averages[$k] = ${$frameSums}[$k] / (${$rightSizes}[$k] + ${$leftSizes}[$k] + 1);
    if($averages[$k] >= $maxAvg){
      $maxAvg = $averages[$k];
      $maxIndex = $k;
    }
  }
  
  my $moment1 = ${$frameSums}[$maxIndex] / (${$rightSizes}[$maxIndex] + ${$leftSizes}[$maxIndex] + 1);
  my $moment2 = ${$frameSquareSums}[$maxIndex] / (${$rightSizes}[$maxIndex] + ${$leftSizes}[$maxIndex] + 1);
  
  my $p = $moment1 / ($moment2 - ($moment1**2) );
  $p = 0.5 if ( ($p < 0.00001) or ($p > 0.99999) );
  my $r = ($p * $moment1) / (1 - $p);
  
  return ($p, $r);
}


###############################################################################
##### Function: initializeFrameSums  ##############################
###############################################################################
sub initializeFrameSums
{
  my ($frameRadii, $frameSums, $frameSquareSums, $hashtags, $startPosition, $rightPositions, $leftPositions) = @_;
  
  for(my $i = 0; $i < scalar(@{$frameRadii}); $i++)
  {
      for(my $j = 0; $j < 1 + 2*${$frameRadii}[$i]; $j++ )
      { 
	  ${$frameSums}[$i] += ${$hashtags}{ $startPosition + $j } if (exists(${$hashtags}{ $startPosition + $j }));
	  ${$frameSquareSums}[$i] += ${$hashtags}{ $startPosition + $j }**2 if (exists(${$hashtags}{ $startPosition + $j })); 
      }
  }
  #Left and right positions are inclusive!!!
  @{$leftPositions}= map { $startPosition } @{$frameRadii};
  @{$rightPositions}= map { $startPosition + 2*$_ } @{$frameRadii};
}

###############################################################################
##### Function: shiftFrameSums  ##############################
###############################################################################
# IMPORTANT: This function is a helper for get lambda and it should only be called by  getlamda
sub shiftFrameSumsToTheRight
{
  my ($frameRadii, $frameSums, $frameSquareSums, $hashtags, $endPosition , $rightPositions, $leftPositions, $currentPosition ) = @_;
  #Left and right positions are both inclusive!!!
  for(my $i = 0; $i < scalar(@{$frameRadii}); $i++)
  {
    if( ( ($currentPosition - ${$leftPositions}[$i]) > ${$frameRadii}[$i] )   and  (( ${$rightPositions}[$i] + 1) < $endPosition ) )
    {
      if(exists(${$hashtags}{ ${$leftPositions}[$i] }) )
      {
	${$frameSums}[$i] -= ${$hashtags}{ ${$leftPositions}[$i] };
	${$frameSquareSums}[$i] -= ${$hashtags}{ ${$leftPositions}[$i] }**2; 
      }
      
      if(exists(${$hashtags}{ ${$rightPositions}[$i] + 1}) )
      {
	${$frameSums}[$i] += ${$hashtags}{ ${$rightPositions}[$i] + 1 }; 
	${$frameSquareSums}[$i] += ${$hashtags}{ ${$rightPositions}[$i] + 1 }**2; 
      }
      
      ${$leftPositions}[$i]++;
      ${$rightPositions}[$i]++;
    }
  }
}

# sub shiftFrameSumsToTheLeft
# {
#   my ($frameRadii, $frameSums, $frameSquareSums, $hashtags, $startPosition , $rightPositions, $leftPositions, $currentPosition ) = @_;
#   #Left and right positions are inclusive!!!
#   for(my $i = 0; $i < scalar(@{$frameRadii}); $i++)
#   {
#     if( ( ($currentPosition - ${$rightPositions}[$i]) > ${$frameRadii}[$i] )  and ( ${$leftPositions}[$i] - 1) > $startPosition )
#     {
#       ${$frameSums}[$i] += ${$hashtags}{ ${$leftPositions}[$i] -1} if(exists(${$hashtags}{ ${$leftPositions}[$i] -1 }) );
#       ${$frameSums}[$i] -= ${$hashtags}{ ${$rightPositions}[$i] } if(exists(${$hashtags}{ ${$rightPositions}[$i] }) );
#       ${$leftPositions}[$i]--;
#       ${$rightPositions}[$i]--;
#     }
#   }
# 
# }


###############################################################################
##### Function: getRPKM ##############################
###############################################################################
sub getRPKM{
  my($file, $count, $RPKMhash)= @_;
  my $i = 0;
  open(my $FH, $file) or die("ERROR: getRPKM(): Couldn't open the file $file");
  
  while(<$FH>)
  {
    chomp($_);
    my @contents = split(/\s+/, $_);
    if(scalar(@contents) >= 3){
      $i++;
      ${$RPKMhash}{$contents[1]} = $contents[2]
    }
  }

  ${$count} = $i;
  close($FH);
}



###############################################################################
##### Function: getRPKM ##############################
###############################################################################
#  Each row of the bedGraph file is as
#  <choromosome number> <start position of the match> <end position of the match + 1> <#reads (or matchs)>

sub getBGVals{
  my($file, $hashTags, $sortedKeys, $count) = @_;
  
  my $i = 0;
  open(my $FH, $file) or die("ERROR: getBGVals(): Couldn't open the file $file");
  while(<$FH>)
  {
    chomp($_);
    my @contents = split(/\s+/, $_);
    if(scalar(@contents) >= 4)
    {
      ${$hashTags}{$contents[1]} = $contents[3];
      $i++;
    }
  }
  @{$sortedKeys} = sort {$a<=>$b} (keys %{$hashTags});
  ${$count} = $i;
  close($FH); 
}


###############################################################################
##### Function: getNBparameters ##############################
###############################################################################
# NB PARAMETER FILE FORMAT
#Interval_name  p_parameter r_parameter RPKM Interval_Count Interval_Length

sub getNBparameters{
  my ($file, $hashParameters, $count) = @_;
  
  open(my $FH, $file) or die("Error: getNBparameters(): Couldn't open the file $file");
  my $i = 0;
  
  while(<$FH>)
  {
    chomp($_);
    my @contents = split(/\s+/, $_);
    if(scalar(@contents) >= 3){
      my $currentParameter = parameter->new( p => $contents[1], r => $contents[2]);
      ${$hashParameters}{$contents[0]} = $currentParameter;
      $i++;
    }
  }
  
  ${$count} = $i;
  close($FH);
}


###############################################################################
##### Function: i s I n I n t e r v a l s        ##############################
###############################################################################
# Returns the interval identifier if the given position is in any of the given intervals.
sub isInIntervals
{
  my ($position, $arrayIntervals) = @_;
  
  # This is the EASY way.
  # We can improve the performence using binarry search later.
  foreach my $interval (@{$arrayIntervals})
  {
    return $interval->identifier if( ($position >= $interval->start) and ($position < $interval->end)  );
  }
  return "";  
}

###############################################################################
##### Function: w r i t e P e a k s             ##############################
###############################################################################
## Output the found eaks to the given filehandle $FH

sub writePeaks
{
  my ($peaks, $FH, $chromosome) = @_;
  
  foreach my $peak (@{$peaks})
  {
    print $FH join("\t", $chromosome, $peak->strand, $peak->interval, $peak->intervalLength , $peak->tagcount, $peak->max_height, $peak->footsize, $peak->start, $peak->end, $peak->RPKM, $peak->weighted_center, $peak->pVal  ), "\n"; 
  }
}


###############################################################################
##### Function: FindParameters(Slow)       ####################################
###############################################################################

# given an array and 
sub estimateParameters
{
  my ($countArray, $numberOfTerms) = @_;
  
  my ($sum, $squareSum) = (0 , 0);
  
  foreach my $elt (@{$countArray})
  {
    $sum += $elt;
    $squareSum += ($elt*$elt);
  }
  
  return calculateParameters($sum, $squareSum, $numberOfTerms); 
}

sub calculateParameters
{
  my ($sum, $squareSum, $numberOfTerms) = @_;
  
  $numberOfTerms = 1 if($numberOfTerms == 0);
  
  my $moment1 = $sum / $numberOfTerms;
  my $moment2 = $squareSum / $numberOfTerms;
  
  my $p = 0.5;
  $p = $moment1 / ($moment2  - $moment1**2) if( ($moment2  - $moment1**2) != 0 );
  $p = 0.0000001 if( $p < 0.0000001 );
  $p = 0.9999999 if( $p > 0.9999999 );
  my $r = ($p * $moment1) / (1 - $p);
  
  return($p, $r);
}



sub getLocalParameters
{
  my ($sums, $squareSums, $frameRadii) = @_;
  
  my ($p , $r) = calculateParameters( ${$sums}[0], ${$squareSums}[0], ${$frameRadii}[0]);
  my ($returnP, $returnR) = ($p, $r);
  
  my $maxExpected = $r*(1 - $p) / $p;
  my $maxExpectedIndex = 0;
  
  #Maximize the expected value 
  #That is return the parameters p,r coming from the frames that gives is the maximum expected value.
  for(my $i = 1; $i < scalar(@{$frameRadii}) ; $i++ )
  {
    ($p , $r) = calculateParameters( ${$sums}[$i], ${$squareSums}[$i], ${$frameRadii}[$i]);
    if($maxExpected <  $r*(1 - $p) / $p)
    {
      $maxExpected =  $r*(1 - $p) / $p;
      $maxExpectedIndex = $i;
      ($returnP, $returnR) = ($p, $r);
    }
  }
  
  return ($returnP, $returnR);
}




# sub getLocalParameters
# {
#   my ($sums, $squareSums, $frameRadii) = @_;
#   
#   my $maxExpected =  ( ${$sums}[0]  / (1 + 2*${$frameRadii}[0]) );
#   my $maxExpectedIndex = 0;
#   
#   for(my $i = 1; $i < scalar(@{$frameRadii}) ; $i++ )
#   {
#     if ( $maxExpected > ( ${$sums}[0]  / (1 + 2*${$frameRadii}[0]) ) )
#     {
#       $maxExpected =  ( ${$sums}[$i]  / (1 + 2*${$frameRadii}[$i]) );
#       $maxExpectedIndex = $i;
#     }
#   }
#   
#   return calculateParameters( ${$sums}[$maxExpectedIndex], ${$squareSums}[$maxExpectedIndex], ${$frameRadii}[$maxExpected]);
# }

#### DEBUG FUNCTIONS #############################################################################

  # This struct is defined above. IntervalArray holds these structs 
  #struct( interval => [ identifier => '$', start => '$', end => '$', strand => '$' , count => '$', sumSquare => '$' , RPKM => '$' ]);                       
### DEBUG PURPOSSES
sub printIntervals
{
  my ($arrayIntervals, $ind) = @_;
  
  print "Interal List\n", "--------------------------------------------------------------------------\n";
  print "Identifier \t start \t end \t p \t r \t count \t sumSquare \t strand \t RPKM\n";
  
  if($ind eq "all")
  {
    foreach my $elt (@{$arrayIntervals})
    {
      print $elt->identifier,"\t",$elt->start,"\t", $elt->end, "\t", $elt->p, "\t", $elt->r , "\t", $elt->count, "\t", $elt->sumSquare, "\t", $elt->strand, "\t", $elt->RPKM, "\n";
    }
  }else
  {
    my $int = ${$arrayIntervals}[$ind];
    print $int->identifier,"\t",$int->start,"\t", $int->end, $int->p, "\t", $int->r , "\t", $int->count, "\t", $int->sumSquare, "\t", $int->strand, "\t", $int->RPKM, "\n";
  }
  print "End of Interval List\n", "--------------------------------------------------------------------------\n";
}

#DEBUG
sub printParameters
{
  my ($hashParameters) = @_;
  print "Parameters:\nInterval\tp\tr\n";
  foreach my $key (keys %{$hashParameters})
  {
    print "$key\t", ${$hashParameters}{$key}->p, "\t", ${$hashParameters}{$key}->r, "\n";
  }
  print "----------------------------------------------\n";
}

# Tests the negative binomial in the command line
sub testNBManually
{
  my ($r ,$p) = @_;
  
  print "Testing NB Manually with the parameters p = $p, r = $r\n";
  print "Returns observing x or more reads where x is read from the standart input\n";
  
  while(<>)
  {
    chomp($_);
    print NegativeBinomial($_, $r, $p, 1), "\n";
  }
}

__END__


=head1 NAME

findPeaks.pl

=head1 SYNOPSIS

  To run with the RNASeq data:
  findPeaks.pl -b <bedGraph file> -r <reference bed file> -c <chromosome> -o <output dir> -l <NB parameter file> -rpkm <RPKM file> [-p <p-cut-off value>] [-gap <gap number>]

  To run without the RNASeq data:
  findPeaks.pl -nornaseq -b <bedGraph file> -r <reference bed file> -c <chromosome> -o <output dir> -radii <frame radii> [-p <p-cut-off value>, default 0.01] [-gap <gap number>, default = 2]

  findPeaks.pl -help

  findPeaks.pl -version


For help, run this program with the -help option.

=head1 OPTIONS

=head2 -b <bedGraph file>

The bedGraph file contains the read counts for each position. It contains the RipSeq information.

=head2 -r <reference bed file>

The bed file that contains the genome annotation data.

=head2 -c <chromosome>

Each run of this program is for a specific chromosome.
It is given in this parameter.

=head2 [-p <p-cut-off value>]

The confidence value used in determining the peaks.
This parameter is optional. Its default value is 0.01.

=head2 -o <output dir>

The directory that contains the outputs.
For each chromosome, there is an output file containing peak information.

=head2 -l <NB parameter file>

This must be specified if -nornaseq is not set.
There are two parameters (p and r) of the statistical test used to determine the peaks.
This file contains these parameters for each interval.

=head2 -rpkm <RPKM file>

The RPKM file contains the rpkm values of the intervals.

=head2 -gap <gap number>

 This parameter is optional. Its default value is 2.
 This is the maximum nuber of gaps between the reads that can be tolerated in finding peaks.
 For example, suppose that all readcounts above 9 are high enough to be called peaks.
 If the RipSeq data in the bedgraph file is as

 position -> count
 100 -> 10
 101 -> 10
 102 -> 0
 103 -> 0
 104 -> 10
 105 -> 10
 
 Then we have a peak containing the interval [100,105].
 Note  that the above data wouldn't be counted as one peak if the gap number was 1 or 0.
 
=head2 -version

Display the version number.
  
=head2 -help

Display this documentation.
  
=head1 DESCRIPTION

This program finds the peaks in the given Rip / Clip Seq data.
The Rip / Seq data should be partitioned into chromosomes, i.e., for each  chromosome,
there must be a separate file.
Though you can use this program for calling  peaks without the NB parameter file (with -nornaseq option), 
namely when you don't have RNA Seq data, it is highly recommended to use it with the NB parameter file.
See the user manual for details on the NB arameter file and RNA Seq data.

=head2 Input

=over 5

=item Rip / ClipSeq data: This is provided in the bedGraph file. 
The bedGraph file is given by the -b parameter.

=item RefSeq data       : Genome annotation information. The directory containing this data is given in the -r parameter.

=item Negative Binomial Parameters     : These values are interval specific. They are used in the statistical test used to find the peaks. Given in the -l parameter.

=back

=head2 Output

The output files contain the peaks found in the Rip / Clip Seq data. 
For each  chromosome, a separate output file is created containing the peak information.
The columns of the output file are as follows.


             *Chrom          : Chromosome number
             *strand         : Strand information.
             *Interval       : Interval Identifier.
             *Length         : The length of the interval.  
             *max_height     : Maximum height for a detected peak.
             *footsize       : Footsize of the peak.
             *start          : Start position of the peak.(Inclusive)  
             *end            : Stop position of the peak. (Exclusive)
             *RPKM           : Interval RPKM      
             *weighted_center: Weigthed center for a peak that is a position which is found by summing 
                               up the heights from the beginning of a peak for each position till to pass the 
                               half of the  # all sums of the heights in the peak.
             *pVal           : pVal of a peak. The minimmum pVal is 1e-45.


=head1 EXAMPLE

It is highly recomended to use this  program with RNASeq data though you can use it without RNASeq data.

=head2 With RNASeq Data

If you have the RNASeq data, then you should obtain a lambda file containing the parameter values of the intervals using the RNASeq data.
In order to obtain the lambda file, first find the counts in the exons, then make the RPKM file and then run the lambda calculation
script. After completing these steps, you are ready to run this program as you will have the lambda file. 
The lambda values are used as the parameter of the statistical test used in determining the peaks.
The Rip Seq Data is partitioned into chromosomes. For each chromosome, run this program separately.

Say you want to find the peaks in the first chromosome. 
The RipSeq data, for chromosome 1, is in the file chr1RipSeq.bg in bedGraph format.
Your bed file is /home/user/region.bed and you want to save the found peaks
in the directory /home/user/found_peaks.
Suppose the full path of your lambda file is /home/user/NBParameterFile.txt.
Then  you can run this script as

findPeaks.pl -b chr1RipSeq.bg -r /home/user/region.bed -c chr1 -p 0.01 -o /home/user/found_peaks -l /home/user/NBParameterFile.txt

=head2 Without RNASeq Data

If you don't have the RNASeq data, then you will not have a NB parameter file.
So, simply run this script as above but without the -l parameter and WITH -nornaseq option.

findPeaks.pl -b chr1RipSeq.bg -r /home/user/refseq -c chr1 -p 0.01 -o /home/user/found_peaks -nornaseq

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
 
