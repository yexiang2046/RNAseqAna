#########################################################################################
#                                       func.pm
#########################################################################################
# 
# Common functions for peakfinder.
#
#########################################################################################
# AUTHORS:
#
# Can Cenik, PhD 
# Alper Kucukural, PhD 
# Hakan Ozadam, PhD
#########################################################################################

package lib::func;

use Class::Struct;
use File::Basename;
use Data::Dumper;
use POSIX qw(ceil floor);
 
#################### STRUCTS ######################

struct( 'countstruct', {
 exon_len => '$',
 intron_len=> '$', 
 threeP_len  => '$',
 fiveP_len  => '$',
 exonic_count  => '$', 
 intronic_count  => '$', 
 threeP_count  => '$', 
 fiveP_count  => '$', 
 threeP_ext_count => '$',
 fiveP_ext_count => '$',
});

########################################################################
######################## F U N C T I O N S #############################
########################################################################

###  c o u n t T a g s  ################################################
#
#	This function takes the ref seq data of a gene  and the read counts info
#	from the wiggle file and determines the counts for each exon
#	
#	
#
# INPUT: i) RefSeq Data ($strand, $exonCount, $exon1, $exon2) coming from
#		one line of the RefSeq file
#	ii) the wig file: counts for each position
#	The wig file information (position => counts) is given in hashtags
#
# OUTPUT: Exon lengths ($lengths, reference to an array holding the exon lengths), 
#	  intron lengths ($intronlengths, reference to a scalar, the total of intron lengths) and 
#	counts for each exon and counts 
#	for the introns

sub countTags 
{
 my ($exonlengths, $intronlengths, $lengths, $conlengths, $tagsAll, $hashtags, $strand, $exonCount, $exon1, $exon2) = @_;
 
  # Exon start positions
  # Note that the exon1 (exon start position) is 0 based.

  my @ex1 = split(/\,/, $exon1);

  #Exon end positions
 my @ex2 = split(/\,/, $exon2);
 my $tagcount = 0;

 if($strand eq "+")
 {
   for (my $i = 0; $i < $exonCount; $i++)
   {
     # if (${$conlengths}[$i] =~ /^$/)
     ${$conlengths}[$i] = 0 if (not exists(${$conlengths}[$i]));
     ${$tagsAll}[$i] = 0 if (not exists(${$tagsAll}[$i]) ); 
     
     #NOTE: 
     # if the positions an exon spans the nucleotides  are i,i+1,,,,,,,k,
     # then,In the RefSeq file, in the exon start field (ex1) there is i-1 and in the exon end field (ex2)
     # So the exon starts are 0 based and the exon ends are 1 based.
     # there is k. So when we consider ex1 - ex2 = k - (i -1) = k-i+1, we get the exon length.
     # Hence there is no need to add 1 to the subraction to determine the exon length.
     # A similar argument also works to determine the intron lengths. 

     push(@{$lengths}, $ex2[$i] - $ex1[$i]);
     ${$exonlengths} += ($ex2[$i] - $ex1[$i]);

     if ($i < $exonCount - 1)
     {
      ${$intronlengths} += ($ex1[$i+1] - $ex2[$i]);
     }

     for(my $j = $ex1[$i] + 1; $j <= $ex2[$i]; $j++) 
     {

	if ( exists(${$hashtags}{$j}) )
        {
          $tagcount = ${$hashtags}{$j};
        }
        else        
        {
            $tagcount = 0;
        }
        
        if($tagcount > 0)
        {
	  ${$conlengths}[$i]++;
        }
        
        ${$tagsAll}[$i] += $tagcount;
      }
   }
 }
 else # We are on the - strand
 {
   for (my $i = $exonCount - 1 ; $i >= 0; $i--)
   {
     push(@{$lengths}, $ex2[$i] - $ex1[$i]);

     my $reverseIndex = $exonCount - $i -1;
     ${$conlengths}[$reverseIndex] = 0 if (not exists(${$conlengths}[$reverseIndex])); 
     ${$tagsAll}[$reverseIndex] = 0 if (not exists(${$tagsAll}[$reverseIndex]) );  


     for(my $j = $ex1[$i] + 1; $j <= $ex2[$i]; $j++)
     {
        if ( exists(${$hashtags}{$j}) )
	{
          $tagcount = ${$hashtags}{$j};
        }
        else        
        {
          $tagcount = 0;
        }

        if($tagcount > 0)
        {
         ${$conlengths}[$reverseIndex]++;
        }
        
        ${$tagsAll}[$reverseIndex] += $tagcount;
     }
   }
 }
}

##############################################################
#########  c h e c k P a r a m e t e rs  #############
##############################################################
# Returns 1 if the parameter file is ok
# Return 0 if there is an error
# The error (if any) is reported in the return message
#
# Input:
#        parametersHash: Reference  to the parameters  
#        read from the parameter file 
#
# Output:
#	 returnMessage: Explains the errors if there are any
#	returns 1 if success
#		0 if failure
################################################################

sub checkParameters
{
  my ($parametersHash, $returnMessage) = @_;
  
  my ($success, $failure) = (1,0);
  my $result = $success;
  
  #Parameters needed whether there is rnaseq data or not
  my @essentialParameters = ("MainDir", "ExecDir", "BedDir", "OutDir", "LogDir", "Lib", "Strand", "frameRadii", "RunInCluster", "StepCheckParams", "StepGetChroms", "StepConvert2Bed", "StepCountMappedReads", "StepSortBed", "StepRegionSeparation", "StepChromSeparation", "StepPrepInput", "StepMakeRPKM", "StepNBParameters", "StepPeak", "StepMakeBed", "StepMakeBedGraph","StepSubmitJobs");
  
  #parameters needed when there is rnaseq data
  my @otherParameters = ("RNASeq", "ClusterTimeOut");
  
  if( not exists(${$parametersHash}{"RunInCluster"}))
  {
    ${$returnMessage} .= "Error: Please specifyif it is going to run in cluster(1) or not(0) in the parameter \"RunInCluster\".\n";
      return $failure;
  }
  
  if(${$parametersHash}{"RunInCluster"})
  {
    if(  (not exists(${$parametersHash}{"Queue"})) || (${$parametersHash}{"Queue"} eq "") )
    {
      ${$returnMessage} .= "Error: Please specify the cluster queue in the parameter \"Queue\".\n";
      return $failure;
    }
  }
  
  # First check the essential parameters
  foreach my $param (@essentialParameters)
  {
    if(  (not exists(${$parametersHash}{$param})) || (not ${$parametersHash}{$param})  )
    {
      ${$returnMessage} .= "Error: Please specify the parameter \"$param\"\n";
      $result = $failure;
    }
  }
  
  return $result if ($result == $failure);
  
  # Check the frame radius values
    my @radii = split(",",${$parametersHash}{"frameRadii"});
    
    if(scalar(@radii) < 1)
    {
      ${$returnMessage} .= "Error: Please check the parameter \"frameRadii\"\n";
      ${$returnMessage} .= "       Radii values must be separate by commas\n";
      ${$returnMessage} .= "       Example: frameRadii = 100,500,1000\n\n";
      $result = $failure;
    }else
    {
      foreach my $radius (@radii)
      {
	if( (not ($radius =~ /^\d+\z/)) or ($radius <= 0 ) )
	{
	  ${$returnMessage} .= "Error: The radius value $radius given in the parameter \"frameRadii\" must be a positive integer\n";
	  $result = $failure;
	}
      }
    }
  
  # If there is RNASeq data, we need the remaining parameters as well
  if( exists ${$parametersHash}{"RNASeq"})
  {
    foreach my $otherParam (@otherParameters)
    {
      if(  (not exists(${$parametersHash}{$otherParam})) || (not ${$parametersHash}{$otherParam})  )
      {
	${$returnMessage} .= "Error: Please specify the parameter \"$otherParam\"\n";
	$result = $failure;
      }
    }
    
    return $result if ($result == $failure);
  }
  # If there is RNASeq data, check the remaining parameters
  
  return $result;
}


##############################################################
###############    r u n S t e p s          ##################
##############################################################
sub runSteps
{
    my ($LOG, $params, $jobnum, $type)=@_;
    my $dispatcher=$type;
     lib::io::wlog($LOG, "type=$type"); 
    if ($type ne "")
    {
      $dispatcher->($LOG, $params, $jobnum, $type);
    }
}


sub StepCheckParams
{
   my ($LOG, $params, $jobnum, $type)=@_;
   my $message="";
   if (checkParameters($params, \$message))
   {
     $message="We couldn't detect any problems in the parameter file!";
   }
   lib::io::wlog($LOG, "The parameter check step:\n".$message); 
   return;
}

sub StepConvert2Bed
{
  my ($LOG, $params, $jobnum, $type)=@_;
  my $jobnuminit=${$jobnum};
  my $parallel=${$params}{"RunInCluster"};
  my $outd=${$params}{"OutDir"};
  
  print "StepConvert2Bed is running\n" if(${$params}{Debug});
  
  #for lib
  convert2Bed($LOG, $params, $jobnum, ${$params}{"LibName"}, ${$params}{"Lib"}, ${$params}{"LibFormat"});
  convert2Bed($LOG, $params, $jobnum, ${$params}{"ControlName"}, ${$params}{"Control"}, ${$params}{"ControlFormat"});
  convert2Bed($LOG, $params, $jobnum, ${$params}{"RNASeqName"}, ${$params}{"RNASeq"}, ${$params}{"RNASeqFormat"});
  lib::io::writeJOBWaitLine($parallel,"$outd/JOBLIST.txt", $jobnum, "wait for all convertion to finish", $jobnuminit, ${$jobnum});

}

sub StepGetChroms
{
  my ($LOG, $params, $jobnum, $type) = @_;

  my $beddir=${$params}{"BedDir"};
  my $outd=${$params}{"OutDir"};
  
  my $com="awk '{print \$1}' $beddir/*.bed | sort -u > $outd/chroms.txt";
  `$com`;
}

sub StepRegionSeparation
{
  my ($LOG, $params, $jobnum, $type)=@_;
  my $parallel=${$params}{"RunInCluster"};
  my $jobnuminit=${$jobnum};
  my $outd=${$params}{"OutDir"};
  
  print "StepRegionSeparation is running\n" if(${$params}{Debug});
  
  if (exists ${$params}{"BedDir"})
  {
      separateRegions($LOG, $params, $jobnum);
      lib::io::writeJOBWaitLine($parallel,"$outd/JOBLIST.txt", $jobnum, "wait for all region separation jobs to finish", $jobnuminit, ${$jobnum});
  }
}

sub StepChromSeparation
{
  my ($LOG, $params, $jobnum, $type)=@_;
  my $parallel=${$params}{"RunInCluster"};
  my $jobnuminit=${$jobnum};
  my $outd=${$params}{"OutDir"};
  
  print "StepChromSeparation is running\n" if(${$params}{Debug});
  
  if (exists ${$params}{"BedDir"})
  {
      separateChroms($LOG, $params, $jobnum, ${$params}{"LibName"});
      separateChroms($LOG, $params, $jobnum, ${$params}{"ControlName"});
      separateChroms($LOG, $params, $jobnum, ${$params}{"RNASeqName"});
      lib::io::writeJOBWaitLine($parallel,"$outd/JOBLIST.txt", $jobnum, "wait for all chrom separation jobs to finish", $jobnuminit, ${$jobnum});
  }
}

sub StepPrepInput
{
  my ($LOG, $params, $jobnum, $type)=@_;
  my $parallel=${$params}{"RunInCluster"};
  my $jobnuminit=${$jobnum};
  my $outd=${$params}{"OutDir"};

  print "StepPrepInput is running\n" if(${$params}{Debug});
  
  if (exists ${$params}{"BedDir"})
  {
      prepInput($LOG, $params, $jobnum, ${$params}{"LibName"});
      prepInput($LOG, $params, $jobnum, ${$params}{"ControlName"});
      prepInput($LOG, $params, $jobnum, ${$params}{"RNASeqName"});
      lib::io::writeJOBWaitLine($parallel,"$outd/JOBLIST.txt", $jobnum, "wait for all input preperation jobs to finish", $jobnuminit, ${$jobnum});
  }
}

sub StepCountMappedReads
{
  my ($LOG, $params, $jobnum, $type)=@_;
  my $jobnuminit=${$jobnum};
  my $parallel=${$params}{"RunInCluster"};
  my $outd=${$params}{"OutDir"};

  print "Count mapped reads in all bed files\n" if(${$params}{Debug});
  
  countMappedReads($LOG, $params, $jobnum, ${$params}{"LibName"});
  countMappedReads($LOG, $params, $jobnum, ${$params}{"ControlName"});
  countMappedReads($LOG, $params, $jobnum, ${$params}{"RNASeqName"});
  lib::io::writeJOBWaitLine($parallel,"$outd/JOBLIST.txt", $jobnum, "wait for all mapped counts to finish", $jobnuminit, ${$jobnum});
}

sub StepSortBed
{
  my ($LOG, $params, $jobnum, $type)=@_;
  my $jobnuminit=${$jobnum};
  my $parallel=${$params}{"RunInCluster"};
  my $outd=${$params}{"OutDir"};

  print "Sort all bed files\n" if(${$params}{Debug});
  
  sortBed($LOG, $params, $jobnum, ${$params}{"LibName"}, ${$params}{"Lib"}, ${$params}{"LibFormat"});
  sortBed($LOG, $params, $jobnum, ${$params}{"ControlName"}, ${$params}{"Control"}, ${$params}{"ControlFormat"});
  sortBed($LOG, $params, $jobnum, ${$params}{"RNASeqName"}, ${$params}{"RNASeq"}, ${$params}{"RNASeqFormat"});
  lib::io::writeJOBWaitLine($parallel,"$outd/JOBLIST.txt", $jobnum, "wait sorting to finish", $jobnuminit, ${$jobnum});
}

sub StepMakeRPKM
{
  my ($LOG, $params, $jobnum, $type)=@_;
  my $outd=${$params}{"OutDir"};
  my $parallel=${$params}{"RunInCluster"};
  if (!(${$params}{"nornaseq"}) && (${$params}{"RNASeq"})) {
    print "Make one RPKM file using all RNASeq libraries\n" if(${$params}{Debug});
  
    makeRPKM($LOG, $params, $jobnum, ${$params}{"RNASeqName"});
  }
}

sub StepNBParameters
{
  my ($LOG, $params, $jobnum, $type)=@_;
  my $outd=${$params}{"OutDir"};
  my $parallel=${$params}{"RunInCluster"};
  my $jobnuminit=${$jobnum};

  if (!(${$params}{"nornaseq"}) && (${$params}{"RNASeq"})) {
   print "Estimate NB parameters for each Library\n" if(${$params}{Debug});
  
   estimateNBParameters($LOG, $params, $jobnum, ${$params}{"LibName"});
   estimateNBParameters($LOG, $params, $jobnum, ${$params}{"ControlName"});
   lib::io::writeJOBWaitLine($parallel,"$outd/JOBLIST.txt", $jobnum, "wait for the Negative Binomial(NB) parameter(p, r) estimation to finish", $jobnuminit, ${$jobnum});
  }
}

sub StepPeak
{
  my ($LOG, $params, $jobnum, $type)=@_;
  my $outd=${$params}{"OutDir"};
  my $parallel=${$params}{"RunInCluster"};
  my $jobnuminit=${$jobnum};

  print "Find peaks for each library\n" if(${$params}{Debug});

  findPeaks($LOG, $params, $jobnum, ${$params}{"LibName"});
  findPeaks($LOG, $params, $jobnum, ${$params}{"ControlName"});
  lib::io::writeJOBWaitLine($parallel,"$outd/JOBLIST.txt", $jobnum, "wait for finding all peaks for all the chromosomes", $jobnuminit, ${$jobnum});
}

sub StepCombinePeaks
{
  my ($LOG, $params, $jobnum, $type)=@_;
  my $outd=${$params}{"OutDir"};
  my $parallel=${$params}{"RunInCluster"};
  my $jobnuminit=${$jobnum};

  print "Combine peaks for each library\n" if(${$params}{Debug});

  combinePeaks($LOG, $params, $jobnum, ${$params}{"LibName"});
  combinePeaks($LOG, $params, $jobnum, ${$params}{"ControlName"});
  lib::io::writeJOBWaitLine($parallel,"$outd/JOBLIST.txt", $jobnum, "wait for combining the peaks", $jobnuminit, ${$jobnum});
}


sub StepCalcFDR
{
  my ($LOG, $params, $jobnum, $type)=@_;
  my $outd=${$params}{"OutDir"};
  my $parallel=${$params}{"RunInCluster"};
  my $jobnuminit=${$jobnum};
  if (${$params}{"Control"}) { 
   print "Calculate FDR for each library\n" if(${$params}{Debug});

   calcFDR($LOG, $params, $jobnum, ${$params}{"LibName"}, ${$params}{"ControlName"});
   lib::io::writeJOBWaitLine($parallel,"$outd/JOBLIST.txt", $jobnum, "wait for calculating FDR for the peaks", $jobnuminit, ${$jobnum});
  }
}

#Convert peaks span to bed file to intersect with inputBed file to get real reads in the peak regions.
sub StepPeaks2Bed
{
  my ($LOG, $params, $jobnum, $type)=@_;
  my $peakdir = ${$params}{"PeakOutDir"};
  my $outd=${$params}{"OutDir"};
  my $parallel=${$params}{"RunInCluster"};
  my $jobnuminit=${$jobnum};
  
  peaks2Bed($LOG, $params, $jobnum, ${$params}{"LibName"});
  peaks2Bed($LOG, $params, $jobnum, ${$params}{"ControlName"});
  lib::io::writeJOBWaitLine($parallel,"$outd/JOBLIST.txt", $jobnum, "wait for converting peaks span to Bed file", $jobnuminit, ${$jobnum});  
}

sub StepMakeBedGraph
{
  my ($LOG, $params, $jobnum, $type)=@_;
  my $peakdir = ${$params}{"PeakOutDir"};
  my $outd=${$params}{"OutDir"};
  my $parallel=${$params}{"RunInCluster"};
  my $jobnuminit=${$jobnum};
  
  makeBedGraph($LOG, $params, $jobnum, ${$params}{"LibName"});
  makeBedGraph($LOG, $params, $jobnum, ${$params}{"ControlName"});
  lib::io::writeJOBWaitLine($parallel,"$outd/JOBLIST.txt", $jobnum, "wait for making bedGraph file/s", $jobnuminit, ${$jobnum});  

}

sub StepWigToBigWig
{
  my ($LOG, $params, $jobnum, $type)=@_;
  my $peakdir = ${$params}{"PeakOutDir"};
  my $outd=${$params}{"OutDir"};
  my $parallel=${$params}{"RunInCluster"};
  my $jobnuminit=${$jobnum};
  
  wigToBigWig($LOG, $params, $jobnum, ${$params}{"LibName"});
  wigToBigWig($LOG, $params, $jobnum, ${$params}{"ControlName"});
  lib::io::writeJOBWaitLine($parallel,"$outd/JOBLIST.txt", $jobnum, "wait for wigToBigWig Conversion", $jobnuminit, ${$jobnum});  

}

sub StepSubmitJobs
{
  my ($LOG, $params, $jobnum, $type)=@_;
  my $outd=${$params}{"OutDir"};
  my $libname=${$params}{"LibName"};
  my $parallel=${$params}{"RunInCluster"};
  
  my $nodes="-n ".${$params}{"Nodes"} if (${$params}{"Nodes"});
  my $queue="-q ".${$params}{"Queue"} if (${$params}{"Queue"});
   
  if ($parallel)
  {
   $com= "perl ".${$params}{"ExecDir"}."/submitJobs.pl -j $outd/JOBLIST.txt -s $libname $queue $nodes" ;
   print $com."\n";
   `$com`;
  }
}

sub openChroms
{
   my ($chroms, $beddir, $chromfile)=@_;
  if (!(-s $chromfile)) {
     my $com="awk '{print \$1}' $beddir/*.bed | sort -u > $chromfile";
    `$com`;
  }
  
   open(IN, $chromfile);
   while ($line=<IN>)
   {
     #print $line."\n";
     chomp($line);
     push(@{$chroms}, $line);
   }
}

sub getRegions
{
  my ($reg, $beddir)=@_;
  my $com="ls $beddir/*.bed 2>/dev/null";
  my $res=`$com`;
  my @regions=split(/\n/, $res);
  foreach my $region(@regions)
  {
    push(@{$reg}, basename($region, ".bed"));
  }
}

sub separateChroms
{
  my ($LOG, $params, $jobnum, $name) = @_;

  my $beddir=${$params}{"BedDir"};
  my $parallel=${$params}{"RunInCluster"};
  my $outd=${$params}{"OutDir"};

  my @regions=();
  if (${$params}{"Region"}) {
     push(@regions, ${$params}{"Region"});
  }
  else
  {
   getRegions(\@regions, $beddir);
  }
  
  my @names=split(/:/, $name);

  for(my $i=0; $i<@names; $i++)
  {
   print "Separate Chromosomes for ".$names[$i]."\n" if (${$params}{Debug});
   foreach my $region(@regions)
   {
     my $regoutd=$outd;
     $regoutd=~s/mainfiles/$region/;
     my $chromoutdir="$regoutd/inputBed/".$names[$i];
     lib::io::makeDir($chromoutdir);

     $com="perl ".${$params}{"ExecDir"}."/separateBed.pl $regoutd/inputBed/".$names[$i]."_".$region.".bed $chromoutdir";
     lib::io::writeJOBList($parallel,"$outd/JOBLIST.txt", $jobnum, $com, $outd, "step_separate_chrom.".$region."_".$names[$i].".sh" );
   }
  }
}

sub separateRegions
{
  my ($LOG, $params, $jobnum) = @_;
  print "Separate regions running...\n" if (${$params}{Debug});
  #for lib
  separateForLibs($LOG, $params, $jobnum, ${$params}{"LibName"}, ${$params}{"Lib"}, ${$params}{"LibFormat"});
  separateForLibs($LOG, $params, $jobnum, ${$params}{"ControlName"}, ${$params}{"Control"}, ${$params}{"ControlFormat"});
  separateForLibs($LOG, $params, $jobnum, ${$params}{"RNASeqName"}, ${$params}{"RNASeq"}, ${$params}{"RNASeqFormat"});
}

sub checkBEDTools
{
  my $res=`intersectBed 2>&1 > /dev/null`;
  if ($res=~/command not found/) {
     print "Please install bedtools to run this feature\n";
     print "BEDTools can be downloaded from the link below\n";
     print "http://code.google.com/p/bedtools/\n";
    return 0;
  }
  return 1;
}

sub checkSAMTools
{
  my $res=`samtools 2>&1 > /dev/null`;
  if ($res=~/command not found/) {
     print "Please install samtools to run this feature\n";
     print "SAMTools can be downloaded from the link below\n";
     print "http://samtools.sourceforge.net\n";
    return 0;
  }
  return 1;
}

sub checkWigToBigWig
{
  my $res=`wigToBigWig 2>&1 > /dev/null`;
  if ($res=~/command not found/) {
     print "To use this option, Please install wigToBigWig into your system\n";
     print "Detailed information can be found below.\n";
     print "https://cgwb.nci.nih.gov/goldenPath/help/bigWig.html\n";
     print "You can download the binary file using the link below.\n";
     print "http://hgdownload.cse.ucsc.edu/admin/exe/\n";
    return 0;
  }
  return 1;
}

sub separateForLibs
{
  my ($LOG, $params, $jobnum, $name, $lib, $format)=@_;

  unless (checkBEDTools) { exit;}
  
  my $beddir=${$params}{"BedDir"};
  print "Please define annotaion directory (--beddir your_bed_dir)\n" if (!(-s $beddir));
  my $parallel=${$params}{"RunInCluster"};
  my $outd=${$params}{"OutDir"};
  my $outdir= ${$params}{"OutDir"}."/inputBed";
  my $regionoutd=$outd;
  my $jobnuminit=${$jobnum};
            
  my @names=split(/:/, $name);
  my @libs=split(/:/, $lib);
  my @formats=split(/:/, $format);
    
  my $f=$formats[0];
  my @suffixlist=();
  push(@suffixlist, ".".(uc $f));
  push(@suffixlist, ".".(lc $f));
  
  my @regions=();
  if (${$params}{"Region"}) {
     push(@regions, ${$params}{"Region"});
  }
  else
  {
   getRegions(\@regions, $beddir);
  }
  
  for (my $i=0; $i<@names; $i++)
  {
    if (defined $formats[$i])
    {
      $f=$formats[$i];
      push(@suffixlist, ".".(uc $f));
      push(@suffixlist, ".".(lc $f));
    }

    foreach my $region(@regions)
    {
     print "Separate region: [".$region."] for [".$names[$i]."]\n" if (${$params}{Debug});
	 
     my $regoutd=$regionoutd;
     $regoutd=~s/mainfiles/$region/;
     my $base = basename($libs[$i], @suffixlist);
     my $dir = dirname($libs[$i]);
     #print "\n\n===>$regoutd/inputbed\n\n";
     lib::io::makeDir("$regoutd/inputBed");

     my $com="intersectBed -a $outdir/".$names[$i].".bed -b $beddir/$region.bed -wo -s > $regoutd/inputBed/".$names[$i]."_".$region.".inter;";
     $com.="awk \'{print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$10\"\\t\"\$5\"\\t\"\$6}\' $regoutd/inputBed/".$names[$i]."_".$region.".inter > $regoutd/inputBed/".$names[$i]."_".$region.".bed";
     print $com."\n";
     lib::io::writeJOBList($parallel,"$outd/JOBLIST.txt", $jobnum, $com, $outd, "step_separate_region.".$names[$i].".$region.sh" );
    }
  }
}

sub prepInput
{
  my ($LOG, $params, $jobnum, $name)=@_;
  my $beddir=${$params}{"BedDir"};
  my $parallel=${$params}{"RunInCluster"};
  my $outd=${$params}{"OutDir"};
  my $jobnuminit=${$jobnum};
  my $s="";
  $s=" -s" if (${$params}{"Strand"});
  
  my @chroms=();
  openChroms(\@chroms, $beddir, "$outd/chroms.txt");
  my @names=split(/:/, $name);
  
  my @regions=();
  if (${$params}{"Region"}) {
     push(@regions, ${$params}{"Region"});
  }
  else
  {
   getRegions(\@regions, $beddir);
  }
  
  for (my $i=0; $i<@names; $i++)
  {
    print "prepInput: ".$names[$i]."\n" if (${$params}{Debug});
    foreach my $region(@regions)
    {
      my $regoutd=$outd;
      my $inputsetdir="$regoutd/inputBed/".$names[$i];   
      $regoutd=~s/mainfiles/$region/;
      my $outsetdir="$regoutd/inputBG/".$names[$i];
      lib::io::makeDir("$outsetdir");

      $inputsetdir=~s/mainfiles/$region/;
      foreach my $chrom(@chroms)
      {
	my $input="$inputsetdir/$chrom.bed";
	my $output="$outsetdir";

        $com="perl ".${$params}{"ExecDir"}."/prepInput.pl -i $input -o $output $s -r $beddir/$region.bed -c $outd/inputBed/".$names[$i];
	#print $com."\n";
        lib::io::writeJOBList($parallel,"$outd/JOBLIST.txt", $jobnum, $com, $outd, "step_prep_".$names[$i]."_".$chrom.".$region.sh" );
      }
    }
  } 
}

sub makeRPKM
{
  my ($LOG, $params, $jobnum, $name)=@_;
  my $beddir=${$params}{"BedDir"};
  my $parallel=${$params}{"RunInCluster"};
  my $outd=${$params}{"OutDir"};
  my @watcriarray=();
  if (${$params}{"Strand"}) {
    push(@watcriarray, "watson");
    push(@watcriarray, "crick");
  }
  else
  { 
    push(@watcriarray, "all");
  }
 
  my @chroms=();
  print "$outd/chroms.txt\n";
  openChroms(\@chroms, $beddir, "$outd/chroms.txt");
  my @names=split(/:/, $name);

  my @regions=();
  
  if (${$params}{"Region"}) {
     push(@regions, ${$params}{"Region"});
  }
  else
  {
   getRegions(\@regions, $beddir);
  }
  foreach my $region(@regions)
  {
    my $regoutd=$outd;
    $regoutd=~s/mainfiles/$region/;
    
    foreach my $watcri (@watcriarray)
    {
      my $jobnuminit=${$jobnum};
      my $outsetdir="$regoutd/rpkm/$watcri";
      lib::io::makeDir("$outsetdir");
      my $preinputcounts="";
      for (my $i=0; $i<@names; $i++)
      {	
        $preinputcounts.="$regoutd/inputBG/".$names[$i]."/".$watcri."counts/\[chrom\].intervalcount:";
      }
      chop($preinputcounts);
      foreach my $chrom (@chroms)
      {
	
        #print "$chrom\n";
	my $com="";
	my $inputcounts=$preinputcounts;
	$inputcounts=~s/\[chrom\]/$chrom/gi;
        if (@names>1)
        {
          $com="perl ".${$params}{"ExecDir"}."/makeRPKM.pl -i $inputcounts -o $regoutd/rpkm/$watcri/$chrom.csv";
	  #print $com."\n";
        }
        else
        {
         if ((-s $inputcounts) || ${$params}{"RunInCluster"})
         { 
	   $com="awk \'{print \$1\"\\t\"\$2\"\\t\"\$6}\' $inputcounts > $regoutd/rpkm/$watcri/$chrom.csv";
         }
         else
         {
           print "Warning: RPKM interval count file doesn't exist $inputcounts\n";
         }
        }
	#print $com."\n";
	lib::io::writeJOBList($parallel,"$outd/JOBLIST.txt", $jobnum, $com, $outd, "step_make_RPKM_$watcri.$chrom.$region.sh" );
      }
      lib::io::writeJOBWaitLine($parallel,"$outd/JOBLIST.txt", $jobnum, "wait for the peak detection to finish", $jobnuminit, ${$jobnum});
       $com="cat $regoutd/rpkm/$watcri/*.csv>$regoutd/rpkm/RPKM$watcri.csv 2>/dev/null";
       $jobnuminit=${$jobnum};
       lib::io::writeJOBList($parallel,"$outd/JOBLIST.txt", $jobnum, $com, $outd, "step_merge_RPKM_$watcri.$region.sh" );
       lib::io::writeJOBWaitLine($parallel,"$outd/JOBLIST.txt", $jobnum, "wait for RPKM merge to finish", $jobnuminit, ${$jobnum});
    }
  }
}

sub estimateNBParameters
{
  my ($LOG, $params, $jobnum, $name)=@_;
  my $beddir=${$params}{"BedDir"};
  my $parallel=${$params}{"RunInCluster"};
  my $outd=${$params}{"OutDir"};
  my @watcriarray=();
  my $estimationRadius = ${$params}{"estimationRadius"};

  if (${$params}{"Strand"}) {
    push(@watcriarray, "watson");
    push(@watcriarray, "crick");
  }
  else
  { 
    push(@watcriarray, "all");
  }
  my @names=split(/:/, $name);
   
  my @regions=();
  if (${$params}{"Region"}) {
     push(@regions, ${$params}{"Region"});
  }
  else
  {
   getRegions(\@regions, $beddir);
  }
  
  foreach my $region(@regions)
  {
    my $regoutd=$outd;
    $regoutd=~s/mainfiles/$region/;
    foreach my $watcri (@watcriarray)
    {
      for (my $i=0; $i<@names; $i++)
      {
        my $rpkmfile.="$regoutd/rpkm/RPKM$watcri.csv";
	my $inputBGcounts = "$regoutd/inputBG/".$names[$i]."/".$watcri."counts";
	my $outsetdir="$regoutd/inputBG/".$names[$i]."/".$watcri."NBparams";
        lib::io::makeDir("$outsetdir");
        my  $outputfile="$outsetdir/NB$watcri.csv";
        $com="perl ".${$params}{"ExecDir"}."/estimateNBParameters.pl -rpkm $rpkmfile -i $inputBGcounts -o $outputfile -radius $estimationRadius";
        print "$com\n";
	lib::io::writeJOBList($parallel,"$outd/JOBLIST.txt", $jobnum, $com, $outd, "step_estimateNBParams_$watcri.$names[$i].$region.sh" );
      }
    }
  }
}

sub findPeaks
{
  my ($LOG, $params, $jobnum, $name)=@_;
  my $beddir=${$params}{"BedDir"};
  my $parallel=${$params}{"RunInCluster"};
  my $outd=${$params}{"OutDir"};
  my $jobnuminit=${$jobnum};
  my $gapNumber = (exists ${$params}{"gapNumber"} and ${$params}{"gapNumber"} >= 0 ) ? ${$params}{"gapNumber"} : 2;

  if (${$params}{"Strand"}) {
    push(@watcriarray, "watson");
    push(@watcriarray, "crick");
  }
  else
  { 
    push(@watcriarray, "all");
  }
  
  my @chroms=();
  openChroms(\@chroms, $beddir, "$outd/chroms.txt");
  my @names=split(/:/, $name);
  
  my @regions=();
  if (${$params}{"Region"}) {
     push(@regions, ${$params}{"Region"});
  }
  else
  {
   getRegions(\@regions, $beddir);
  }
  
  for (my $i=0; $i<@names; $i++)
  {
    print "findPeaks: ".$names[$i]."\n" if (${$params}{Debug});
    foreach my $region(@regions)
    {
      my $regoutd=$outd;
      $regoutd=~s/mainfiles/$region/;
      my $inputsetdir="$regoutd/inputBG/".$names[$i]; 
      foreach my $watcri (@watcriarray)
      {
	my $rpkmfile.="$regoutd/rpkm/RPKM$watcri.csv";
	my $NBparamfile.="$regoutd/inputBG/".$names[$i]."/".$watcri."NBparams/NB$watcri.csv";
	my $output="$regoutd/foundPeaks/$names[$i]/$watcri";
	lib::io::makeDir("$output");
        foreach my $chrom(@chroms)
        {
	  print "Find peaks for $region $chrom $watcri \n" if(${$params}{Debug});
	  my $input="$inputsetdir/$watcri/$chrom.bg";
	  if (!(${$params}{"nornaseq"}) && (${$params}{"RNASeq"})) {
            $com="perl ".${$params}{"ExecDir"}."/findPeaks.pl -b $input -o $output -rpkm $rpkmfile -l $NBparamfile -r $beddir/$region.bed -c $chrom -gap $gapNumber -p 0.01";
	  }
	  else
	  {
	    $com="perl ".${$params}{"ExecDir"}."/findPeaks.pl -b $input -o $output -nornaseq -r $beddir/$region.bed -c $chrom -gap $gapNumber -p 0.01";
	  }
  	  #print $com."\n";
          lib::io::writeJOBList($parallel,"$outd/JOBLIST.txt", $jobnum, $com, $outd, "step_findpeaks_".$names[$i]."_".$watcri."_".$chrom.".$region.sh" );
	}
      }
    }
  } 
}

sub combinePeaks
{
  my ($LOG, $params, $jobnum, $name)=@_;
  my $beddir=${$params}{"BedDir"};
  my $parallel=${$params}{"RunInCluster"};
  my $outd=${$params}{"OutDir"};
  my $jobnuminit=${$jobnum};

  if (${$params}{"Strand"}) {
    push(@watcriarray, "watson");
    push(@watcriarray, "crick");
  }
  else
  { 
    push(@watcriarray, "all");
  }
  
  my @names=split(/:/, $name);
  
  my @regions=();
  if (${$params}{"Region"}) {
     push(@regions, ${$params}{"Region"});
  }
  else
  {
   getRegions(\@regions, $beddir);
  }
  
  for (my $i=0; $i<@names; $i++)
  {
    print "combinePeaks: ".$names[$i]."\n" if (${$params}{Debug});
    my $foundpeaksdir=$outd;
    $foundpeaksdir=~s/mainfiles/foundpeaks/g;
    $foundpeaksdir.="/".$names[$i]."/peaks";
    #print $foundpeaksdir."\n";
    lib::io::makeDir("$foundpeaksdir");
    foreach my $region(@regions)
    {
      my $regoutd=$outd;
      $regoutd=~s/mainfiles/$region/;

      my $peakfiles="$regoutd/foundPeaks/$names[$i]/*";
      $com="grep -h -v '^#' $peakfiles/*.peaks > $foundpeaksdir/$region.reg.peaks";
      lib::io::writeJOBList($parallel,"$outd/JOBLIST.txt", $jobnum, $com, $outd, "step_combinepeaks_".$names[$i].".$region.sh" );

    }
    #$jobnuminit=${$jobnum};
    
    lib::io::writeJOBWaitLine($parallel,"$outd/JOBLIST.txt", $jobnum, "wait for combining peaks to finish for each region", $jobnuminit, ${$jobnum});
    $com="grep -h -v '^#' $foundpeaksdir/*.reg.peaks > $foundpeaksdir/".$names[$i].".peaks";
    lib::io::writeJOBList($parallel,"$outd/JOBLIST.txt", $jobnum, $com, $outd, "step_combinepeaks_".$names[$i].".sh" );
  } 
}

sub calcFDR
{
  my ($LOG, $params, $jobnum, $libname, $controlname)=@_;
  my $beddir=${$params}{"BedDir"};
  my $parallel=${$params}{"RunInCluster"};
  my $outd=${$params}{"OutDir"};
  my $jobnuminit=${$jobnum};
  
  if (${$params}{"Strand"}) {
    push(@watcriarray, "watson");
    push(@watcriarray, "crick");
  }
  else
  { 
    push(@watcriarray, "all");
  }
  
  my @libnames=split(/:/, $libname);
  my @controlnames=split(/:/, $controlname);
  
  my @regions=();
  if (${$params}{"Region"}) {
     push(@regions, ${$params}{"Region"});
  }
  else
  {
   getRegions(\@regions, $beddir);
  }
  
  for (my $i=0; $i<@libnames; $i++)
  {
    print "Calculate FDR: ".$libnames[$i]."\n" if (${$params}{Debug});
    my $peaksdir=$outd;
    $peaksdir=~s/mainfiles/foundpeaks/g;
    my $foundlibpeaksdir="$peaksdir/".$libnames[$i]."/peaks";

    foreach my $region(@regions)
    {
      my $regoutd=$outd;
      $regoutd=~s/mainfiles/$region/;
      my $controllibname="";
      my $com="";
      if (@controlnames>1) {
	my $cattxt="";
	for(my $j=0; $j<@controlnames; $j++)
	{
	  $cattxt.=" $peaksdir/".$controlnames[$j]."/peaks/$region.reg.peaks ";
	}
	$controllibname="$peaksdir/".$libnames[$i]."/peaks/control.$region.reg.peaks";
	$com="grep -h -v '^#' $cattxt>$controllibname;";
      }
      else
      {
	$controllibname="$peaksdir/".$controlnames[0]."/peaks/$region.reg.peaks";
      }
      
      $com.="perl ".${$params}{"ExecDir"}."/calcFDR.pl -l $foundlibpeaksdir/$region.reg.peaks -c $controllibname";
      #print $com."\n";
      lib::io::writeJOBList($parallel,"$outd/JOBLIST.txt", $jobnum, $com, $outd, "step_calcFDR_".$libnames[$i].".$region.sh" );

    }

    if (@controlnames>1) {
      my $cattxt="";
      for(my $j=0; $j<@controlnames; $j++)
      {
	  $cattxt.=" $peaksdir/".$controlnames[$j]."/peaks/$controlnames[$i].peaks ";
      }
      $controllibname="$peaksdir/".$libnames[$i]."/peaks/control.peaks";
      $com="grep -h -v '^#' $cattxt>$controllibname;";
    }
    else
    {
	$controllibname="$peaksdir/$controlnames[0]/peaks/$controlnames[0].peaks";
    }
    $com.="perl ".${$params}{"ExecDir"}."/calcFDR.pl -l $foundlibpeaksdir/$libnames[$i].peaks -c $controllibname";
    print $com."\n";
    lib::io::writeJOBList($parallel,"$outd/JOBLIST.txt", $jobnum, $com, $outd, "step_calcFDR_".$libnames[$i].".sh" );
  } 
}

sub peaks2Bed
{
  my ($LOG, $params, $jobnum, $name)=@_;
  my $parallel=${$params}{"RunInCluster"};
  my $outd=${$params}{"OutDir"};
  my @names=split(/:/, $name);
  my $peaksdir=$outd;
  my $inputbeddir=$outd."/inputBed";
  $peaksdir=~s/mainfiles/foundpeaks/g;

  for (my $i=0; $i<@names; $i++)
  {
    my $foundpeaksdir.="$peaksdir/$names[$i]/peaks";
    my $foundpeaksinterdir = "$peaksdir/$names[$i]/inter";
    my $foundpeaksbeddir = "$peaksdir/$names[$i]/bed";
    my $foundpeaksspandir = "$peaksdir/$names[$i]/span";
    lib::io::makeDir("$foundpeaksinterdir");
    lib::io::makeDir("$foundpeaksbeddir");
    lib::io::makeDir("$foundpeaksspandir");
    
    my $com = "awk \'{if(\$0!~/^#/){print \$1\"\\t\"\$8\"\\t\"\$9\"\\t\"\$3\"\\t0\\t\"\$2}}\' $foundpeaksdir/".$names[$i].".peaks > $foundpeaksspandir/".$names[$i].".span.bed;\n";
    $com.="sort -k1,1b -k2,2n  $foundpeaksspandir/".$names[$i].".span.bed >  $foundpeaksspandir/".$names[$i].".span.sorted.bed;\n";
    $com.="intersectBed -a $inputbeddir/".$names[$i].".bed -b $foundpeaksspandir/".$names[$i].".span.sorted.bed -wo -s > $foundpeaksinterdir/".$names[$i].".inter;\n";
    $com.= "awk \'{print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$10\"\\t0\\t\"\$6}\' $foundpeaksinterdir/".$names[$i].".inter > $foundpeaksbeddir/".$names[$i].".pre.bed;\n";
    $com.="sort -k1,1b -k2,2n $foundpeaksbeddir/".$names[$i].".pre.bed >  $foundpeaksbeddir/".$names[$i].".bed\n";
 
    lib::io::writeJOBList($parallel,"$outd/JOBLIST.txt", $jobnum, $com, $outd, "step_peaks2Bed_".$names[$i].".sh" );
  }
}
sub makeBedGraph
{
  my ($LOG, $params, $jobnum, $name)=@_;
  my $parallel=${$params}{"RunInCluster"};
  my $outd=${$params}{"OutDir"};
  my @names=split(/:/, $name);
  my $peaksdir=$outd;
  $peaksdir=~s/mainfiles/foundpeaks/g;

  for (my $i=0; $i<@names; $i++)
  {
    my $foundpeaksbedfile = "$peaksdir/$names[$i]/bed/".$names[$i].".bed";
    my $foundpeaksoutdir = "$peaksdir/$names[$i]/bedGraph";
    lib::io::makeDir("$foundpeaksoutdir");
    my $com="perl ".${$params}{"ExecDir"}."/makeBedGraph.pl -i $foundpeaksbedfile -o $foundpeaksoutdir -n ".$names[$i];
    #print $com;
    lib::io::writeJOBList($parallel,"$outd/JOBLIST.txt", $jobnum, $com, $outd, "step_makeBedGraph_".$names[$i].".sh" );
  }
}


sub wigToBigWig
{
  my ($LOG, $params, $jobnum, $name)=@_;
  my $parallel=${$params}{"RunInCluster"};
  my $outd=${$params}{"OutDir"};
  #print "OUTDIR=".${$params}{"OutDir"}."\n";
  unless (checkWigToBigWig) { exit;}
  # fetch chrom sizes if this script doesn't work please create a tab separated file
  # put chrom and size of the  chrom into the columns.
  # Ex: chr1  247249719
  #     chr2  242951149  
  #     chr3  199501827  
  
  #print  "$outd/chrom.sizes\n";
  if (!(-s "$outd/chrom.sizes")) {
    my $db=${$params}{"db"};
    my $bashCommand = ${$params}{"MainDir"} . "/scripts/fetchChromSizes.sh " . "$db > $outd/chrom.sizes";
    `$bashCommand`; 
  }
  my @names=split(/:/, $name);
  my $peaksdir=$outd;
  
  $peaksdir=~s/mainfiles/foundpeaks/g;
  my @arr=("", "watson", "crick");
  #my @arr=("watson");
  for (my $i=0; $i<@names; $i++)
  {
    for(my $j=0; $j<@arr; $j++)
    {
      my $watcri="";
      if ($arr[$j] !~/^$/)
      {
	$watcri=".".$arr[$j];
      }
      
      my $bedGraphFile = "$peaksdir/$names[$i]/bedGraph/".$names[$i]."$watcri.bg";
      
      my $bigWigFile = "$peaksdir/$names[$i]/bedGraph/".$names[$i]."$watcri.bw";
      my $com="wigToBigWig -clip -itemsPerSlot=1 $bedGraphFile  $outd/chrom.sizes $bigWigFile";
      print $com."\n";
      lib::io::writeJOBList($parallel,"$outd/JOBLIST.txt", $jobnum, $com, $outd, "step_wigToBigWig_".$names[$i]."$watcri.sh" );
      
    }
  }
}
sub convert2Bed
{
  my ($LOG, $params, $jobnum, $name, $lib, $format)=@_;
  my $parallel=${$params}{"RunInCluster"};
  my $outd=${$params}{"OutDir"};
  my $jobnuminit=${$jobnum};
  my $outdir= ${$params}{"OutDir"}."/inputBed";

  lib::io::makeDir($outdir);

  my @names=split(/:/, $name);
  my @libs=split(/:/, $lib);
  my @formats=split(/:/, $format);
  
  my $f=$formats[0];
  for (my $i=0; $i<@libs; $i++)
  {
    if (!(-s $libs[$i])) {
       print $libs[$i]." doesn't exist\nPlease check this file location and run again\n";
    }
    if (${$params}{Debug}) {
       print "Convert2Bed: ".$libs[$i]."\n";
    }
    if (defined $formats[$i])
    {
      $f=$formats[$i];
    }
    else
    {
      $libs[$i]=~/\.([^.]+)$/;
      $f=$1;
    }
    my $base = basename($libs[$i], $f);
    my $dir = dirname($libs[$i]);
    my $com="";
    if ($f=~/SAM/i)
    {
      $com=sam2bed($params, $libs[$i], "$outdir/".$names[$i].".bed");
    }
    if ($f=~/BED/i)
    {
      
      if (!(-s "$outdir/".$names[$i].".bed"))
      {
	$com="ln -s $libs[$i] $outdir/".$names[$i].".bed";
        #print $com."\n"
      }  
    }
    if ($f=~/BOW/i) {
      $com=bowtie2bed($libs[$i], "$outdir/".$names[$i].".bed");
    }
    if ($f=~/BAM/i) {
      $com=bam2bed($params, $libs[$i], "$outdir/".$names[$i].".bed");
    }
    
    #print "COM:".$com."\n";
    if ($com !~/^$/) {
      lib::io::writeJOBList($parallel,"$outd/JOBLIST.txt", $jobnum, $com, $outd, "step_convert2bed_for_".$names[$i]."_from_".$f.".sh" );
    }    
  }
}

sub countMappedReads
{
  my ($LOG, $params, $jobnum, $name)=@_;
  my $parallel=${$params}{"RunInCluster"};
  my $outd=${$params}{"OutDir"};
  my $jobnuminit=${$jobnum};
  my $outdir= ${$params}{"OutDir"}."/inputBed";
  my @names=split(/:/, $name);
  for (my $i=0; $i<@names; $i++)
  {
    if (${$params}{"Strand"})
    {
      $com="grep \"[[:space:]]+\" $outdir/".$names[$i].".bed|wc -l > $outdir/".$names[$i]."_watson_mapped.count";
      lib::io::writeJOBList($parallel,"$outd/JOBLIST.txt", $jobnum, $com, $outd, "step_count_watson_".$names[$i].".sh" );
      $com="grep \"[[:space:]]-\" $outdir/".$names[$i].".bed|wc -l > $outdir/".$names[$i]."_crick_mapped.count";
      lib::io::writeJOBList($parallel,"$outd/JOBLIST.txt", $jobnum, $com, $outd, "step_count_crick_".$names[$i].".sh" );
    }
    else
    {
      $com="wc -l $outdir/".$names[$i].".bed > $outdir/".$names[$i]."_all_mapped.count";
      lib::io::writeJOBList($parallel,"$outd/JOBLIST.txt", $jobnum, $com, $outd, "step_count_all_".$names[$i].".sh" );
    }
  }
}

sub sortBed
{
  my ($LOG, $params, $jobnum, $name, $lib, $format)=@_;
  my $parallel=${$params}{"RunInCluster"};
  my $outd=${$params}{"OutDir"};
  my $jobnuminit=${$jobnum};
  my $outdir= ${$params}{"OutDir"}."/inputBed";
    
  my @names=split(/:/, $name);
  my @libs=split(/:/, $lib);
  my @formats=split(/:/, $format);
    
  my $f=$formats[0];
  for (my $i=0; $i<@names; $i++)
  {
    if (${$params}{Debug}) {
      print "Bed files will be sorted: ".$names[$i]."\n";
    }
    if (defined $formats[$i])
    {
      $f=$formats[$i];    
    }
    #print $libs[$i]."\n";
    my $base = basename($libs[$i], $f);
    my $dir = dirname($libs[$i]);
    my $com="mv $outdir/".$names[$i].".bed $outdir/".$names[$i].".presort.bed;sort -k1,1b -k2,2n $outdir/".$names[$i].".presort.bed > $outdir/".$names[$i].".bed";
    lib::io::writeJOBList($parallel,"$outd/JOBLIST.txt", $jobnum, $com, $outd, "step_sortBed_for_".$names[$i].".sh" );
  }
}
    
sub bowtie2bed
{
 my ($input, $output) =@_;
 
 my $com="awk \'{offset=0; if(\$4~/^[+-]\$/){offset=2}; if(\$3~/^[+-]\$/){offset=1} ;print \$(3+offset)\"\\t\"\$(4+offset)\"\\t\"\$(4+offset)+length(\$(5+offset))\"\\t\"\$1\"\\t\"\$(7+offset)\"\\t\"\$(2+offset)  }\' $input > $output";
 #print $com."\n";
 return $com;
}

sub bam2bed
{
 my ($params, $input, $output) =@_;
 
 unless (checkBEDTools) { exit;}
   
 my $com= "bedtools bamtobed -split -i $input > $output";  	
 
 return $com;
}

sub sam2bed
{
 my ($params, $input, $output) =@_;
 my $dir = dirname($output);
 my $base = basename($input, (".SAM", ".sam"));
 my $inputbam="$dir/$base.bam";
 unless (checkSAMTools) { exit;}

 #Please change this line to convert sam files to bam files according to your sam file
 my $com="samtools view -bS $input > $inputbam;";
 $com.=bam2bed($params, $inputbam, $output);

 return $com;
}

return 1;
