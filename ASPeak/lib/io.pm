#########################################################################################
#                                       io.pm
#########################################################################################
# 
# Input and output functions for peakfinder.
#
#########################################################################################
# AUTHORS:
#
# Can Cenik, PhD 
# Alper Kucukural, PhD 
# Hakan Ozadam, PhD
#########################################################################################


############## LIBRARIES AND PRAGMAS ################
package lib::io;

use Class::Struct;
use Data::Dumper;
use File::Basename;
use Cwd 'abs_path';


use POSIX qw(ceil floor);

#################### STRUCTS ######################
struct( 'Region', {
   chrom         => '$',
   start         => '$', 
   end           => '$',
   name          => '$',
   strand        => '$',
} );


struct( 'peakParams', {
   WigDir           => '$',
   OutDir           => '$',
   OutWigDir        => '$',
   CountWigDir      => '$',   
   CountOutDir      => '$',
   LambdaCountDir   => '$',
   LambdaOutDir     => '$',
   PeakInputWigDir  => '$',
   PeakLambdaFile   => '$',
   PeakOutDir       => '$',
   PeakWigDir       => '$',
   RefDir           => '$',
   RPKM             => '$',
   GeneExon         => '$',
   Wig              => '$',
   Html             => '$',  
   Count            => '$',
} );


struct( 'BEDFILE', {  
   start            => '$',
   end              => '$',
} );
#################### CONSTANTS ###################### 

#################### FUNCTIONS ######################
sub makeDir
{
  my ($path) = @_;
  if (!(-s $path)) {
    `mkdir -p $path`;
  }
  unless (-d $path)
  {
    die "Cannot create directory '$path'";
  }
}

sub readPipeParams
{
  my ($params, $paramfile, $commandParams)=@_; 
  my $bar = eval(Dumper($commandParams));
  
  #print Dumper($commandParams), Dumper($bar); # pretty print (no array indices)

  
  open(PAR, $paramfile) || die "ioError opening: $paramfile\n";

  my %par=();
  my @order=();
  
  setMainDir(\%par, \@order);
    
  while (my $line=<PAR>)
  {
    chomp($line);
    if ($line!~/^#/)
    {
        if($line=~/^([A-Za-z0-9]+)[\s\t]+=[\s\t]+(.*)/)
        { 
          my $param=$1;
          my $val=$2;
	  my $res=eval("\$commandParams->".$param);
	  
	  if (exists $par{$param})
          {          
            $val=$par{$param}.":$val";       
          }
          else
  	  {
            push(@order, $param); 
          }
	  if ($res!~/^$/)
	  {
	    $val=$res;
	  }
	  if ($param eq "OutDir")
	  {
	    $val.="/[region]";
	  }
	  else
	  {
            if ($val=~/\$([^\/]*)\//)
            {
              my $reppar=$1;
              my $rep=$par{$reppar};
              $val=~s/\$$reppar/$rep/gi;
            }
	  }
	  if ($commandParams->Region =~/^$/)
	  { 
	     $val=~s/\[region\]/mainfiles/gi;
	  }
	  else
	  {
	    my $region=$commandParams->Region;
	    $val=~s/\[region\]/$region/gi;
	  }
          $par{$param}=$val;     
        } 
    }
  }
  
  addExtraCommandParams(\%par,\@order, $commandParams);
  
  $par{"order"}=\@order;
  %{$params}=%par;
  close(PAR);
}

sub addExtraCommandParams
{
  my ($par, $order, $commandParams) = @_;
  
  my @commands = ("LibName", "Lib", "LibFormat", "ControlName", "Control", "ControlFormat", "RNASeqName", "RNASeq" , "RNASeqFormat", "BedDir", "Nodes", "Strand", "Region",  "OutDir", "nornaseq", "RunInCluster");

  foreach my $param (@commands)
  {
    if (!exists ${$par}{$param})
    {
      my $res=eval("\$commandParams->".$param);
      if ($res!~/^$/)
      {
	print "$param -> $res\n";
	if ($param=~/OutDir/) {
	  $res.="/mainfiles";
	}
	${$par}{$param}=$res;
	push(@{$order}, $param);
      }   
    }
  }
  
  setFormatAndName("Lib", $par, $order);
  setFormatAndName("Control", $par, $order);
  setFormatAndName("RNASeq", $par, $order);
}

sub setMainDir
{
  my($par, $order) = @_;
  
  my($mainScriptName, $scriptDirectory, $suffix) = fileparse(abs_path($0)); # Get the directory of the working script
  my @scriptPathPieces = split(/\//,$scriptDirectory); # Make an array of subdirectories
  pop(@scriptPathPieces); # Since this script is under "/dir1/dir2/.../dirN/scripts" directory, it pops "scripts" from the array
  my $baseDir = join("/", @scriptPathPieces); #This is the main directory containing ASPeak files. Hence  it is our main directory.
  
  push(@{$order}, "MainDir");
  ${$par}{"MainDir"} = $baseDir;
}


#If format and name is empty, this will set using name of the file to the name, extension of the files to the extension
sub setFormatAndName
{
  my ($libname, $par, $order) = @_;
  my $checkExist=!exists ${$par}{$libname."Name"};
  if (exists ${$par}{$libname}) {
    my $lname=${$par}{$libname};
    my @allarr=split(/:/, $lname);
    my $c=0;
    foreach my $lib (@allarr)
    {
      my ($base, $path, $ext)=fileparse($lib, qr/\.[^.]*/);
      $ext=~s/\.//;
      my $col="";
      if ($c<@allarr-1) {
	$col=":";
      }
      if ($checkExist)
      {
	if ($c==0) {push(@{$order}, $libname."Name");}
	
        ${$par}{$libname."Name"}.=$base.$col;
      }
      if (!exists ${$par}{$libname."Format"})
      {
	if ($c==0) {push(@{$order}, $libname."Format");}
	
        ${$par}{$libname."Format"}.=$ext.$col;
      }
      $c++;
    }
  }
}


sub writeParams
{
    my ($LOG, $params) = @_;
   
    print $LOG "\n############ PARAMETERS #############################\n#\n";
    my $order=${$params}{"order"};
    foreach my $param (@{$order})       
    {
      if($param!~/^$/)
      {
       
       my $txt = sprintf("# %-20s = %-40s\n", $param, ${$params}{$param}); 
       my @steps=split(/:/, ${$params}{$param});
       if (@steps>1)
       { 
          $txt=sprintf("# %-20s ===>\n", $param);
          for  (my $i=0; $i<@steps; $i++)
          { 
            $txt .= sprintf("# %20s = %-40s\n", ($i+1), $steps[$i]); 
          }
       }
       print $LOG $txt;
      }
    } 
    
    print $LOG "#\n############ PARAMETERS END ##########################\n";   
}

sub writeExp
{
  my ($cmd, $ARGV, $paramfile, $params) =@_;
  
  my @mainfilesDirPieces = split(/\// , abs_path(${$params}{"OutDir"}) );
  pop(@mainfilesDirPieces); # the innermost directory is mainfiles. So we get rid of it here.
  my $outDir = join("/", @mainfilesDirPieces);
  
  if (${$params}{"Explanations"})
  {
  foreach my $c (@{$ARGV})
  { 
    $cmd.= " ".$c;
  }
  if (${$params}{"RunInCluster"})
  {
     print "########### Cluster Runs ###########\n";
     print "  1. ASPeak will be run in cluster. To check JOBLIST please consult JOBLIST.txt file.\n\n";
     print "  ".$outDir."/mainfiles/JOBLIST.txt\n";
     print "  2. STDOUT of the runs will be hold in the files under ".${$params}{"MainDir"}."/sge\n";
     print "  3. The scripts will be run can be found ".$outDir."/mainfiles/scripts\n";
     if (${$params}{"StepSubmitJobs"})
     {
         print "   4. Jobs will be submitted automatically.\n ";  
         print "      If you want to submit the jobs manually set StepSubmitJobs parameter to 0 and run the command below.\n\n";
         print "   $cmd --type StepSubmitJobs\n";
     }
     else
     {
         print "  4. These jobs will not be submitted. To submit the jobs run the command below.\n\n";
         print "   $cmd --type StepSubmitJobs\n\n";
         print "  5. If you want automatic submition please set StepSubmitJobs parameter to 1\n and run the command below\n\n";
         print "   $cmd \n\n";
     }
  }
  print "\n########### Running Any Single Step ###########\n";
  print "  1.To run any single step. Please, run the command below with the name of the step you want to run\n\n";
  print "   $cmd --type [StepName]\n";
  print "  2. Available steps:\n";
  print `grep Step $paramfile`;
 
  if (${$params}{"StepPeak"})
  { 
   print "\n########### Peak Outputs ###########\n";
   print "  1. Called peaks will be written to the directory below for each chromosome.\n";  
   print "  Peak Out Dir: ". $outDir ."/foundpeaks/" .${$params}{"LibName"}. "/peaks\n";
   if (${$params}{"StepMakeBedGraph"})
   { 
     print "  2. BedGraphs of the called peaks will be available in the directory below. If you set StepWigToBigWig, the output also will be in the directory below\n";  
     print "  Peak bedGraph Dir: ". $outDir ."/foundpeaks/" . ${$params}{"LibName"} . "/bedGraph\n";
   }
  }

  print "\n########### LOG Directory ###########\n";
  print "  1. LOG files can be found under the directory below.\n";  
  print "  Log Dir: ". $outDir ."/mainfiles/LOGS\n";
  }
}

sub checkLibDir
{
  my $params=$_[0];
  # If bedgraph files haven't merged make a link for input libdir to CountLibDir
  if (${$params}{"LibDir"} !~/:/)
  {

   if ( !(-s ${$params}{"LibDir"}))
   {
    my $com="ln -s ".${$params}{"InputLibDir"}." ".${$params}{"LibDir"};
    print $com."\n";
    `$com`;
   }
  }
}

sub getGenomeFile
{
  my ($allchrom, $genomedir, $chrom)=@_;
  open(F, "$genomedir/$chrom.fa") or die "Error opening: $!";
  while (my $line=<F>)
  {
    chomp($line);
    if ($line!~/>/)
    {
     ${$allchrom}.=$line;
    }
  }
}


sub getBGFile
{
 my ($vals, $c, $BGfile, $chr) = @_;
 open(BG, $BGfile);
 print ":::".$BGfile."\n";
 my $count=0;
 while(my $r=<BG>)
 {   
   chop($r);
   my @p=split(/[\t\s]+/, $r);
   if ($p[3]!~/^$/)
   {
    ${$vals}{$p[1]}=$p[3];
    $count++;    
   }   
 }
 ${$c}=$count;
 
}

####################################
## F U N C T I O N:  getBGVals   ##
####################################
## INPUT ARGUMENTS:
#
#  $BGfile: counts of chunks for each position
#  Each row of the bedGraph file is as
#  <choromosome number> <start position of the match> <end position of the match + 1> <#reads (or matchs)>
#
# $hastags: a reference to the hash that will hold the 

sub getBGVals
{
 my ($BGfile, $hashtags, $sorted, $totalsum, $totalcount) = @_;
 

 open(BG, $BGfile) or die("Can not open the bedGraph file $BGfile\nAborting...\n");

 while(my $r=<BG>)
 {   
   chomp($r);
   my @p=split(/[\t\s]+/, $r);
   
   
   #only get the lines having exactly 4 entries 
   if(scalar(@p) == 4)
   {
    ${$hashtags}{$p[1]}=$p[3];
    ${$totalsum}+=$p[3];
    ${$totalcount}+=1;
   }
 }
 #sort the start position of the reads in ascending order 
 @{$sorted} = sort { $a <=> $b } keys %{$hashtags};

 close(BG);

}


# open bedGraph File for Merge Operation
sub getBGValsForMerge
{
 my ($vals, $c, $BGfile, $filecount) = @_;
 open(BG, $BGfile);
 my $count=0;
 while(my $r=<BG>)
 {   
   chop($r);
   my @p=split(/[\t\s]+/, $r);
   if ($p[3]>2)
   {
    ${$vals}{$p[1]}{$filecount}=$p[3];
    $count++;    
   }
   
 }
 ${$c}=$count;
}


#############################################################
## FUNCTION: w l o g ######################################## 
# Writes the log files stamping each line with current time
# INPUT: 
### $LOG: File handlle of the logfile
### $text: Lie to be written
#############################################################

sub wlog
{
  my ($LOG, $txt) = @_;
  my $timestamp = localtime(); 
  print $LOG "$timestamp: $txt\n";
}

sub getRPKM
{
 my ($rpkmfile, $c, $hashRPKM, $sortedRPKM, $geneexon) = @_;
 open(RPKM, $rpkmfile);
 my $i=0;
 while(my $r=<RPKM>)
 {   
   chop($r);
   my @p=split(/[\t\s]+/, $r);

   ${$hashRPKM}{$p[0]."_".$p[1]."_".$p[2]}=$p[$geneexon];
   push(@{$sortedRPKM}, $p[0]."_".$p[1]."_".$p[2]);
   $i++;
 }
 ${$c}=$i;
}

#############################################################
## FUNCTION: g e t L a m d a ################################ 
##
## Read the contents of the lambdafile consisting of lambda values
## for each exon into the hash hashLambda.
## Note that $hashLambda is a reference to an actual hash. 
## 
## INPUT: $lambdafile, count reference c and a hash reference $hashLambda 
### Each row of the $lambdafile is: 
### chr_Gene_exon Lambda_Value Exon_RPKM Exon_Count Exon_Length
### Let %hashLambdaActual be the actual hash, referenced by $hashLambda, used in the system.
### Then the elements of it are of the form
### $hashLambdaActual{chr_gene_exon} = lambdaValue
###
### $c holds the number of lambda vallues read. 
#############################################################
sub getLambda
{
 my ($lambdafile, $c, $hashLambda) = @_;
 open(LAMBDA, $lambdafile) || die("Can not open the lambda file: $lambdafile");
 my $i=0;

 while(my $lam=<LAMBDA>)
 {   
   chop($lam);
   my @p=split(/[\t\s]+/, $lam);
   ${$hashLambda}{$p[0]} = $p[1];
   $i++;
 }
 ${$c}=$i;
 close(LAMBDA);
}

#############################################################
## FUNCTION: g e t R P K M B o t h  #########################
##
## Read the RPKM file into the hashes referenced by $hashRPKM, 
## $hashtransRPKM and the list referenced by $sortedRPKM 
##
## $c holds the line count
##
##  
sub getRPKMBoth
{
 my ($rpkmfile, $c, $hashRPKM, $hashtransRPKM, $sortedRPKM) = @_;
 open(RPKM, $rpkmfile) or die("Can not open the RPKM file: $rpkmfile\n Exiting.");
 my $i=0;
 while(my $r=<RPKM>)
 { 
   chop($r);
   my @p=split(/[\t\s]+/, $r);

   ${$hashRPKM}{$p[0]."_".$p[1]."_".$p[2]}=$p[5];
   ${$hashtransRPKM}{$p[0]."_".$p[1]}=$p[6];
   push(@{$sortedRPKM}, $p[0]."_".$p[1]."_".$p[2]);
  $i++;
 }
 ${$c}=$i;

  close(RPKM); 

}


#####################################################################
## FUNCTION: g e t P e a k R e g i o n s  ###########################
##
## It reads the bed file into the variable referenced by $regionhash.
## beddir contains the bed files for each region.
## So this function reads the region file for the given region $region
### 
### INPUT: 
### $beddir: The directory containing the bed files
### $chr: the specified chromosome
###
### Each line fo the region bed file contains
### ($chrom, $start, $end, $name, $strand)
### 
### These attributes are read into the $gene struct.
### All these $gene's are collected in the array @$genes.  
###
sub getPeakRegions
{
 my ($regionhash, $beddir, $region) = @_;
 my $regionfile=do { local(@ARGV, $/) = "$beddir/$region.bed"; <>};
 #$gene1=~s/\r//g;
 my @lines=split /[\n\r]/,$regionfile;
 foreach my $line (@lines)
 {
    my @ele=split /[\t\s]+/,$line;
    my  ($chrom, $start, $end, $name, undef, $strand)=@ele[0..5];
    my $reg=new Region;

    $reg->chrom($chrom);
    $reg->start($start);
    $reg->end($end);
    $reg->name($name);   
    $reg->strand($strand);
    push(@{$regionhash}, $reg);
 }
}

sub getGenes
{
  my ($genes, $refdir, $chr) = @_;

  my $refSeqFile =  $refdir."/".$chr.".ref";

  open(GENEHANDLE, $refSeqFile ) or die(" Error: Can not open the RefSeq file: $refSeqFile aborting.\n");

  my $i = 0;

  while(<GENEHANDLE>)
  {
    chomp($_);
    my @tempGene = split(/\s+/, $_);
    
    # There has to be 11 items or more coming from each line in the refseq file
    # We check this here and ignore the line otherwise
    if ($tempGene[0]=~/^N/ && @tempGene==11)
    {
      ${$genes}[$i++] = $_;
    }
    elsif($tempGene[0]=~/^\d+$/) 
    {
      # This refseq is not the format in the manual but it can be the format downloaded from ucsc, 
      # We are going to try to skip first and 10th column which we are not using them 
      # and will try to run this step accordingly.
      my $index=$i++;
      for (my $j=1; $j<11; $j++)
      {
        ${$genes}[$index] .= $tempGene[$j]."\t";
      } 
      ${$genes}[$index] .= $tempGene[12]."\t";
      
    }
  }
}

# "run" function manages the execution of the program according to the parameter runPBS.
# If runPBS is 1, then the program is executed in the cluster with parameters 
# If runPBS is 0, then the program is executed standalone.  
sub run
{
  my ($LOG, $cmd, $dir, $runPBS, $name, $queue) = @_;
  my $com=$cmd;

  if ($runPBS == 1)
  {
    $com="$dir/runPBS.pl -com \"$cmd\" -name $name -cpu 1 -q $queue";
  }
  else
  {
     $com = $cmd;
  }

  print "Executing: $com\n\n";
  wlog($LOG, $com);
  my $res=`$com`;
  wlog($LOG, $res);  
}

sub getRegionChrom
{
 my ($regionhash, $regionfile, $chr) = @_;
 open(BED, $regionfile) or die("\nError: Can not open the Bed file $regionfile\n");
 
 while(my $line=<BED>)
 {   
  chomp($line);
  my @ele=split /[\s\t]+/,$line;
  my  ($chrom, $start, $end, $name, undef, $strand)=@ele[0..5];
  if ($chrom eq $chr) {
      my $reg=new Region;
      $reg->chrom($chrom);
      $reg->start($start);
      $reg->end($end);
      $reg->name($name);   
      $reg->strand($strand);
      push(@{$regionhash}, $reg);
  }
 }
}

sub getBedFile
{
my ($bedfile, $arr) = @_;
open(BED, $bedfile) or die("\nError: Can not open the Bed file $bedfile\n");
 
while(my $r=<BED>)
{   
  chomp($r);
  my @p=split(/[\t\s]+/, $r);
  if(@p>1)
  {
    my $bed=new BEDFILE;
    $bed->start($p[1]);
    $bed->end($p[2]);
    push (@{$arr},$bed);
   }
 }
}

#######################################################
######## g e t Dir C h r o m o s o m e s  #############
########################################################
#### Gets the chromosomes by looking at the files at the given directory
# Input:
#         $BGdir: Files directory
#         $extension: Extension of the files
#
# Output:
#	  @chromosomes: A list containing the chromosomes
################################################################
sub getDirChromosomes()
{
  my ($BGDir, $extension) = @_;
  
  my $tailLength = length($extension);
  
  opendir(my $DIRHANDLE, $BGDir) or die("Error! Can not open the directory $BGDir.\n Aborting.\n");
  my @BGFiles = grep {(substr($_, -$tailLength) eq $extension)  && -f "$BGDir/$_" } readdir($DIRHANDLE);
  closedir($DIRHANDLE);
  
  my @chromosomes;
  
  foreach my $file (@BGFiles){
    if($file =~ /(chr[0-9a-zA-Z]+)/){
      push(@chromosomes, substr($file, 0, length($file)-$tailLength));
    }
  }
  
  return @chromosomes;
}

sub findBGTotal
{
  my ($BGFile) = @_;
  my $res=`wc -l $BGFile`;
  $res=~/^(\d+)/;
  my $total=$1;
  return $total;
}

sub writeJOBList
{
my ($parallel, $JOBLIST, $jobnum, $com, $outd, $script)=@_;
  if ($parallel)
  {
    makeDir("$outd/scripts/");
    makeBashScript("$outd/scripts", $script, $com);
    my $JOB;
    if (${$jobnum}==1) {
       open ($JOB, ">$JOBLIST");
    }
    else
    {
      open ($JOB, ">>$JOBLIST");
    }
    print $JOB ${$jobnum}."\t$outd/scripts/$script\n";
    close($JOB);
  }
  else
  {
     `$com`;
  }
  ${$jobnum}++;
}

sub writeJOBWaitLine
{
my ($parallel, $JOBLIST, $jobnum, $txt, $start, $end)=@_;
 if ($parallel)
 {
  if (!defined $start)
  {
    $start=0;
  }
  #print "$start\t$end\n";
  if (${$jobnum} !~/^$/ && ${$jobnum}>$start)
  {
    open ($JOB, ">>$JOBLIST");
    #print ${$jobnum}."\techo \"$txt\"\t".$start."\n"  ;
    print $JOB ${$jobnum}."\techo \"$txt\"\t".$start;
    for(my $i=($start+1); $i<$end; $i++)
    {
      print $JOB ",$i";
    }
    print $JOB "\n";
    close($JOB);
    ${$jobnum}++;
  }
 }
}
sub makeBashScript
{
my ($dir, $script, $cmd)=@_;
my $bindir=`pwd`;
chomp($bindir);
my $txt="
#!/bin/bash
## Job name
#PBS -N $script
## Declare job re-runnable
#PBSe-V
#PBS -S /bin/bash
## Output files
#PBS -o 'sge/$script.\$JOB_ID.log' -j y
#PBS -j y
#PBS -q serial

#PBS -l nodes=1:ppn=1

# Log Job Info for debugging
date '+TS[START]: \%Y-\%m-\%d \%k:\%M:\%S.\%N'
echo StartTime is `date`
#echo Working directory is \$PBS_O_WORKDIR
#cd \$PBS_O_WORKDIR
echo Directory is `pwd`
echo Running on host \$PBS_O_HOST / `hostname`
echo Job \$PBS_O_JOBID - \$PBS_JOBNAME in Queue \$PBS_QUEUE

## Force all files to be group writable
umask 0002
## And attempt to get output readable outside
unset DISPLAY

##==== PRE =========================================================
cd $bindir
##==================================================================
date '+TS[JOB_START]: \%Y-\%m-\%d \%k:\%M:\%S.\%N'
set +e
## =============== COMMAND IS RUNING HERE! =========================
$cmd

## Save exit status
RETVAL=\$?
set -e
date '+TS[JOB_END]: \%Y-\%m-\%d \%k:\%M:\%S.\%N'
##==== POST ========================================================
##==================================================================

echo EndTime is `date`
date '+TS[END]: \%Y-\%m-\%d \%k:\%M:\%S.\%N'

# return the executed command's exit status
exit \$RETVAL
";

open (OUT, ">$dir/$script");
print OUT "$txt";
close(OUT);
`chmod 755 $dir/$script`;
}


############################################################

return 1;

