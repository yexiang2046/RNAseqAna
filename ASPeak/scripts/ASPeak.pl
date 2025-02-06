#!/usr/bin/env perl
#                                       ASPeak.pl
#########################################################################################
# 
# This program finds the peaks in a given Clip or Rip Seq data.
# This is a wrapper script for the system. 
# It processes the command line arguments, the parameter file and 
# calls the other scripts accordingly.  
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
  eval{require Math::CDF;};
 
 if($@)
 {
    print("Please install the perl module Math::CDF\n");
    exit;
 }
}

# This block executes before compiling the rest of the script.
# It finds the library file paths. 
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
 use Scalar::Util qw(looks_like_number);
 use strict;
 use lib ".";
 require lib::io;
 require lib::func;
#################### STRUCTS  #######################
struct( 'CommandParams', {
   LibName                 => '$',
   Lib                     => '$',
   LibFormat               => '$',
   ControlName             => '$',
   Control                 => '$',
   ControlFormat           => '$',
   RNASeqName              => '$',
   RNASeq                  => '$',
   RNASeqFormat            => '$',
   BedDir                  => '$',
   Region                  => '$', 
   OutDir                  => '$',
   Strand                  => '$',
   Nodes                   => '$',
   Queue                   => '$',
   RunInCluster            => '$',
   nornaseq                => '$',
   estimationRadius        => '$',
   gapNumber               => '$',
} );

#################### CONSTANTS ###################### 
my $dir=`pwd`;
chomp($dir);
my @types= ("StepCheckParams", "StepGetChroms", "StepConvert2Bed", "StepCountMappedReads", "StepSortBed", "StepRegionSeparation", "StepChromSeparation", "StepPrepInput", "StepMakeRPKM", "StepNBParameters", "StepPeak", "StepCombinePeaks", "StepCalcFDR", "StepPeaks2Bed", "StepMakeBedGraph", "StepWigToBigWig", "StepSubmitJobs");

#################### VARIABLES ###################### 
## Get command line options and initialize values
my (
	$paramfile, #parameter file
	$type,
	$libname,
	$lib,
	$libformat,
	$controlname,
	$control,
	$controlformat,
	$rnaseqname,
	$rnaseq,
	$rnaseqformat,
	$strand,
	$beddir,
	$region,
	$cluster,
	$nodes,
	$queue,
	$outdir,
	$nornaseq,
	$jobnum,
	$estimationRadius,
	$gapNumber,
    	$help,
	$print_version,
);

my $VERSION = '2.0.0';
################### PARAMETER PARSING ####################
my @pars = @ARGV;

### Quick help
unless (@ARGV) { 
	# when no command line options are present
	# print SYNOPSIS
	pod2usage( {
		'-verbose' => 0, 
		'-exitval' => 1,
	} );
}

# Command line options
GetOptions( 
	'param=s'            => \$paramfile,         # parameter file
	'type=s'             => \$type,              # run type
	'libname=s'          => \$libname,           # libray name/s 
	'lib=s'              => \$lib,               # library file/s 
	'libformat=s'        => \$libformat,         # (sam, bam, bed, or bowout format accepted as inputs)
	'controlname=s'      => \$controlname,       # control name/s 
	'control=s'          => \$control,           # control file 
	'controlformat=s'    => \$controlformat,     # (sam, bam, bed, or bowout format accepted as inputs)
	'rnaseqname=s'       => \$rnaseqname,        # rnaseq file identifier
	'rnaseq=s'           => \$rnaseq,            # rnaseq file 
	'rnaseqformat=s'     => \$rnaseqformat,      # (sam, bam, bed, or bowout format accepted as inputs)
	'strand'             => \$strand,            # Strand specific peak calling
	'cluster'            => \$cluster,           # If it is going to be submitted to the cluster     
	'nodes=s'            => \$nodes,             # The number of max processors/cores/nodes that can run at the same time in a cluster
	'queue=s'            => \$queue,             # The name of the queue
	'beddir=s'           => \$beddir,            # Bed Directory for region annotations
	'region=s'           => \$region,            # Region Name
	'outdir=s'           => \$outdir,            # Output directory
	'nornaseq'           => \$nornaseq,          # If you don't have RNASeq data set nornaseq
	'jobnum=s'           => \$jobnum,            # Job start num
	'estimationradius=i' => \$estimationRadius,  # Parameter estimation radius.
	'gapnumber=s'        => \$gapNumber,         # Number of gaps tolerated in peaks
	'help'               => \$help,              # request help
	'version'            => \$print_version,     # print the version
) or die " unrecognized option(s)!! please refer to the help documentation\n\n";

# Print version
if ($print_version) {
	print "Pipeline main script, version $VERSION\n\n";
	exit;
}

# Print help
if ($help) {
	# print entire POD
	pod2usage( {
		'-verbose' => 2,
		'-exitval' => 1,
	} );
}

### Check for requirements
# check param file
if ($paramfile=~/^$/ && $lib=~/^$/) {
   $paramfile = shift @ARGV or
		die "  OOPS! No library file specified! \n use --help\n";

}

# This block resolves the relative path issues.
  my($mainScriptName, $scriptDirectory, $suffix) = fileparse(abs_path($0)); # Get the directory of the working script
  my @scriptPathPieces = split(/\//,$scriptDirectory); # Make an array of subdirectories
  pop(@scriptPathPieces); # Since this script is under "/dir1/dir2/.../dirN/scripts" directory, it pops "scripts" from the array
  my $mainDir = join("/", @scriptPathPieces); #This is the main directory containing ASPeak files. Hence  it is our library patth 


unless ($paramfile) {
   $paramfile =  $mainDir."/scripts/default.txt";
}

#print "PARAMFILE: $paramfile\n";

# set type to ""
unless ($type) {
	$type = "";
}
 

#    Read ParamFile 
 my %params=();
 my $logname="ASPeak";
 if ($type ne "")
 { 
   $logname.=".$type";
 }
 
 my $commandParams=new CommandParams;
 if ($libname!~/^$/){ $commandParams->LibName($libname)};
 if ($lib!~/^$/){ $commandParams->Lib($lib)};
 if ($libformat!~/^$/){ $commandParams->LibFormat($libformat)};
 if ($controlname!~/^$/){ $commandParams->ControlName($controlname)};
 if ($control!~/^$/){ $commandParams->Control($control)};
 if ($controlformat!~/^$/){ $commandParams->ControlFormat($controlformat)};
 if ($rnaseqname!~/^$/){ $commandParams->RNASeq($rnaseqname)};
 if ($rnaseq!~/^$/){ $commandParams->RNASeq($rnaseq)};
 if ($rnaseqformat!~/^$/){ $commandParams->RNASeqFormat($rnaseqformat)};
 if ($beddir!~/^$/){ $commandParams->BedDir($beddir)};
 if ($region!~/^$/){ $commandParams->Region($region)};
 if ($nodes!~/^$/){ $commandParams->Nodes($nodes)};
 if ($queue!~/^$/){ $commandParams->Queue($queue)};
 if ($strand){ $commandParams->Strand(1)};
 if ($nornaseq){ $commandParams->nornaseq(1)};
 if ($cluster){ $commandParams->RunInCluster(1)};
 if ($estimationRadius > 0){ $commandParams->estimationRadius($estimationRadius)};
 if ($gapNumber !~ /^$/ and looks_like_number($gapNumber)){ $commandParams->gapNumber($gapNumber)};
 
 $outdir=~s/\/mainfiles//g;
 if ($outdir!~/^$/){ $commandParams->OutDir($outdir)};
 
 lib::io::readPipeParams(\%params, $paramfile, $commandParams);
 
#  print "The main directory is ", $params{"MainDir"}, "\n";
#  print "The executable directory is ", $params{"ExecDir"}, "\n";
#  exit;
 
 my $outd=$params{"OutDir"};
 my $logdir=$params{"LogDir"};
 
 # if there is no value for the gapNumber coming from the command line arguments
 # set it to its default value in the parameter file
 # if the value in the paramete file is not ok, set it to 2
 unless( $gapNumber !~ /^$/ and looks_like_number($gapNumber) )
 {
  $gapNumber = ( (exists $params{"gapNumber"}) and ( looks_like_number( $params{"gapNumber"} ) ) and (  $params{"gapNumber"} >= 0 ) )  ? $params{"gapNumber"} : 2; 
 }
 
 
 lib::io::makeDir($logdir);
 # set jobnum to ""
 unless ($jobnum) {
	$jobnum = 1;
 }

 my $num=1;
 if (-s "$logdir/$logname.1")
 {
   $num = `ls $logdir/$logname.*|awk '{split(\$1, a, "."); print a[3]; }'|sort -n|tail -n 1`;
   $num++;
 }
 open(my $LOG, ">$logdir/$logname.$num");
  
 my $cmd=$0;

 lib::io::wlog($LOG, "program started at: ");    
 lib::io::wlog($LOG, "The command used is:");
 lib::io::wlog($LOG, "$cmd\n");
 lib::io::wlog($LOG, "Parameters specified:");
 lib::io::wlog($LOG, "input paramfile: $paramfile");
 lib::io::wlog($LOG, "type: $type");
 if ($type eq "" || $type eq "auto")
 {
  lib::io::writeParams($LOG, \%params); 
  lib::io::writeExp($cmd, \@pars, $paramfile, \%params);
 }
 lib::io::wlog($LOG, "Log Dir: $logdir");

################### MAIN PROGRAM ####################
#    
# It runs selected job types in parameters file or given type with --type option.
#    


my $repeatRun=$cmd." --param $paramfile";

if (!$type)
{
   if (-s "$outd/JOBLIST.txt")
   {
     `rm $outd/JOBLIST.txt`;
   }
   my $comparams="";
   $comparams=" -s " if ($strand);
   $comparams.=" --param $paramfile" if ($paramfile);
   $comparams.=" --nodes $nodes" if ($nodes);
   $comparams.=" -re $region" if ($region);
   $comparams.=" --queue $queue" if ($queue);
   $comparams.=" --outdir $outdir" if ($outdir);
   $comparams.=" --lib $lib" if ($lib);
   $comparams.=" --libname $libname" if ($libname);
   $comparams.=" --libformat $libformat" if ($libformat);
   $comparams.=" --control $control" if ($control);
   $comparams.=" --controlname $controlname" if ($controlname);
   $comparams.=" --controlformat $controlformat" if ($controlformat);
   $comparams.=" --rnaseq $rnaseq" if ($rnaseq);
   $comparams.=" --rnaseqname $rnaseqname" if ($rnaseqname);
   $comparams.=" --rnaseqformat $rnaseqformat" if ($rnaseqformat);
   $comparams.=" --beddir $beddir" if ($beddir);
   $comparams.=" --nornaseq" if ($nornaseq);
   $comparams.=" --cluster" if ($cluster);
   $comparams.=" --estimationRadius $estimationRadius" if ($estimationRadius > 0);
   $comparams.=" --gapnumber $gapNumber" if ($gapNumber >= 0);
   
   foreach my $t(@types)
   {
     if (exists $params{$t} && $params{$t}>0)
     {
          # set jobnum to 1
          unless ($jobnum) {
          $jobnum = 1;
       }
       my $com=$cmd." $comparams --type $t --jobnum $jobnum";
       #print $com."\n" if ($params{"Debug"});
       if ($params{"Debug"}) {
         print "Running $t...\n";
       }       
       lib::io::wlog($LOG,$com);
       my $jobout=`$com`;
       #print $jobout."\n"  if ($params{"Debug"});
       
       $jobout=~/jobnum=(\d+)\n/;
       $jobnum=$1;
     }
   }
}
    
lib::func::runSteps($LOG, \%params, \$jobnum, $type);
lib::io::wlog($LOG, "program ended");
close($LOG);

if ($params{"RunInCluster"})
{
   print "jobnum=$jobnum\n";
}

######################################################################################
# POD DOCUMENTATION
#######################################################################################

=head1 NAME

ASPeak.pl

=head1 SYNOPSIS

ASPeak.pl -lib <RIP/CLIP-Seq library file> -rnaseq <rnaseq file> -beddir <annotation directory> -outdir <output directory> 

Alternative run with optional command line parameters;

ASPeak.pl -type <StepCheckParams| StepGetChroms| StepSeparate| StepCount| StepRNASeq| StepLambda| StepPeak| StepMakeBedGraph| StepSubmitJobs>
		-param <parameterfile> -nornaseq 
		-libname <libname> -lib <libraryfile> -libformat <bam|sam|bowout|bed>
		-controlname <controlname> -control <controlfile> -controlformat <bam|sam|bowout|bed>
		-rnaseqname <rnaseqname> -rnaseq <rnaseqfile> -rnaseqformat <bam|sam|bowout|bed>
		-outdir <output directory>
		-beddir <bed directory for region annotations>
		-region <region name> [-gapnumber maxtolaratedgaps]
              
ASPeak.pl -help
              
ASPeak.pl -version

=head1 OPTIONS

=head2 -param <parameter file> [optional]

The parameter file contains system parameters and the paths of input and output documents and folders. If all optional parameters are defined in the parameters file, the other parameters can be optional.

=head2 -type <run type> [optional]
 type is on ptional parameter. If it is empty it executes all the steps that set to 1 in parameter file.  
 StepCount        : Finds counts for each exon in the CLIP / RIP Seq data.
 StepRNASeq       : Finds the RPKM vaues of the exons in the RNASeq data.
 StepLambda       : Calculates lambda values of the exons.
 StepPeak         : Finds peaks.
 StepMakeBedGraph : Makes wiggle  plots of the found peaks in bedGraph format.

=head2 -libname library name [optional]

Library identifier. If there are more than one library name, please put ":" in between the names. Example: lib1:lib2
 
=head2 -lib library file [optional]

The library file is given with this parameter. This is optional parameter.  If there are more than one library file, please put ":" in between filenames.
Please enter filenames with their full path. If you defined this parameter in parameter file, you don't have to give this parameter from command line.
The parameters given by command prompt always overwrite to the parameters red from parameter file. Example: Lib/Full/Path/lib1.sam:Lib/Full/Path/Lib2.sam


=head2 -libformat library file format <bam|sam|bowout|bed> [optional]

File format has to be given either in the parameter file or from the command prompt
We accept bam, sam, bowtie alignment files(please make the extension bowout) and bed file as inputs. Example: [sam:sam] or only sam if both are the same. 

=head2 -controlname control library name [optional]

Control library identifier. If there are more than one control name, please put ":" in between the names. Example: control1:control2

=head2 -control control file [optional]

The control file is given with this parameter. This is optional parameter. If there are more than one control file, please put ":" in between filenames.
Please enter filenames with their full path.
If you defined this parameter in parameter file, you don't have to give this parameter from command line.
The parameters given by command prompt always overwrite to the parameters red from parameter file.


=head2 -controlformat control file format <bam|sam|bowout|bed> [optional]

File format has to be given either in the parameter file or from the command prompt
We accept bam, sam, bowtie alignment files(please make the extension bowout) and bed file as inputs.  Example: [bed:bed] or only bed if both are the same.

=head2 -rnaseqname rnaseq library name [optional]

RNASeq library identifier. If there are more than one RNASeq, please put ":" in between the identifiers. Example: RNASeq1:RNASeq2


=head2 -rnaseq rnaseq file 

The rnaseq file is given with this parameter. This is optional parameter. If there are more than one rnaseq file, please put ":" in between filenames.
Please enter filenames with their full path.
If you defined this parameter in parameter file, you don't have to give this parameter from command line.
The parameters given by command prompt always overwrite to the parameters red from parameter file.

=head2 -rnaseqformat library file format <bam|sam|bowout|bed>

File format has to be given either in the parameter file or from the command prompt
We accept bam, sam, bowtie alignment files(please make the extension bowout) and bed file as inputs. Example: [bowout:bowout] or only bed if both are the same.

=head2 -outdir output directory [optional]

All output will be generated under this directory

=head2 -beddir annotation files for each region [optional]

Annotation files are in bed 6 format can be downloaded from ucsc genome browser. Please consult manual to learn how to download a exonic, intronic or any other regions. 

=head2 -region region name [optional]

Region Annotation: The name of the region has to be filename without extention and path. Only this region will be run when it is entered from the command line.

=head2 -help

Displays this documentation.

=head2 -version

Prints the version  number.

=head1 DESCRIPTION

This program finds the protein binding sites in a CLIP- or RIP-Seq experiment.
This is done by calling the peaks in the CLIP / RIP Seq data.
It utilizes an abundance sensitive peak calling algorithm.
More explicitly, it uses a statistical test to distinguish the peaks from the background
noise (possibly coming from experimental errors, sequencing errors and etc.) in the CLIP / RIP Seq data. 
This test uses the poisson distribution as the statistic.
This program estimates the lambda value of the poisson distribution for each exon by averaging the counts of
exons (in the RIP Seq data)  having similar RPKM values (in the RNASeq data).
Here we define the abundance of an exon as its RPKM value in the RNA Seq data.

Calling the peaks in the CLIP / RIP Seq data is done in several steps.
This is a wrapper script for the whole system. 
It processes the command line arguments, the parameter file and 
calls other scripts accordingly.  

This program requires a run type and a parameter file.
The parameter file contains the parameters and location of the input and output files.
You should have received a sample parameter file with ASPeak.
You can modify it for your data set.

For most users, we recommend using this script in the auto mode.  

     [Script path]/ASPeak.pl -type auto -param [File Path]/parameters.txt

Then you can check the status by running this script in the status mode.

     [Script path]/ASPeak.pl -type status -param [File Path]/parameters.txt

=head2 Input

=head3 Parameter file

=head3 BedDir data (Region Annotations with BED format) 

=head3 CLIP / RIP-Seq data

=head3 RNASeq data
     
=head2 Output

=head3 Peaks

The output file has the peak information whose columns, from left to right, are like below.

	      Column	      		
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

=head3 Wiggle Plots

When used in auto mode or StepMakeBedGraph mode, i.e., "type StepMakeBedGraph",
peaks are also plotted in bedGraph files for visualization purposes.
    
    Column 1: Chromosome
    Column 2: Start position of the peak
    Column 3: End position of the peak
    Column 4: Average tag count of the peak

=head1 Running ASPeak    
    
Have your Ref-Seq, CLIP / RIP-Seq and RNASeq data ready.
For the data formats, please see the user manual provided with ASPeak.
Note that CLIP / RIP-Seq and RNASeq data must be in bedGraph format.
Positions of the the reads must be in ascending order.
Also, files must be separated according to the chromosomes and named accordingly.
Namely, each file should hold exactly one chromosome data and 
should be given the name of the chromosome with the ".bg" extension.
(e.g. chr2.bg, chrY.bg).

Modify the parameter file according to your needs.
Change the parameter "RefSeq" to the path of the Ref-Seq data.
Change the parameter "WigDir" to the path of the CLIP / RIP Seq data.
Change the parameter "RNASeq" to the parameter of the RNA Seq data. 

Make sure that the parameter "nornaseq" is set to 0. If you have RNASeq data.

=head2 Auto Mode

For most users, we recommend running ASPEak without any type.
That is

     [Path of the script]/ASPeak.pl -param parameterFile

Especially with large datasets, it may take a while (a couple of hours or days depending on your system)
to finish. 

     
=head2 Manual Mode

It is also possible to run ASPEak in a step-by-step basis manually.
Note that the order of the steps is important as input of a later step is output of former steps.
Before proceeding to the next step, make sure that the current step has finished successfully.
You can check this by reading the log files, in the output directory, and seeing the output files.

After going through the above preparation, begin with the count run type.
This will find the exon counts in the CLIP / RIP seq data. 

     [Path of the script]/ASPeak.pl -type StepCount -param parameterFile
     
Next, run it with the type rpkm.
This will find the RPKM values of the exons in the RNASeq data. 

     [Path of the script]/ASPeak.pl -type StepRNASeq -param parameterFile
     
Then run with lambda parameter.
This will make the lambda file which holds the lambda values of the exons.
     
     [Path of the script]/ASPeak.pl -type StepLambda -param parameterFile
    
Now everything is ready to call the peaks. This is done with the peak parameter.  
     
     [Path of the script]/pipline.pl -type StepPeak -param parameterFile
     
If you wish, you can obtain wiggle plots of these peaks in bedGraph format 
     
     [Path of the script]/ASPeak.pl -type StepMakeBedGraph -param parameterFile

=head2 Using ASPeak without RNASeq Data

Though it is strongly recommended to use ASPeak with  RNASeq data,
it can also be used without RNASeq data.

Please be aware that, without RNASeq data, lambda values will be determined locally
from the CLIP / RIP-Seq data. Consequently, you will loose the advantage of 
finding the peaks using the abundance sensitive algorithm of ASPeak.

Prepare your data and parameter file.
In the parameter file, set the parameters for corresponding regions to 0.

Then you can run ASPeak by

     [Path of the script]/ASPeak.pl -param parameterFile

Or, if you want to run it manually, then skip the steps: count, rpkm and lambda.
Just start  with the peak step.

     [Path of the script]/ASPeak.pl -type peak -param parameterFile
     
If you wish, you can obtain wiggle plots of these peaks in bedGraph format
     
     [Path of the script]/ASPeak.pl -type MakeBedGraph -param parameterFile

     
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


__END__

