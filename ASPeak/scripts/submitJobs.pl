#!/usr/bin/env perl
#########################################################################################
#                                      job_submitter.pl
#########################################################################################
#
# This program submit the jobs in JOB_LIST.txt file. S 
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
 use POSIX;
 use Class::Struct;
 use File::Basename;
 use Getopt::Long;
 use Pod::Usage;
 use strict;
 require lib::io;
 require lib::func;

#################### CONSTANTS ###################### 
my $VERSION = '1.0.0';

#################### VARIABLES ######################
my $file          = "";
my $set           = "";
my $pe            = 1; #The number  nodes in a machine
my $q             = "";
my $processors    = 100; #Maximum number of nodes that are going to be used in a cluster at a time.
my $help          = "";
my $print_version = "";
my $pe_options    = "";
my $sge           = "./sge";
my $tmp           = "./tmp";


############## LIBRARIES AND PRAGMAS ################RAMETER PARSING ####################
my $cmd = $0." ".join(" ",@ARGV); ####command line copy 

GetOptions(
        'j=s'       => \$file, #
        'set=s'     => \$set,
        'p=s'       => \$pe,
	'n=s'       => \$processors,
        'q=s'       => \$q,
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

pod2usage( {'-verbose' => 0, '-exitval' => 1,} ) if ( ($file eq "") or ($set eq ""));

if ($pe eq "") {$pe=1;}

# let's create the location to store the output from SGE
lib::io::makeDir("$sge");
lib::io::makeDir("$tmp");


open (IN,$file) || die "Please specify a file: $!\n";
my $x=0;
my @qsub=();
my $hold_jobs="";
while(my $line=<IN>)
{
   chomp($line);
   #print $line."\n";
   $line=~s/[\s\t]+/ /gi;
   $line=~s/^[\t\s]+//;
   my @arr=split(/[\t\s]+/, $line);

   if ($line=~/echo/)
   {
      my @arr2=split(/"/, $line);
      if (@arr2!=3)
      {
	print "Error on input:\n$line: malformed quotes \"\n";
        exit;
      }
      $qsub[$x]="qsub -hold_jid $x";
      $arr2[2] =~ s/^[\s]+//g;
      @arr=split(/[\t\s]+/,$arr2[2]);
      $hold_jobs=$arr[0];
      $qsub[$x]="NORUN_NEEDED qsub -hold_jid $arr[0]";
      if ($arr[1])
      {
         $qsub[$x].=",$arr[1]";
         $hold_jobs.=",$arr[1]";
      }
    }
    else
    {
        if ($hold_jobs)
        {
          $qsub[$x]="qsub -cwd -V -hold_jid $hold_jobs $arr[1]";
        }	
        else
	{
	   $qsub[$x]="$arr[1]";
	}
     }
   $x++;
}

my $y=1;
my @qsub_job_number=();
foreach my $jobline (@qsub)
{
  # parse out the command to execute and run it
  if ($jobline!~/NORUN_NEEDED/)
  {
    my $command=$jobline; 
    my @job_name_split=split(/step_/, $jobline);
    my $job_full_name=$job_name_split[1];
    $job_full_name=~ s/\.sh//g;
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
    my $tod=$mon."_".$year."_".$yday; $tod.="_$hour$min$sec";
    if (!$job_full_name) { $job_full_name=$set."_$y";}
    #if ($q eq "ANY")
    #{
    #  my $com="qstat -g c|grep chassis|awk \'{print \$1\"\\t\"\$4}\'|sort -k2,2nr|head -n 1|awk \'{print \$1}\'";
    #  #print $com;
    #  $q=`$com`;
      
    #  chomp($q);
    #}
    #print "QUEUE:$q\n";
    
    my $queue="";
    
    if ($q ne "" && $q ne "ANY"){$queue="-q $q";}
    my $pipe="qsub $pe_options -cwd -V -N $job_full_name $queue -o $sge/sge.$job_full_name"."_$tod -j y -S /bin/bash";
    if ($jobline =~/-hold_jid/)
    {
      my @job_split=split(/[\s\t]+/, $jobline);
      my $ids=$job_split[4];
      my $ec=$job_split[5];
      $pipe.=" -hold_jid ";
      my @jids=split(/\,/,$ids);
      foreach my $id (@jids)
      {
        $pipe.=$qsub_job_number[$id].",";
      }
      $pipe=~ s/, +/ /g;
      $pipe=~ s/,$//g;
      $command=" ".$ec;
     }
     print "Attempting to submit job\n\t$pipe $command\n";
     my $a=`qstat|wc -l`;
     chomp($a);  
     while ($a>$processors)
     {
       $a=`qstat|wc -l`;
       sleep(1);
       chomp($a);  
       print "[".$a."]";
     }

     open(SGE_QSUB, "|$pipe > $tmp/$$.out 2>&1") || die "Can't open the pipe to submit jobs to";
     print SGE_QSUB $command, "\n";
     close SGE_QSUB;
     open (IN, "$tmp/$$.out") || die "Can't open $tmp/$$.out for QSUBd job, even though everything appeared to work fine";
     my $line=<IN>;
     close IN;
     $line =~ /Your job (\d+)/i;
     my $jobnumber=$1;
     if ($jobnumber) 
     {
       unlink("$tmp/$$.out");
     }
     else 
     {
  	print STDERR "WARNING: No job number. This message was received:\n$line";
		die "Error on submission or Job number";
     }
     $qsub_job_number[$y]=$jobnumber;
     print "\t - Job submitted successfully as job id# $jobnumber\n\n";
     open(FILE_OUTPUT,">>$sge/sge.$job_full_name\_$tod");
     print FILE_OUTPUT "Starting $job_full_name with Job ID $jobnumber\n";
     print FILE_OUTPUT "Command submitted\n\t$pipe $command\n";
     print FILE_OUTPUT "\t - Job submitted successfully as job id# $jobnumber\n\n";
     close(FILE_OUTPUT);
  }
  $y++;
}

__END__
	
=head1 NAME

jobSubmitter.pl

=head1 SYNOPSIS

jobSubmitter.pl -j <JOBlist> -s <library name> 

jobSubmitter.pl -help

jobSubmitter.pl -version  

For help, run this script with -help option.

=head1 OPTIONS

=head2 -j <JOBlist>

It contains the jobs that are going to be run on a cluster system.

=head2 -s <library name>
 
The name of the library.

=head2 [-p procs]

The number of processors.

=head2 [-q queue]

The name of the queue.

=head2 -help

Display this documentation.

=head2 -version

Display the version

=head1 DESCRIPTION

It submits the jobs in JOBLIST.txt file.

=head2 Input

The input is the JOBLIST.txt file.
=over 6
File contents;

line_number {spaces or tab} script_name

and or

line_number {spaces or tab} echo "wait for my job to finish" JOB_#s seperated by comma's to wait for before running this job

File example:

1       ./scripts/step_count.sh

2       echo "wait for matching jobs to finish" 1

3       ./scripts/step_lambda.sh

4       echo "wait for matching jobs to finish" 1,3

5       ./scripts/step_peakchr1.sh

6       ./scripts/step_peakchr2.sh

7       echo "wait for post_matching_by_chr jobs to finish"     1,3,5,6


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
