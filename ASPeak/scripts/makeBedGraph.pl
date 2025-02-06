#!/usr/bin/env perl

#########################################################################################
#                                       makeBedGraph.pl
#########################################################################################
# 
# This program converts bed files to BedGraph using the reads from bed files under the
# peak regions.
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
 my $outdir           = "";
 my $name             = "";
 my $help             = "";
 my $print_version    = "";
 my $version          = "1.0.0";
################### PARAMETER PARSING ####################

my $cmd=$0." ".join(" ",@ARGV); ####command line copy

GetOptions( 
	'input=s'        => \$input,
	'outdir=s'       => \$outdir,
	'name=s'         => \$name,
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

pod2usage( {'-verbose' => 0, '-exitval' => 1,} ) if ( ($input eq "") or ($outdir eq "") or ($name eq "") );	


################### MAIN PROGRAM ####################
#    Read input file calculate the # of reads centered each nucleotide and write it into output file

 
 open(IN, $input);
  
 my %bgvals=();
 my %diffisoforms=(); #Find all different isoforms for different center point. Divide read count to total isoforms for this position
 my $watcri="a";
 my $chrom="";
 my $old_chrom="";
 
 open(my $OUT, ">$outdir/$name.bg");
 open(my $OUTWATSON, ">$outdir/$name.watson.bg");
 open(my $OUTCRICK, ">$outdir/$name.crick.bg");
 while(my $line=<IN>)
 {   
   chomp($line);
   print $line."\n";
   my @a=split(/[\t\s]+/,$line );
   $chrom=$a[0];
   if ($old_chrom ne $chrom && $old_chrom !~ /^$/)
   {
       print "WRITING the chrom:$old_chrom\n";
       writeToFiles($old_chrom, \%bgvals, \%diffisoforms, $OUT, $OUTWATSON, $OUTCRICK);
       %bgvals=();
       %diffisoforms=();
   }
   for (my $pos=$a[1]; $pos<$a[2]; $pos++)
   {
      $watcri="w";
      if ($a[5] eq "-") {
	   $watcri="c";
      }
      #print $pos."\n";
      $bgvals{$pos}{$watcri}++;
      $diffisoforms{$pos}{$watcri}{$a[3]}=1;
   }
 
   $old_chrom=$chrom;
 }
       print "chrom\n";
writeToFiles($chrom, \%bgvals, \%diffisoforms, $OUT, $OUTWATSON, $OUTCRICK);
 
close(IN);
close($OUT);
close($OUTWATSON);
close($OUTCRICK);

sub  writeToFiles
{
 my ($chrom, $bgvals, $diffisoforms, $OUT, $OUTWATSON, $OUTCRICK)=@_;
    
 my @arr=("c", "w");
  
 foreach my $pos (sort{$a<=>$b} keys %{$bgvals})
 {
   my $totcount=0;
   foreach my $wca (@arr)
   { 
     if (exists ${$bgvals}{$pos}{$wca})
     {
       my $countpos=${$bgvals}{$pos}{$wca};
       my %diff=%{${$diffisoforms}{$pos}{$wca}};
       my $size = scalar keys %diff;
       my $count = floor($countpos/$size);
       my $O = $OUTWATSON;
       $O = $OUTCRICK if ($wca eq "c");
        
       print $O "$chrom\t$pos\t".($pos+1)."\t".$count."\n";
       $totcount+=$count;
     }
   }
   print $OUT "$chrom\t$pos\t".($pos+1)."\t".$totcount."\n";
  }
 }
  
__END__


=head1 NAME

makeBedGraph.pl

=head1 SYNOPSIS

makeBedGraph.pl -i input file <BED format> -o output file <BG format> -n name

makeBedGraph.pl -help

makeBedGraph.pl -version

For help, run this script with -help option.

=head1 OPTIONS

=head2 -i  input file <BED format> 

It converts whole bed file together.

=head2 -o output file <BG format>

BedGraph file that includes only called peak regions. 

=head2 -help

Display this documentation.

=head2 -version

Display the version

=head1 DESCRIPTION

 This program converts bed files to BedGraph using the reads from bed files under the
 peak regions.

=head1 EXAMPLE

./makeBedGraph.pl -p your_input_dir/bed_file.bed -o your_output_dir/bed_file.bg 

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
