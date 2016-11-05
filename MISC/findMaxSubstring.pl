#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use MyGeom;
use MyConfigs;
use PDB;
use ConfigPDB;
use Math::Geometry ;
use Math::Geometry::Planar;


use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($idx,$infile,$postfix,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$protein);
my ($ignorefile,@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "postfix=s"=>\$postfix ,
            "infile=s"=>\$infile ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "idx=i"=>\$idx ,
            "verbose=i"=>\$verbose ,
            "cutoff=f"=>\$cutoff ,
           );
require String::LCSS;

usage( "Need to give a protein pdb id -option -protein  ") if(!defined $protein);
usage( "Need to give a protein pdb id -option -outfile  ") if(!defined $outfile);

my $ofh = util_open_or_append($outfile);
my ($A,$firstlinea) = util_readfasta("FASTADIR_NT/$protein.ALL.1.fasta");
my ($B,$firstlineb) = util_readfasta("FASTADIR_NT/$protein.ALL.1.fasta.comp.fasta");


my ($longest,$DIFF) = ProcessFwdAndReverse($A,$B);
 
if($DIFF eq 0){
	print "Palindromic\n";
    my ($l1) = ($A =~ /(.*)$longest/);
    my $rev = util_getComplimentaryString($l1) ;
    my $DIFF = ProcessFwdAndReverse($l1,$rev);


}

sub ProcessFwdAndReverse{
	my ($a,$b) = @_ ;
   my $totlen = length($a);

  my @result = String::LCSS::lcss ( $a, $b);
  #my $longest = String::LCSS::lcss ( $a,$b);
  my $longest = $result[0];
  my $len = length($longest);
  my $rev = util_getComplimentaryString($longest) ;


  my ($l1) = ($a =~ /(.*)$longest/);
  my ($l2) = ($a =~ /(.*)$rev/);
  my $N1 = length($l1);
  my $N2 = length($l2);
  my $diff = ($N1 - $N2);
  print  "$protein $totlen $len $longest $rev $N1 $N2 $diff\n";
  return ($longest,$diff);
}

  #print "$result[0] ($result[1],$result[2])\n";

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
