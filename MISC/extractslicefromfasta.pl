#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use MyGeom;
use PDB;
use ConfigPDB;

use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($start,$end,$infile,$outfile,$which_tech,$listfile,$protein);
my (@expressions);
my $howmany ;
my $verbose = 1 ;
my $strict = 0 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "strict=i"=>\$strict ,
            "start=i"=>\$start ,
            "end=i"=>\$end ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_open_or_append($outfile);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
usage( "Need -start  ") if(!defined $start);
usage( "Need -end  ") if(!defined $end);


my $CNT = 0 ; 
my $info = {};


my ($str,$firstline) = util_readfasta($infile);

my $LEN = length($str);
print "len = $LEN , start = $start, end = $end \n";
chomp $firstline ;

if(!$strict){
	$start = 1 if($start < 0);
	$end = $LEN if($end > $LEN);
}



print $ofh  "${firstline}_$start--$end\n";


my $retstr = util_extractSliceFromFasta($str,$start,$end);
print $ofh "$retstr \n" ;

my $NNN = length($retstr);

my $diff = ($end - $start) + 1;
my @l = split "", $str ; 
my $N = @l ;
print "$N is the length of the whole, and you want $diff of it \n" if($verbose);



chmod 0777, $outfile ;

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
