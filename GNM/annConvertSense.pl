#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use MyGeom;
use PDB;
use ConfigPDB;
use MyGNM;
use Math::Geometry ;
use Math::Geometry::Planar;

## disable perl's warning mechanism
no warnings 'recursion';

use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($isNT,$findcharstring,$forWGS,$checkforsame,$infile,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$protein);
my ($ignorefile,@expressions,$trs,$blastdir,$oldblastcutoff);
my $strict = 1 ;
my $verbose = 0 ;

my $blastcutoff ;
my $percentlength ;
my $percentmatched ;
my $percentidentity ;
my $expectlimit ;

my $postfix = "blast.nt";
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "postfix=s"=>\$postfix ,
            "trs=s"=>\$trs ,
            "checkforsame"=>\$checkforsame ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "listfile=s"=>\$listfile ,
            "blastdir=s"=>\$blastdir ,
            "ignorefile=s"=>\$ignorefile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "blastcutoff=i"=>\$blastcutoff ,
            "oldblastcutoff=i"=>\$oldblastcutoff ,
            "strict=i"=>\$strict ,
            "forWGS=i"=>\$forWGS ,
            "verbose=i"=>\$verbose ,
            "percentlength=i"=>\$percentlength ,
            "percentmatched=i"=>\$percentmatched ,
            "percentidentity=i"=>\$percentidentity ,
            "findcharstring=i"=>\$findcharstring ,
            "isNT=i"=>\$isNT ,
            "expectlimit=f"=>\$expectlimit ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -infile ") if(!defined $infile);

$outfile = "$infile.reversesense";

my $ofh = util_write($outfile);
my $ifh = util_read($infile);

my $infomap = {};
while(<$ifh>){
     next if(/^\s*$/);
     next if(/^\s*#/);
	 die if(/only one/);

     my ($nm,@l) = split ;
	 foreach my $i (@l){
	 	 $infomap->{$i} = {} if(!defined $infomap->{$i});
		 $infomap->{$i}->{$nm} = 1 ;
	 }
}

foreach my $k (sort keys %{$infomap}){
	my $table = $infomap->{$k};
	my @l = (sort keys %{$table});
	print $ofh "$k @l \n";
}



sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
