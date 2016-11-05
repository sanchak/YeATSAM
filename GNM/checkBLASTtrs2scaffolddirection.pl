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
use MyGNM;
use Math::Geometry ;
use Math::Geometry::Planar;


use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($idx,$infile,$postfix,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$protein);
my ($blastdir,$ignorefile,@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "blastdir=s"=>\$blastdir ,
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
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);


if(defined $infile){
	my $getdir = GNM_DirectionBlast($infile);
	if($getdir eq 0){
        my $ofh = util_write($outfile);
	}
}
else{

my $ofh = util_write($outfile);
usage( "Need to give a output file name => option -blastdir ") if(!defined $blastdir);
my @files = <$blastdir/*.blast> ;
my $N = @files ;
print "There are $N blast ....\n";


my $CNTWrongdir = 0 ;
my $CNTNohits = 0 ;
foreach my $infile(@files){
    my $info = {};

	my $getdir = GNM_DirectionBlast($infile);
	if(!$getdir){
		print $ofh "$infile\n";
		$CNTWrongdir++;
	}
	elsif($getdir eq -1){
		print $ofh "$infile\n";
		$CNTNohits++;
	}
}
print "CNTNohits $CNTNohits $CNTWrongdir CNTWrongdir \n";
}



sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
