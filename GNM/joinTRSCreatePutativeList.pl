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
my ($mergedir,$ignorefile,@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "mergedir=s"=>\$mergedir ,
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
usage( "Need to give a output file name => option -mergedir ") if(!defined $mergedir);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);

my $ofh = util_write($outfile);


foreach my $idx (0..100){
	##my $infile = "$mergedir/data.MERGED.OVERLAP". $idx ;
	my $infile = "$mergedir/data.MERGED.OVERLAP". $idx .".out";
	next if(! -e $infile);
    my $info = {};
    my $ifh = util_read($infile);
    while(<$ifh>){
		my @l = split ;
		die if(@l ne 4);
		my $a = $l[1];
		my $b = $l[2];
		print $ofh "$a\n";
		print $ofh "$b\n";
    }
    close($ifh);
}



sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
