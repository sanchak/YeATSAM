#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use MyGeom;
use PDB;
use ConfigPDB;
use Math::Geometry ;
use Math::Geometry::Planar;


use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($idx,$infile,$f1,$f2,$outfile,$cutoff,$which_tech,$listfile,$protein);
my ($ignorefile,@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "idx=i"=>\$idx ,
            "infile=s"=>\$infile ,
            "f1=s"=>\$f1 ,
            "f2=s"=>\$f2 ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "verbose=i"=>\$verbose ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);
usage( "Need to give a input file name => option -f1 ") if(!defined $f1);
usage( "Need to give a input file name => option -f2 ") if(!defined $f2);
usage( "Need to give a input file name => option -idx ") if(!defined $idx);

#my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
 my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP,$BLASTOUT,$BLASTDB) = util_SetEnvVars();

my ($list1,$mean,$sd,$assign) = util_ExtractOneIdxFromFile($f1,$idx,"junk");
my ($list2,$mean,$sd,$assign) = util_ExtractOneIdxFromFile($f2,$idx,"junk");

my @l1 = @{$list1};
my @l2 = @{$list2};
die if(@l1 ne @l2);


while(@l1){
	my $a = shift @l1 ;
	my $b = shift @l2 ;
	my $diff = util_format_float($a - $b,1) ;
	print $ofh "$diff = diff \n";
}

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
