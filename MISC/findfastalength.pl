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
my ($start,$end,$name,$infile,$outfile,$which_tech,$listfile,$protein);
my (@expressions);
my $howmany ;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "name=s"=>\$name ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "start=i"=>\$start ,
            "end=i"=>\$end ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);

my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP,$BLASTOUT,$BLASTDB) = util_SetEnvVars();
if(defined $outfile){
   usage( "Need to give a output file name => option -name ") if(!defined $name);
   $infile = "$FASTADIR/$name.ALL.1.fasta";
}
else{
   usage( "Need to give a output file name => option -infile ") if(!defined $infile);
}

my ($str,$firstline) = util_readfasta($infile);
my $len = length($str);


if(defined $outfile){
    my $ofh = util_append($outfile);
    print $ofh "$name $len\n";
}
else{
    print "$infile length=$len\n";
}

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
