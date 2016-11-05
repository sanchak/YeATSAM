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
 no warnings "all";


use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($blastdir,$force,$sname,$justchecking,$checkforsame,$infile,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$protein);
my ($promoters,$scaffoldfasta,$ignorefile);
my @trss ;
my $howmany = 100000 ;
my $verbose = 1 ;

my $percentlength = 10;
my $percentmatched = 70;
my $percentidentity = 30;


GetOptions(
            "justchecking=i"=>\$justchecking ,
            "force"=>\$force ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "sname=s"=>\$sname ,
            "blastdir=s"=>\$blastdir ,
            "scaffoldfasta=s"=>\$scaffoldfasta ,
            "checkforsame"=>\$checkforsame ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            "outfile=s"=>\$outfile ,
            "promoters"=>\$promoters ,
            "trss=s"=>\@trss,
            "howmany=i"=>\$howmany ,
            "verbose=i"=>\$verbose ,
            "percentlength=i"=>\$percentlength ,
            "percentmatched=i"=>\$percentmatched ,
            "percentidentity=i"=>\$percentidentity ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);

my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP,$BLASTOUT,$BLASTDB) = util_SetEnvVars();
my $HARD_DIFFFORMATCH = 5;
my $HARD_LONGDISTANCEFACTOR = 10;


usage( "Need to give a input file name => option -sname ") if(!defined $sname);
usage( "Need to give a input file name => option -blastdir ") if(!defined $blastdir);
usage( "Need to give a input file name => option -justchecking ") if(!defined $justchecking);
usage( "Need to give a input file name => option -scaffoldfasta ") if(!defined $scaffoldfasta);
usage( "Need to give a input file name => option -trss ") if(!@trss);


if(!$justchecking){
    $outfile = "INFO/$sname.info";
    die "Already done $outfile" if(-e $outfile && !defined $force);
}

my $ofh = util_write($outfile);

print "Processing $sname \n";





my ($strScaffold) = util_readfasta($scaffoldfasta);
my $scafflen = length($strScaffold);


my $ERRSTATE = 0 ;
while(@trss){
   my $trs = shift @trss ;
   my $infile = "$blastdir/$trs.$sname.blast";
   if(! -e $infile){
	warn " $infile does not exist ";
	$ERRSTATE = 1 ;
	next ;
   }

   my $str = GNM_MapOneTRS2Scaffold($trs,$sname,$infile,$strScaffold,$scafflen,$ofh,$justchecking,$HARD_DIFFFORMATCH,$HARD_LONGDISTANCEFACTOR);

   if(!defined $str){
       if( -e $infile){
	    $ERRSTATE = 1 ;
   	    system("unlink $infile");
       }
   }

}
if($ERRSTATE){
       if(-e $outfile){
   	     system("unlink $outfile");
       }
}
else{
	print " All well for $sname \n";
}


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}

