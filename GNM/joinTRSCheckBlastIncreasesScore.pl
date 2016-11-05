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
my ($infile,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$blastdir);
my ($ignorefile,@expressions);
my $howmany = 100000 ;
my $verbose = 0 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "blastdir=s"=>\$blastdir ,
            "infile=s"=>\$infile ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
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
my $ofhdone = util_write($outfile. ".done");
my $ofhwrong = util_write($outfile. ".wrong");
my $ofhnotdone = util_write($outfile. ".notdone");
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
my $ifh = util_read($infile);

usage( "Need to give a blastdir pdb id -option -blastdir  ") if(!defined $blastdir);
my $CNT = 0 ; 
#my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
 my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP,$BLASTOUT,$BLASTDB) = util_SetEnvVars();

my $PWD = cwd;

my $ofhI = util_write("merge.include");
my $ofhE = util_write("merge.exclude");

my $info = {};
while(<$ifh>){
	chomp ;
     next if(/^\s*$/);
	 next if(/^\s*#/);
	 my $V = $_;

	 my ($a,$b,$m) = split ; 

	 my ($blastscorea) = parseSinglefile("$blastdir/$a.ann");
	 my ($blastscoreb) = parseSinglefile("$blastdir/$b.ann");
	 my ($blastscorem) = parseSinglefile("$blastdir/$m.ann");
	 next if(!($blastscorea || $blastscoreb || $blastscorem));
	 if(($blastscorea < 10 && $blastscoreb < 10) || !($blastscorea || $blastscoreb || $blastscorem)){
	      print $ofhnotdone "$V \n";
	      next ;
	 }

	 my $DIFF ;
	 my $AB = $blastscorea + $blastscoreb ;
	 my $PERCENT = .2 * $AB ;

	 ## if one match is very less, be strict about the whole match
	 if($blastscorea < 10 || $blastscoreb < 10){
	    $DIFF = $AB + $PERCENT ;
	 }
	 else{
	    $DIFF = $AB - $PERCENT ;
	 }


	 if($blastscorem > $DIFF ){
	      print " $V $blastscorea $blastscoreb $blastscorem $blastscorem $DIFF \n" if ($verbose);
	      print $ofhE "$a\n";
	      print $ofhE "$b\n";
	      print $ofhI "$m\n";
	      print $ofhdone "$V \n";
	 }
	 elsif($blastscorea > 60 && $blastscoreb > 60){
	      # this is actually a wrong merge - this means that the merged match does not improve 
	      print " Wroing $V $blastscorea $blastscoreb $blastscorem $blastscorem $DIFF \n" if($verbose);
	      print $ofhwrong "$V\n";
	 }
	 else{
	      print $ofhnotdone "$V \n";
	 }
}

sub parseSinglefile{
	my ($file) = @_ ; 
	my $ifh = util_read($file);
	my $blastscore  = 0 ;
	#print "$file \n";
	while(<$ifh>){
		next if(/Warning/);
		my @l = split ;
		my $N = @l - 1;
		$blastscore = $l[$N -1];
	}
	return $blastscore ;
}
print STDERR "Output written in $outfile\n";


chmod 0777, $outfile ;

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
