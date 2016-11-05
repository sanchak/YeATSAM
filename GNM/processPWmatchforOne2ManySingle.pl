#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use MyGeom;
use MyGNM;
use PDB;
use ConfigPDB;
use Math::Geometry ;
use Math::Geometry::Planar;


use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($tag,$genomefile,$infile,$f1,$f2,$outfile,$cutoff,$getcommands,$which_tech,$listfile,$protein);
my ($isNT,$outdir,$ignorefile,@expressions);
my $howmany = 100000 ;
my $verbose = 0 ;
my $maxscore = 200 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "f1=s"=>\$f1 ,
            "f2=s"=>\$f2 ,
            "listfile=s"=>\$listfile ,
            "genomefile=s"=>\$genomefile ,
            "ignorefile=s"=>\$ignorefile ,
            "tag=s"=>\$tag ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "maxscore=i"=>\$maxscore ,
            "isNT=i"=>\$isNT ,
            "getcommands=i"=>\$getcommands ,
            "verbose=i"=>\$verbose ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);

usage( "Need to give a listfile -option -f1  ") if(!defined $f1);
usage( "Need to give a listfile -option -tag  ") if(!defined $tag);
usage( "Need to give a listfile -option -isNT  ") if(!defined $isNT);
usage( "Need to give a outfile pdb id -option -outfile  ") if(!defined $outfile);
my $CNT = 0 ; 
#my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
 my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP,$BLASTOUT,$BLASTDB) = util_SetEnvVars();


my $BLASTCOMMAND = $isNT ? "blastn" : "blastp";
my $map2name = {};

my $N = @expressions ;

my $cnt = 0 ;
my $bestid ;
my $bestscore = 0;
my $foundGREAT = 0 ;

my $lastexpr ; 
foreach my $expr (@expressions){
	$lastexpr = $expr ;

	$cnt++;

	my $commandstring = "$BLASTCOMMAND -query $f1 -subject $expr -out $outfile";
	print "command = $commandstring\n"  if($verbose);
	system("$commandstring");
	my ($info,$querylength,$Subjectname,$queryname,$blastscore,$expect) = GNM_PARSEBLAST_BESTVAL($outfile);
	die "- querylength not defined in $outfile" if(!defined $querylength);
	next if(!defined $blastscore);
	if($blastscore > $bestscore){
		$bestscore = $blastscore ;
		$bestid = $expr ;
	}
	if($blastscore > $maxscore){
		$foundGREAT = 1 ;
	    last ;
	}
			
}

if(!$foundGREAT && defined $bestid &&  $lastexpr ne $bestid){
	print "Did not find a great, so doing bestid $bestid with bestscore $bestscore again, lastexpr was $lastexpr\n" if($verbose);
	system("$BLASTCOMMAND -query $f1 -subject $bestid -out $outfile");
}

print "There were $N blasts to run, will stop at first $maxscore, stopped at $cnt\n";


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
