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
my ($genomefile,$infile,$f1,$f2,$outfile,$cutoff,$getcommands,$which_tech,$listfile,$protein);
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
            "outfile=s"=>\$outfile ,
            "outdir=s"=>\$outdir ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "isNT=i"=>\$isNT ,
            "getcommands=i"=>\$getcommands ,
            "verbose=i"=>\$verbose ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);

usage( "Need to give a input file name => option -infile ") if(!defined $infile);
usage( "Need to give a listfile -option -f1  ") if(!defined $f1);
usage( "Need to give a listfile -option -f2  ") if(!defined $f2);
usage( "Need to give a listfile -option -isNT  ") if(!defined $isNT);
usage( "Need to give a outdir pdb id -option -outdir  ") if(!defined $outdir);

## genomefile is to get the annotation
usage( "Need to give a output file name => option -genomefile ") if(!defined $genomefile && !$getcommands);
my $map2name = {};
#$map2name = util_mapID2fullStringInFastaFile($genomefile) ;

my $CNT = 0 ; 
#my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
 my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP,$BLASTOUT,$BLASTDB) = util_SetEnvVars();


my $BLASTCOMMAND = $isNT ? "blastn" : "blastp";
my $ifh = util_read($infile);


my $TARGETS = {};
while(<$ifh>){
	my ($i,@l) = split ;

	my $bestmatch = 0 ;
	my $cnt = 0 ;

	my $f2files = "";
	foreach my $j (@l){
		 next if($i eq $j); ## Can happen for self ...
		 my $outfile = "$outdir/$i.blast.nt";
	     system("$BLASTCOMMAND -query $f1/$i.ALL.1.fasta -subject $f2/$j.ALL.1.fasta -out $outfile");

        my ($info,$querylength,$Subjectname,$queryname,$blastscore,$expect) = GNM_PARSEBLAST_BESTVAL($outfile);
        die "- querylength not defined in $outfile" if(!defined $querylength);
        next if(!defined $blastscore);
        if($blastscore > $maxscore){
            last ;
        }


	}


}

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
