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
use MyGNM;


use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($FILE,$idx,$infile,$listfile,$postfix,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$protein);
my ($ignorefile,@listfiles);
my $howmany = 100000 ;
my $verbose = 0 ;
my $proteincpercent = 60 ;

## there should be $cntoftissues > $countpercentUP - and the rest should be lower than $countpercentDOWN
my $countpercentUP = 50 ;
my $cntoftissues ;

my $countpercentDOWN = 20 ;

GetOptions(
            "FILE=s"=>\$FILE ,
            "protein=s"=>\$protein ,
            "postfix=s"=>\$postfix ,
            "infile=s"=>\$infile ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "ignorefile=s"=>\$ignorefile ,
            "outfile=s"=>\$outfile ,
            "listfiles=s"=>\@listfiles,
            "idx=i"=>\$idx ,
            "verbose=i"=>\$verbose ,
            "cntoftissues=i"=>\$cntoftissues ,
            "proteincpercent=i"=>\$proteincpercent ,
            "cutoff=f"=>\$cutoff ,
           );
my $PWD = cwd;
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
usage( "Need to give a output file name => option -cutoff ") if(!defined $cutoff);
usage( "Need to give a output file name => option -cntoftissues ") if(!defined $cntoftissues);
#usage( "Need to give a input file name => option -infile ") if(!defined $infile);
#my $ifh = util_read($infile);
#usage( "Need to give a input file name => option -fastafile ") if(!defined $fastafile);
#my ($tablefasta,$infoSeq2PDB,$mapChainedName2Name) = util_parsePDBSEQRESNEW($fastafile,0,0);

 my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP,$BLASTOUT,$BLASTDB) = util_SetEnvVars();

my $MINCOUNT = 0 ;
my $MINLENGHT = 100 ;
my $MAXAADIFF = 45 ;
my $CNT = 0 ;

my $countStr = Config_getCountStr();
# for reading fasta
#my ($fe,$x) =  util_readfasta($infile);
#chomp $x ;


my ($tableTRS) = util_maketablefromfile("list.transcriptome.clean");
$CNT = (keys %{$tableTRS});
print "There are $CNT transcripts \n";


my ($tableRRNA) = util_maketablefromfile("DATA/RRNA.anno.1000.anno.real");

my ($tableTRSlength) = util_maketablefromfile("DATA/length.trs");
my ($tableScaffLength) = util_maketablefromfile("DATA/length.scaff");
my ($tableRepeatAA) = util_maketablefromfile("DATA/repeatAA.txt");
my ($infoAnno) = GNM_parseANNFile("genome.ann.realgood",$tableRRNA);
$CNT = (keys %{$infoAnno});
print "There are $CNT annotate TRS from genome.ann.realgood \n";
my ($tableNotanno) = util_subtract2tables($tableTRS,$infoAnno);
$CNT = (keys %{$tableNotanno});
print "There are $CNT non annotated TRS \n";

## 
$FILE = "COUNTS/bwa_counts_run1.txt.0.CLEAN.normalized.usingMeanSD.splicemerge";
$FILE = "COUNTS/bwa_counts_run1.txt.0.CLEAN.normalized.usingMeanSD.splicemerge" ;  ## this is used for housekeeping
$FILE = "COUNTS/bwa_counts_run1.txt.0.CLEAN.normalized.usingMeanSD.splicemerge.homologymerge";
$FILE = "COUNTS/bwa_counts_run1.txt.0.CLEAN.normalized.usingMeanSD.splicemerge.homologymerge.500";
$FILE = "COUNTS/bwa_counts_run1.txt.0.CLEAN.normalized.usingMeanSD";
print "Doing $FILE\n";
my ($tableCountSummed,$allinfocounts) = GNM_parseCountSummed("$FILE.anno",$MINCOUNT);
$CNT = (keys %{$tableCountSummed});
print "There are $CNT trs with cnt > $MINCOUNT \n";
my @sortedCounts = sort { $tableCountSummed->{$b} <=> $tableCountSummed->{$a}} (sort keys %{$tableCountSummed});

#save the complete count as a table with the TRS
my ($tableCount) = util_maketablefromfile_firstentry("$FILE");


my $sorting = {};
my @TISSUELIST = Config_getTissueList();

system("mkdir -p UNBIAS.$cutoff");
foreach my $t (@TISSUELIST){
	my $ofhsingle = util_write("UNBIAS.$cutoff/$t.uniq");
	my $ofhsinglec = util_write("UNBIAS.$cutoff/$t.uniq.c");
	my $ofhsinglea = util_write("UNBIAS.$cutoff/$t.uniq.a");
}

my $countonetrs = {};

foreach my $nm (sort keys %{$allinfocounts}){
	my $full = $tableCount->{$nm} ;

	my $TOTAL=	     $allinfocounts->{$nm}->{TOTAL} ;
	my $NUMBER = $allinfocounts->{$nm}->{NUMBER} ;
	my $MEAN = $allinfocounts->{$nm}->{MEAN} ;
	my $SD = $allinfocounts->{$nm}->{SD} ;
	my $RATION = $allinfocounts->{$nm}->{RATION} ;

	 my ($nmxx,@l) = split " ",$full ; 
	 my $N = @l ;
	 my $max = util_get_max(\@l);
	 my $UPPERCENT = ($max * $countpercentUP)/100 ; 
	 my $DOWNPERCENT = ($max * $countpercentDOWN)/100 ; 

	 my $cnt = 0 ;
	 my $upcnt = 0 ;
	 my $downcnt= 0 ;
	 my $singleidx = -1 ;
	 # Count number of tissues which has more than 10% of the max
	 foreach my $idx (0..19){
	 	my $i = $l[$idx];
	 	if($i >= $UPPERCENT){
			$singleidx = $idx ;
			$upcnt++;
		}
		elsif ($i <= $DOWNPERCENT){
			$downcnt++;
		}
	 }


	## if the protein is annotated
    if(exists $infoAnno->{$nm}){
		my $annline =  $infoAnno->{$nm}->{FULLLINE};
		my $PERCENT = $infoAnno->{$nm}->{PERCENT};

	    if($MEAN > $cutoff && $PERCENT > $proteincpercent ){
			my $AAA = 20 - $cntoftissues;
	        if($upcnt eq $cntoftissues && $downcnt eq $AAA){
	 	       my $NM = $TISSUELIST[$singleidx];
			   $countonetrs->{$nm} = $NM ;

	        }
			
		    $sorting->{$full} = $TOTAL;
	    }
	}
}


my $ofh = util_write($outfile);
my @sorted = sort { $sorting->{$b} <=> $sorting->{$a} } (keys %{$sorting});
foreach my $full (@sorted){
	my @l = split " ",$full;
	my $nm = $l[0];


	if(exists $countonetrs->{$nm}){
		my $NM = $countonetrs->{$nm} ;
		my $annline =  $infoAnno->{$nm}->{FULLLINE};
		my $ofhsingle = util_append("UNBIAS.$cutoff/$NM.uniq");
		my $ofhsinglea = util_append("UNBIAS.$cutoff/$NM.uniq.a");
		my $ofhsinglec = util_append("UNBIAS.$cutoff/$NM.uniq.c");
		print $ofhsingle "$full \n";
		print $ofhsinglec "$full \n";
		print $ofhsingle "$annline \n";
		print $ofhsinglea "$annline \n";
	}
	print $ofh "$full \n";
}

#system("findPWDistance.pl -ig ~/blank -cutoff 0 -inf COUNTS/bwa_counts_run1.txt.0.CLEAN.normalized.usingMeanSD.anno -anno -out PW");



sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
