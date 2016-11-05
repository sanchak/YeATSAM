#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use MyGeom;
use MyGNM;
use MyConfigs;
use PDB;
use ConfigPDB;
use Math::Geometry ;
use Math::Geometry::Planar;
use Memory::Usage;
use Algorithm::Combinatorics qw(combinations) ;
my $mu = Memory::Usage->new();
$mu->record('');
use AAConfig;



use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($fastafile,$idx,$infile,$postfix,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$protein);
my ($ignorefile,@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "fastafile=s"=>\$fastafile ,
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
die "Dont recognize command line arg @ARGV " if(!@ARGV);

my ($targetID,$fdir1,$fdir2,@blastouts) = @ARGV ;

my ($b1,$b2,$b3) = @blastouts ;
$b3 = "" if(!defined $b3);
#my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
 my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP,$BLASTOUT,$BLASTDB) = util_SetEnvVars();
my $aaconfig = new AAConfig("$SRC/aa.config");
my  ($tableThree2One,$tableOne2Three,@sortedSingle) = Config_AACodes();

system("mkdir -p MERGEUNIQ");

my $trs1 = $b1 ;
my $trs2 = $b2 ;
my $trs3 = $b3 ;
$trs1 =~ s/.blast.nt//;
$trs2 =~ s/.blast.nt//;
$trs3 =~ s/.blast.nt//;

$trs1 =~ s/BLASTOUT_AA.//;
$trs2 =~ s/BLASTOUT_AA.//;
$trs3 =~ s/BLASTOUT_AA.//;


my $TRSREALNAME = $trs1 ;
$TRSREALNAME =~ s/.ORF.*//;


my $f1 = "$fdir1/$trs1.ALL.1.fasta" ;
my $f2 = "$fdir1/$trs2.ALL.1.fasta" ;
my $f3 = "$fdir1/$trs3.ALL.1.fasta" ;
my $ftarge = "$fdir2/$targetID.ALL.1.fasta" ;

die "one of these not there .. $f1 .. $f2 .. $ftarge" if(! -e $f1 || ! -e $f2 || ! -e $ftarge);
die "not there $f3" if($b3 ne "" && ! -e $f3);

## Read original blast scores 
my (@l1) = GNM_PARSEBLAST_BESTVAL($b1);
my (@l2) = GNM_PARSEBLAST_BESTVAL($b2);
my @l3 ;
@l3 = GNM_PARSEBLAST_BESTVAL($b3) if($b3 ne "");


my $NN = @l1 ;
my $score1 = $l1[$NN -2];
my $score2 = $l2[$NN -2];
my $score3 = 0 ;
$score3 = $l3[$NN -2] if($b3 ne "");

my $SUMALLTHREE = $score1+$score2+$score3 ;
my $SUMeoightpercent = 0.8 * $SUMALLTHREE ;
my $SUMninentypercent = 0.9 * $SUMALLTHREE ;
print "scores are $score1, $score2 and $score3\n" if($verbose);

if ($score1 < 100 || $score2 < 100){
	exit ;
}

my ($str1,$firstline1) = util_readfasta($f1);
my ($str2,$firstline2) = util_readfasta($f2);
my ($str3) = "" ;
($str3) = util_readfasta($f3) if($b3 ne "");

my @SEQ ;
my $A = $str1 . "ZZZ" . $str2 . "ZZZ" . $str3 ; push @SEQ, $A;
my $B = $str2 . "ZZZ" . $str1 . "ZZZ" . $str3 ; push @SEQ, $B;
my $C = $str1 . "ZZZ" . $str3 . "ZZZ" . $str2 ; push @SEQ, $C;
my $D = $str2 . "ZZZ" . $str3 . "ZZZ" . $str1 ; push @SEQ, $D;
my $E = $str3 . "ZZZ" . $str2 . "ZZZ" . $str1 ; push @SEQ, $E;
my $F = $str3 . "ZZZ" . $str1 . "ZZZ" . $str2 ; push @SEQ, $F;

my $ONE = "$trs1.REM.$trs2.REM.$trs3";
my $TWO = "$trs2.REM.$trs1.REM.$trs3";
my $THREE = "$trs1.REM.$trs3.REM.$trs2";
my $FOUR = "$trs2.REM.$trs3.REM.$trs1";
my $FIVE = "$trs3.REM.$trs2.REM.$trs1";
my $SIX = "$trs3.REM.$trs1.REM.$trs2";

my @NAMES ;
push @NAMES, $ONE ;
push @NAMES, $TWO ;
push @NAMES, $THREE ;
push @NAMES, $FOUR ;
push @NAMES, $FIVE ;
push @NAMES, $SIX ;



my @SCORES;
foreach my $idx (0..5){
	my $NAME = $NAMES[$idx];
	my $SEQ = $SEQ[$idx];
    my $tmp1 = "MERGEUNIQ/$NAME.ALL.1.fasta";
    my $ofh1 = util_write($tmp1);
    print $ofh1 ">$NAME\n";
    print $ofh1 "$SEQ\n";
    close($ofh1);
	if(! -e "MERGEUNIQ/$NAME.blast.nt"){
        system("myblastcompare2Fastafiles.csh $tmp1 $ftarge MERGEUNIQ/$NAME.blast.nt P");
	}
    my (@l3) = GNM_PARSEBLAST_BESTVAL("MERGEUNIQ/$NAME.blast.nt");
    my $score3 = $l3[$NN -2];
	push @SCORES, $score3;
}


my $REALNAME = $trs1 ;
$REALNAME =~ s/.ORF.*//;
my $ofhout = util_write("MERGEUNIQ//$REALNAME.out");

### Only two 
if($b3 eq ""){

   my $score3 = $SCORES[0];
   my $score4 = $SCORES[1];
   
   if($score3 < $SUMeoightpercent  &&  $score4 < $SUMeoightpercent){
   	  print $ofhout "$score3 < $SUMeoightpercent  &&  $score4 < $SUMeoightpercent \n";
   	  exit ;
   }
   
   
   
   my $maxofboth = $score1 > $score2 ? $score1 : $score2 ;
   if($score3 <= $maxofboth && $score4 <= $maxofboth ){
             print $ofhout "$ONE BOTH SIDES? maxofboth = $maxofboth LESS ======== scores are $score1 $score2  later $score3 $score4 \n";
   	exit ;
   }
   if($score3 > $SUMeoightpercent  &&  $score4 > $SUMeoightpercent){
        if($score3 > $SUMninentypercent  &&  $score4 > $SUMninentypercent){
             print $ofhout "$ONE BOTH SIDES? ======== scores are $score1 $score2  later $score3 $score4 \n";
   		  exit ;
        }
   	if($score3 > $SUMninentypercent){
   		print $ofhout "$ONE GREAT $score1 $score2  $score3 = 3 \n";
   	}
   	else{
   		print $ofhout "$TWO GREAT $score1 $score2 $score4 = 4 \n";
   		}
   }
   else{
   	if($score3 > $SUMeoightpercent){
   		print $ofhout "$ONE GREAT $score1 $score2  $score3 = 3 \n";
   	}
   	else{
   		print $ofhout "$TWO GREAT $score1 $score2 $score4 = 4 \n";
   	}
   }

}
else{
	my $maxscore = 0 ;
	my $maxidx = 0 ;
    foreach my $idx (0..5){
        my $s = $SCORES[$idx];
		if($s > $maxscore){
			$maxscore = $s ;
			$maxidx = $idx ;
		}
    }
	my $NAME = $NAMES[$maxidx];
	my $percent = 0.3 * $score1;
	my $thresh = $score1 + $percent ;
	print "thresh=$thresh, percent = 0.3 * $score1 = $percent \n";
	if($maxscore > $SUMninentypercent){
	    print "$NAME maxscore = $maxscore, maxidx = $maxidx $score1 $score2 $score3 SUMALLTHREE=$SUMALLTHREE, SUMninentypercent = $SUMninentypercent\n";
   		print $ofhout "#$NAME GREAT maxscore = $maxscore, maxidx = $maxidx last scores $score1 $score2 $score3\n";
   		print $ofhout "\\cp -r MERGEUNIQ/$NAME.ALL.1.fasta FASTADIR_ORFBEST/$TRSREALNAME.MERGESAME.ALL.1.fasta \n";
		
		my $ofhremove = util_open_or_append("MERGEUNIQ/remove");
		my $ofhadd = util_open_or_append("MERGEUNIQ/add");
		print $ofhadd "$TRSREALNAME.MERGESAME";
		print $ofhremove "$TRSREALNAME";
	}
	else{
		print $ofhout "#$NAME Nothing doing: since $maxscore < $thresh\n";
	}
}

