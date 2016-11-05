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
use MyGNM;


use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($justchecking,$checkforsame,$infile,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$protein);
my ($promoters,$scaffoldfasta,$ignorefile,@expressions,$trs);
my $howmany = 100000 ;
my $verbose = 0 ;

my $percentlength = 10;
my $percentmatched = 70;
my $percentidentity = 30;

GetOptions(
            "justchecking"=>\$justchecking ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "trs=s"=>\$trs ,
            "scaffoldfasta=s"=>\$scaffoldfasta ,
            "checkforsame"=>\$checkforsame ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            "outfile=s"=>\$outfile ,
            "promoters"=>\$promoters ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "verbose=i"=>\$verbose ,
            "percentlength=i"=>\$percentlength ,
            "percentmatched=i"=>\$percentmatched ,
            "percentidentity=i"=>\$percentidentity ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);

my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP,$BLASTOUT,$BLASTDB) = util_SetEnvVars();
my $HARD_LONGDISTANCEFACTOR = 10;
my $HARD_DIFFFORMATCH = 5;

print "-percentidentity $percentidentity -percentmatched $percentmatched -percentlength $percentlength\n" if($verbose);

usage( "Need to give a input file name => option -infile ") if(!defined $infile);
usage( "Need to give a input file name => option -trs ") if(!defined $trs);


#my ($info,$querylength) = util_PARSEBLAST($infile);
my ($info,$querylength,$Subjectname,$queryname) = GNM_PARSEBLAST($infile);

die "$trs - querylength not defined in $infile" if(!defined $querylength);

### Just to get the Subjectname
my $ifh = util_read($infile);
my $Subjectname ; 
while(<$ifh>){
	if(/Subject=/){
		chomp;
		s/Subject=//;
		$Subjectname = $_ ;
		last ;
	}
}


my $done = {};
my @others ;
my $found100percent = 0;
my @allQ ;
my @allS ;


## Sort based on blastscore
my $sort = {};
foreach my $k (@{$info}){
	my $org = $k ;
    my ($name,$description,$blastscore,$subjectlength,$iden_percent,$subjectmatched,$querystart,$queryend,$subjectstart,$subjectend,$expect) = split " ",$k ;
	$sort->{$k} = $blastscore ;
}
my @sorted = sort {$sort->{$b} <=> $sort->{$a}} (keys %{$sort});


## keep a count of where we start
my $START ;
my $ignoredfar = 0 ;
my $TOTALMATCHEDINSUB = 0 ;

my $MATCHEDScale = {};
foreach my $k (@sorted){
	my $org = $k ;
    my ($name,$description,$blastscore,$subjectlength,$iden_percent,$subjectmatched,$querystart,$queryend,$subjectstart,$subjectend,$expect) = split " ",$k ;
	print "$name,$description,$blastscore,$subjectlength,$iden_percent,l=$subjectmatched,$querystart,$queryend,$subjectstart,$subjectend,$expect\n" if($verbose);
	$TOTALMATCHEDINSUB = $TOTALMATCHEDINSUB + $subjectmatched ;

	if (!defined $START){
	   if($subjectstart > $subjectend){
	   	   system("touch $trs.isrev");
		   die "Error: The first is reversed\n";
	   }
	   $START = $subjectstart ;
	}
	else{
	   if($subjectstart > $subjectend){
		   print "Info: ignoring since it is reversed, and not the first entry\n" if($verbose);
		   next ;
	   }
	}

	### Ignore something which is really far
	if(abs($START - $subjectstart) > $HARD_LONGDISTANCEFACTOR*$querylength){
		print "Info: Ignored since far away \n" if($verbose);
		$ignoredfar++;
		next ;
	}


	push @allQ, $querystart;
	push @allQ, $queryend;
    my $diffQ = $queryend - $querystart + 1 ;

	push @allS, $subjectstart;
	push @allS, $subjectend;

	# Annotate a scale
	foreach my $i ($subjectstart..$subjectend){
		$MATCHEDScale->{$i} = 1 ;
	}


	## if we find a match which is long, why bother further
	if(abs($diffQ -$querylength) < $HARD_DIFFFORMATCH){
		print "Info: $diffQ -$querylength, equal or almost equal match found \n";
		last ;
	}
}




print "$querylength is querylength and TOTALMATCHEDINSUB is $TOTALMATCHEDINSUB \n";

## this was just meant for writing $trs.isrev
## please note to check this
if(defined $justchecking){
	exit ;
}

my @sortQ = sort {$a <=> $b} @allQ ;
my @sortS = sort {$a <=> $b} @allS ;

my $N = @sortQ ;
die if ($N ne @sortS);

my $qS = $sortQ[0];
my $qE = $sortQ[$N-1];

my $sS = $sortS[0];
my $sE = $sortS[$N-1];

my $diffQ = $qE - $qS + 1 ;
my $diffS = $sE - $sS + 1  ;

my $diffMatchedQ = $querylength - $diffQ ;
my $percentmatched = int(100 *($TOTALMATCHEDINSUB/$querylength));
if($diffMatchedQ < $HARD_DIFFFORMATCH){
	print "Exact match (diff=$diffMatchedQ) for $trs and $Subjectname\n";
}

my $diffMatchedQ2S = $diffS - $querylength ;


my $introncount = 0 ;
{
   unlink ("INTRONS/$trs.intron");
   die if(!exists $MATCHEDScale->{$sS});
   die if(!exists $MATCHEDScale->{$sE});
   my $STATE = 0 ;
   my $intronstart = 0;
   foreach my $i ($sS..$sE){
	   if($STATE eq 0 && !exists $MATCHEDScale->{$i}){
	   	   # switch to intron state 
	   	   $STATE  = 1;
		   $intronstart = $i ;
	   }
	   else{
	        if($STATE eq 1 && exists $MATCHEDScale->{$i}){
	   	        # switch to exon state 
				$STATE = 0 ;
				my $intronend = $i - 1;
	            system ("extractslicefromfasta.pl -out INTRONS/$trs.intron -inf $scaffoldfasta -sta $intronstart -end $intronend");
				$introncount++;
		    }
	   }
   }
}
print $ofh "$trs $Subjectname $introncount $diffMatchedQ $percentmatched $querylength $qS $qE $sS $sE $diffQ $diffS $diffMatchedQ $diffMatchedQ2S $ignoredfar \n";


if(defined $promoters){
print "========== Getting promoters ============\n";
if(defined $scaffoldfasta){
    my $strS = "";
    my $strQ = "";
	my $trsfasta = "$FASTADIR/$trs.ALL.1.fasta";
    ($strQ) = util_readfasta($trsfasta);
    my $lenQ = length($strQ);
    ($strS) = util_readfasta($scaffoldfasta);
    my $lenS = length($strS);

	#system ("extractslicefromfasta.pl -out $outfile.query -inf $trsfasta -sta $qS -end $qE");
	_ExtractSlice(500,$lenS);
	_ExtractSlice(1000,$lenS);
	_ExtractSlice(2000,$lenS);

}
}

sub _ExtractSlice{
	my ($diff,$lenScaffold) = @_ ;
	my $preStart = $sS - $diff ;
	my $preEnd = $sS - 1 ;
	system("mkdir -p PROMOTERS");
	if($preStart >= 0){
	     system ("extractslicefromfasta.pl -out PROMOTERS/$trs.$diff.5tick.fasta -inf $scaffoldfasta -sta $preStart -end $preEnd");
	}

	my $postStart = $sE + $diff ;
	my $postEnd = $sE + 1 ;
	# note the reversal of start and end
	if($postStart < $lenScaffold){
	system ("extractslicefromfasta.pl -out PROMOTERS/$trs.$diff.3tick.fasta -inf $scaffoldfasta -sta $postEnd -end $postStart");
	}
}

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}

