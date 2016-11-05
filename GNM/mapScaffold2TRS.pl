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
my ($trs,$sname,$FASTADIR,$SCAFFOLDDIR,$blastdir,$force,$justchecking,$checkforsame,$infile,$p1,$p2,$cutoff,$which_tech,$listfile,$protein);
my ($promoters,$scaffoldfasta,$ignorefile);
my $howmany = 100000 ;
my $buffersize = 1000 ;
my $verbose = 1 ;

my $percentlength = 10;
my $percentmatched = 70;
my $percentidentity = 30;


GetOptions(
            "force"=>\$force ,
            "sname=s"=>\$sname ,
            "trs=s"=>\$trs ,
            "protein=s"=>\$protein ,
            "FASTADIR=s"=>\$FASTADIR ,
            "SCAFFOLDDIR=s"=>\$SCAFFOLDDIR ,
            "infile=s"=>\$infile ,
            "blastdir=s"=>\$blastdir ,
            "scaffoldfasta=s"=>\$scaffoldfasta ,
            "checkforsame"=>\$checkforsame ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            "promoters"=>\$promoters ,
            "howmany=i"=>\$howmany ,
            "buffersize=i"=>\$buffersize ,
            "verbose=i"=>\$verbose ,
            "percentlength=i"=>\$percentlength ,
            "percentmatched=i"=>\$percentmatched ,
            "percentidentity=i"=>\$percentidentity ,
            "cutoff=f"=>\$cutoff ,
           );

my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP,$BLASTOUT,$BLASTDB) = util_SetEnvVars();
my $HARD_DIFFFORMATCH = 5;
my $HARD_LONGDISTANCEFACTOR = 10;


## hardcoded 
$blastdir = "BLASTOUT_SCAFFOLD/" if(!defined $blastdir);
my $infodir = "SCAFF.INFO/";
my $introndir = "$infodir/INTRONS";
system("mkdir -p $infodir");
system("mkdir -p $introndir");

my $ERRSTATE = 0 ;
print "Checking for BLAST OUTS \n";

my $infile = "$blastdir/$trs.$sname.blast";
if(! -e $infile){
    die " $infile does not exist ";
}
  
$scaffoldfasta = "$SCAFFOLDDIR/$sname.ALL.1.fasta" ;
my ($infoScaff,$infoSeq2PDB,$mapChainedName2Name) = util_parsePDBSEQRESNEW($scaffoldfasta,0,0);
my ($strScaffold) = $infoScaff->{$sname};


my $strScaffold = $infoScaff->{$sname};

my $infile = "$blastdir/$trs.$sname.blast";

my ($str,$extractcommand,$NAMMEE) = LOCAL_MapOneTRS2Scaffold($trs,$sname,$infile,$strScaffold,$HARD_DIFFFORMATCH,$HARD_LONGDISTANCEFACTOR);

my $outfile1 = "$infodir/$sname.$trs.info";
my $outfile2 = "$infodir/$trs.info";
my $outfile3 = "$infodir/$trs.$sname.extract.csh";
my $ofh1 = util_write($outfile1);
my $ofh2 = util_write($outfile2);
my $ofh3 = util_write($outfile3);
print $ofh1 "$str\n";
print $ofh2 "$str\n";

print $ofh3 "$extractcommand\n";





exit ;

## Keep this function here...
sub LOCAL_MapOneTRS2Scaffold{
    my ($trs,$sname,$infile,$strScaffold,$HARD_DIFFFORMATCH,$HARD_LONGDISTANCEFACTOR) = @_ ;
    my ($info,$querylength,$Subjectname,$queryname) = GNM_PARSEBLAST($infile);
	my $isrev = 0 ;
	if($Subjectname =~ /_rev/){
		$Subjectname =~ s/_rev//;

		$isrev = 1 ;
	}
	$Subjectname =~ s/ //g;
	die "$sname ne $Subjectname" if($sname ne $Subjectname);
	die "$trs ne $queryname"  if($trs ne $queryname);
    
    my 	$counter = 0;
    die "$trs - querylength not defined in $infile" if(!defined $querylength);
    
    my $done = {};
    my @others ;
    my $found100percent = 0;
    my @allQ ;
    my @allS ;
    
    ## Sort based on blastscore
                                   
    my $sort = {};
    my $DONEBSCORE = {};
    foreach my $k (@{$info}){
        my $org = $k ;
        my ($name,$description,$blastscore,$subjectlength,$iden_percent,$subjectmatched,$querystart,$queryend,$subjectstart,$subjectend,$expect) = split " ",$k ;
        while(exists $DONEBSCORE->{$blastscore}){
            $blastscore = $blastscore - 0.1;
        }
        $DONEBSCORE->{$blastscore} = 1 ;
        $sort->{$k} = $blastscore ;
    }
    my @sorted = sort {$sort->{$b} <=> $sort->{$a}} (keys %{$sort});

    
    
    ## keep a count of where we start
    my $ignoredfar = 0 ;
    my $ignoredrev = 0 ;
    my $TOTALMATCHEDINSUB = 0 ;
    
    my $MATCHEDScale = {};
    my $bestblastscore ;
    my $MAPSTARTOFBOTH = {};
    
    my $ERRSTATE = 0 ;
    my $START ;
    while(@sorted){
		my $k = shift @sorted ;
        my $org = $k ;
        my ($name,$description,$blastscore,$subjectlength,$iden_percent,$subjectmatched,$querystart,$queryend,$subjectstart,$subjectend,$expect) = split " ",$k ;
        print "$name,$description,$blastscore,$subjectlength,$iden_percent,l=$subjectmatched,$querystart,$queryend,$subjectstart,$subjectend,$expect\n" if($verbose);
        $bestblastscore = $blastscore if(!defined $bestblastscore);
		$counter++;
    
        if (!defined $START && $counter eq 1){
           if($subjectstart > $subjectend){
               # this was a rev seq problem - but sometimes a tied score can bring the plus/plus fwd
			   ## make sure there is a tie

			   $k = shift @sorted ;
			   my $xxblastscore ;
               ($name,$description,$xxblastscore,$subjectlength,$iden_percent,$subjectmatched,$querystart,$queryend,$subjectstart,$subjectend,$expect) = split " ",$k ;
			   if($xxblastscore eq $blastscore){
			   }
			   else{
                   die  "Error: kk Should not have got reverse here for $trs and $infile.";
			   }
           }
           $START = $subjectstart ;
        }
        else{
           if($subjectstart > $subjectend){
               print "Info: ignoring since it is reversed, and not the first entry\n" if($verbose);
               $ignoredrev++;
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
        $MAPSTARTOFBOTH->{$querystart} = $subjectstart ;
    
        # Annotate a scale
        foreach my $i ($subjectstart..$subjectend){
            $MATCHEDScale->{$i} = 1 ;
        }
    
    
        $TOTALMATCHEDINSUB = $TOTALMATCHEDINSUB + $subjectmatched ;
        ## if we find a match which is long, why bother further
        if(abs($diffQ -$querylength) < $HARD_DIFFFORMATCH){
            print "Info: $diffQ -$querylength, equal or almost equal match found \n" if($verbose);
            last ;
        }
    }
    if($ERRSTATE){
        die "Killing due to previous errors";
    }
    
    
    print "$querylength is querylength and TOTALMATCHEDINSUB is $TOTALMATCHEDINSUB \n" if($verbose);
    
    
    my @sortQ = sort {$a <=> $b} @allQ ;
    my @sortS = sort {$a <=> $b} @allS ;
    
    my $N = @sortQ ;
    die if ($N ne @sortS);
    
    my $qS = $sortQ[0];
    my $qE = $sortQ[$N-1];
    
    my $sS = $sortS[0];
    my $sE = $sortS[$N-1];
    
    
    ## ensure the starts are the same 
    die if(!exists $MAPSTARTOFBOTH->{$qS});
    my $MAPPEDEND = $MAPSTARTOFBOTH->{$qS};
    my $orderscrewed = 0;
    if($MAPPEDEND ne $sS){
        print "Order screwed up \n" if($verbose);
        $orderscrewed = 1;
    }
    
    my $foundNNNinProximity= 0 ;
                                                            
    {
       my @llll = split "",$strScaffold;
       my $NNN = @llll  ;
       my $postend = $sE+500 ;
       $postend = $NNN  if($postend >= $NNN);
       my $poststr = util_extractSliceFromFastaString($strScaffold,$sE,$postend);
       $foundNNNinProximity= 1 if($poststr =~ /NNNNNNNN/);
    
       my $prestart = $sS-500 ;
       $prestart = 0 if($prestart < 0);
       $poststr = util_extractSliceFromFastaString($strScaffold,$prestart,$sS);
       $foundNNNinProximity= 1 if($poststr =~ /NNNNNNNNN/);
    }
    
    my $diffQ = $qE - $qS + 1 ;
    my $diffS = $sE - $sS + 1  ;
    
    my $diffMatchedQ = $querylength - $diffQ ;

    ## this can be misleading, esp when orderscrewed is true, but is necesssary when there are X's...
    my $percentmatched = int(100 *($diffQ/$querylength));
    ## this is the real thing matched
    my $percentmatchedreal = int(100 *($TOTALMATCHEDINSUB/$querylength));
    if($diffMatchedQ < $HARD_DIFFFORMATCH){
        print "Exact match (diff=$diffMatchedQ) for $trs and $Subjectname\n";
    }
    
    my $diffMatchedQ2S = $diffS - $querylength ;
    
    
    my $introncount = 0 ;
    my $exoncount = 1 ;

    my $maxintronlength = 0 ;
    my $foundNNNinIntron = 0 ;
    if(1)
    {
       unlink ("$introndir/$trs.intron");
       die if(!exists $MATCHEDScale->{$sS});
       die if(!exists $MATCHEDScale->{$sE});
       my $STATE = 0 ;
       my $intronstart = 0;
       foreach my $i ($sS..$sE){
           if($STATE eq 0 && !exists $MATCHEDScale->{$i}){
                  # switch to intron state 
                  $STATE  = 1;
              $exoncount++;
               $intronstart = $i ;
           }
           else{
                if($STATE eq 1 && exists $MATCHEDScale->{$i}){
                    # switch to exon state 
                    $STATE = 0 ;
                    my $intronend = $i - 1;
                    my $diff = abs($intronstart - $intronend);
                    $maxintronlength = $diff if($diff > $maxintronlength);
                    my $intronstr = util_extractSliceFromFastaString($strScaffold,$intronstart,$intronend);
                    $foundNNNinIntron = 1 if($intronstr =~ /NNNNNNNNN/);
                    my $ofhintron = util_open_or_append("$introndir/$trs.intron");
                    print $ofhintron ">$sname.$intronstart.$intronend\n";
                    print $ofhintron ">$intronstr\n";
                    $introncount++;
                }
           }
       }
    }

    my $str = "$trs $Subjectname %=$percentmatched Ex=$exoncount In=$introncount qL=$querylength qs=$qS qE=$qE sS=$sS sE=$sE diffQ=$diffQ diffS=$diffS diffMatchedQ=$diffMatchedQ diffMatchedQ2S=$diffMatchedQ2S $ignoredfar $ignoredrev maxintronlength=$maxintronlength orderscrewed=$orderscrewed NNinProximity=$foundNNNinProximity NNinIntron=$foundNNNinIntron rev=$isrev";

	my $FILEFASTA = "$Subjectname.ALL.1.fasta.comp.fasta";
	if(!$isrev){
	     $FILEFASTA = "$Subjectname.ALL.1.fasta";	
	}

	my $extractStart = $sS - $buffersize;
	my $extractEnd = $sE + $buffersize;
	my $NAMMEE = "$trs.$Subjectname.$sS.$sE.EXTRACT";
	my $extractcommand = "newfile.csh $infodir/$NAMMEE.ALL.1.fasta \n";
	$extractcommand = $extractcommand . "extractslicefromfasta.pl -in SCAFFOLDDIR/$FILEFASTA -outf $infodir/$NAMMEE.ALL.1.fasta -star $extractStart -end $extractEnd -stric 0 \n\n";
	$extractcommand = $extractcommand . " getorf $infodir/$NAMMEE.ALL.1.fasta ORFTMP/$NAMMEE.orf  \n 
	fixORFnames.pl -inf ORFTMP/$NAMMEE.orf -out ORF/$NAMMEE.orf \n 
	findrepeatedorfs.pl -trs $NAMMEE -orfdir ORF/ -write 4 \n";
	$extractcommand = $extractcommand . "orfreallength.pl -in FASTADIR_JUSTONLONGEST/$NAMMEE.ALL.1.fasta -out FASTADIR_JUSTONLONGEST/$NAMMEE.LONGESTM.ALL.1.fasta \n ";
	$extractcommand = $extractcommand . "yeatsCoExLTPPattern.pl -outf TMP/$NAMMEE -fas FASTADIR_JUSTONLONGEST/$NAMMEE.LONGESTM.ALL.1.fasta \n";
    return ($str,$extractcommand,$NAMMEE) ;
}


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}

