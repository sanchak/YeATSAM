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
my ($uniqfix,$mergetable,$exclude,$idx,$infile,$p1,$p2,$outfile,$cutofflen,$which_tech,$annlist);
my ($ignorefile,@expressions);
my $useEval = 0 ;
my $verbose = 0 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "mergetable=s"=>\$mergetable ,
            "uniqfix=s"=>\$uniqfix ,
            "infile=s"=>\$infile ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "annlist=s"=>\$annlist ,
            "ignorefile=s"=>\$ignorefile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "exclude=i"=>\$exclude ,
            "idx=i"=>\$idx ,
            "useEval=i"=>\$useEval ,
            "verbose=i"=>\$verbose ,
            "cutofflen=i"=>\$cutofflen ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);

usage( "Need to give a listfile -option -annlist  ") if(!defined $annlist);
usage( "Need to give a ignorefile -option -ignorefile  ") if(!defined $ignorefile);
my $CNT = 0 ; 
#my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
 my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP,$BLASTOUT,$BLASTDB) = util_SetEnvVars();

my $PWD = cwd;


# These are the MAKER identified scaffolds
my $MAKERTABLE = {};
if(defined $ignorefile){
   ($MAKERTABLE) = util_maketablefromfile($ignorefile);
}

# These are the putative mergeable scaffolds
my $MERGEDTABLE = {};
if(defined $mergetable){
   ($MERGEDTABLE) = util_maketablefromfile($mergetable);
}

# These are the putative mergeable scaffolds
my $UNIQFIX = {};
if(defined $uniqfix){
   ($UNIQFIX) = util_maketablefromfile($uniqfix);
}

## this is the TRS annotations
my @list= util_read_list_sentences($annlist);
my $anntable = {};
map { chomp;  my @l = split ; $anntable->{$l[0]} = $_ ; } @list ;

my $done = {};


my $totalnotinmaker = 0 ;
my $ofhlog = util_write($outfile);

my @values =       qw(5000 4000 3000 2000 1000 500 1 );
my @percentvals = qw(80   80    80   80   80   90 90 ); ## make things strict for smaller scaffolds

if(1){
   foreach my $val (@values){
      my $proteinpercent = shift @percentvals ;
      my $not = ProcessOneIter($infile,$outfile,$val,$done,$proteinpercent);
      $totalnotinmaker = $totalnotinmaker + $not ;
   }
}
else{
    my $not = ProcessOneIter($infile,$outfile,-1,$done,60);
   $totalnotinmaker = $totalnotinmaker + $not ;
}
print "Total NotinMaker = $totalnotinmaker\n";


sub ProcessOneIter{
    my ($infile,$outfile,$CUTOFFscafflen,$DONE,$proteinpercent) =@_ ;
    my $ofh = util_write("$outfile.$CUTOFFscafflen");
    
    my $ifh = util_read($infile);
    
    my $info = {};
    my $ignpercentmatched = 0 ;
    my $ignscafflen = 0 ;
    my $ignnonanno = 0 ;
    my $ignSmallScaff = 0 ;
    while(<$ifh>){
        chomp ;
        next if(/^\s\s*$/);
        next if(/^\t*$/);
        next if(/^\s*#/);
        my ($trs , $scafflen , $bestblastscore , $exoncount , $introncount , $diffMatchedQa , $percentmatched , $querylength , $qS , $qE , $sS , $sE , $diffQ , $diffS , $diffMatchedQ , $diffMatchedQ2S , $ignoredfar , $ignoredrev , $maxintronlength , $orderscrewed , $foundNNNinProximity , $foundNNNinIntron , $Subjectname, $isrev) = split ;

       die if(!($percentmatched =~ /%/));
       $percentmatched =~ s/%//;

	   # the rev string has two versions - mine and emboss
       $Subjectname =~ s/_rev//;
    
       next if(exists $DONE->{$Subjectname});
    

	   # if CUTOFFscafflen < 0, then do all - and do not ignore MAKER (see below)
       if($scafflen < $CUTOFFscafflen){
           $ignscafflen++;
           next ;
       }
    
       ### ignore MAKER models - 
       my $TTTT = $Subjectname;
       die "kkk $_ " if(!defined $TTTT);
       ## hack - just choose the heading, making it more conservative
       $TTTT =~ s/:.*//;
       if(exists $MAKERTABLE->{$TTTT} && $CUTOFFscafflen > 0){
           next ;
       }
    
       ## Be a little more flexible if there is a long intron, or unknown in the intron/proximity
       my $percentcutoff = $maxintronlength > 2000 || $foundNNNinIntron || $foundNNNinProximity ? 70 : 90 ;
       if($percentmatched < $percentcutoff){
           $ignpercentmatched++;
    	   next ;
       }

	   my $PERCENTM = ($diffMatchedQ2S/$querylength)*100;
       if($PERCENTM < -80){
	   	   #die "$_ kkkkkkkk";
           $ignSmallScaff++;
           next ;
       }
       if(! exists $anntable->{$trs}){
           $ignnonanno++;
    	   next ;
       }
    
       $info->{$Subjectname} = {} if(!defined $info->{$Subjectname});
       $info->{$Subjectname}->{$trs} = $_ ;
    }
    
    close($ifh);
    
    my $NNotinMaker = (keys %{$info});
    
    my $orderscrewedN = 0 ;
    my $foundNNNinProximityN = 0 ;
    my $foundNNNinIntronN = 0 ;
    my $longINTRON = 0 ;
    my $CNTFOUND = 0;
    foreach my $scaff (keys %{$info}){
         my $table = $info->{$scaff} ; 
         my $first = 0;
         my $N =     (keys %{$table});
         my $CNTFOUNDTRS = 0;
         my $scafflenCC  ;
         foreach my $trs (keys %{$table}){
             my $line = $table->{$trs};
            my ($trs , $scafflen , $bestblastscore , $exoncount , $introncount , $diffMatchedQa , $percentmatched , $querylength , $qS , $qE , $sS , $sE , $diffQ , $diffS , $diffMatchedQ , $diffMatchedQ2S , $ignoredfar , $ignoredrev , $maxintronlength , $orderscrewed , $foundNNNinProximity , $foundNNNinIntron , $Subjectname, $isrev) = split " ",$line ;
        $scafflenCC = $scafflen;
    
    	$isrev = 1 if(defined $isrev || $Subjectname =~ /_rev/);
        $Subjectname =~ s/_rev//;
    
        if($foundNNNinIntron){
           $foundNNNinIntronN++ ;
        }
        elsif($foundNNNinProximity){
             $foundNNNinProximityN++ ;
        }
        elsif($maxintronlength > 2000){
             $longINTRON++ ;
        }
        elsif($orderscrewed){
           $orderscrewedN++ ;
        }
    
        die if(! exists $anntable->{$trs});
    
        my $ann = $anntable->{$trs};
        my @vals = split " ",$ann ;
        my $EVAL = $vals[@vals -1];
        my $PROTEINPERCENT = $vals[@vals -3];
		## this looks at the protein perentage
		## if it was merged, lower the threshold
		if(exists $MERGEDTABLE->{$trs} || exists $UNIQFIX->{$trs} ){
            if($PROTEINPERCENT < 40){
			   next ;
		    }
			else{
			      print "Found $trs with percent $PROTEINPERCENT that is putative merged or uniqfix\n" if($verbose);
			}
		}
        elsif($PROTEINPERCENT < $proteinpercent){
			next ;
		}
    
        $CNTFOUNDTRS++; 
        if($CNTFOUNDTRS eq 1){
            print $ofh  "$scaff $scafflenCC\n" ;
            $CNTFOUND++;
         }


		if(exists $MERGEDTABLE->{$trs} || exists $UNIQFIX->{$trs} ){
            print $ofh "\t#####\n";
		}
        print $ofh "\t$ann \n";
        print $ofh "\t\t$line \n";
    
    	} # iterate on TRS
    	if($CNTFOUNDTRS){
    	   print $ofh "\n==========================\n"  ;
           $DONE->{$scaff} = 1 ;
    	}
    }
    
	util_printAndWrite2Script( "CUTOFFscafflen=$CUTOFFscafflen NotinMaker=$CNTFOUND ignscafflen=$ignscafflen foundNNNinProximity=$foundNNNinProximityN foundNNNinIntron=$foundNNNinIntronN longINTRON=$longINTRON orderscrewed=$orderscrewedN ignpercentmatched $ignpercentmatched,  ignnonanno =$ignnonanno ignSmallScaff=$ignSmallScaff\n",$ofhlog);
	return $CNTFOUND;
    
}


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
