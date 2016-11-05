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
my $diffCCandCdotC = 25 ;
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
            "diffCCandCdotC=i"=>\$diffCCandCdotC ,
            "verbose=i"=>\$verbose ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
usage( "Need to give a input file name => option -fastafile ") if(!defined $fastafile);
print "Reading fastafile $fastafile\n";
my ($tablefasta,$infoSeq2PDB,$mapChainedName2Name) = util_parsePDBSEQRESNEW($fastafile,0,0);


my $tableallowed = {};
$tableallowed->{8} = 1 ;
$tableallowed->{9} = 1 ;
$tableallowed->{11} = 1 ;
$tableallowed->{12} = 1 ;
$tableallowed->{13} = 1 ;
$tableallowed->{14} = 1 ;
$tableallowed->{19} = 1 ;
$tableallowed->{20} = 1 ;

system("rm -f $outfile.*");

my $ofhCC = util_write("$outfile.CC");
my $ofhCCfa = util_write("$outfile.CC.fa");
my $ofhprobinmid = util_write("$outfile.probinmid");
my $ofhPattern = util_write("$outfile.pattern");
my $ofhdebug = util_write("$outfile.debug");
#my $PATTERN = "CC.G.(K|R).*(D|E)(R|K)...C.C.(R|K).*CG";
my $PATTERN = "CC.G.(K|R).*(D|E)(R|K)...C.C";
#$PATTERN = "CC.G.(K|R).*C.C";
print "Checking for pattern =~ /$PATTERN/), diffCCandCdotC=$diffCCandCdotC \n";


sub printSeqInfo{
	my ($i,$seq,$id,$N) = @_ ;
	my $ofh1 = util_open_or_append("$outfile.$id");
	my $ofh2 = util_open_or_append("$outfile.fa.$id");
	if(defined $N){
	    print $ofh1 "$i $N\n";
	}
	else{
	    print $ofh1 "$i\n";
	}
	print $ofh2 ">$i \n";
	print $ofh2 "$seq \n";
}

my $ignore = {};
my $igncnt = 0 ;
foreach my $i (sort keys %{$tablefasta}){
	my $seq = $tablefasta->{$i} ;
	my $ORIGSEQ = $seq ;
	my $ORIGLEN = length($ORIGSEQ);
	my @C = ($seq =~ /C/g);
	my @CC = ($seq =~ /CC/g);
	my @CCC = ($seq =~ /CCC/g);
	my @CdotC = ($seq =~ /C.C/g);

	if(@CCC){
		printSeqInfo ($i,$ORIGSEQ,"ign.0.CCCexists"); $igncnt++;
		$ignore->{$i} = 1 ; 
		next ;
	}
	if(!@C){
		printSeqInfo ($i,$ORIGSEQ,"ign.1.noCatall"); $igncnt++;
		$ignore->{$i} = 1 ; 
		next ;
	}
	if(!@CC){
		printSeqInfo ($i,$ORIGSEQ,"ign.2.noCC"); $igncnt++;
		$ignore->{$i} = 1 ; 
		next ;
	}
	if(@C < 8){
		my $N = @C;
		printSeqInfo ($i,$ORIGSEQ,"ign.3.Clt8",$N); $igncnt++;
		$ignore->{$i} = 1 ; 
		next ;
	}

	if(@CC eq 2){
		$seq =~ s/CC/XC/;
		if($i eq "aaTC20225.SCAFFOLD11286.EXTRACT" ){
		     $seq =~ s/CC/XC/;
		     $seq =~ s/XC/CC/;
			 print "HACK for $i\n";
			 print "seq=$seq\n";
		}
		printSeqInfo ($i,$ORIGSEQ,"warn.1.CCnumeq2"); $igncnt++;
	}

	my $patternCCandCDotC = "CC.*C.C";
	if($seq !~ /$patternCCandCDotC/){
		  printSeqInfo ($i,$ORIGSEQ,"ign.4.noCCandCCdotC"); $igncnt++;
		  $ignore->{$i} = 1 ; 
		  next ;
	}
	if(@CdotC > 2 || @CdotC eq 0){
		printSeqInfo ($i,$ORIGSEQ,"ign.5.CdotCnotoneOrTwo"); $igncnt++;
		$ignore->{$i} = 1 ; 
		next ;
	}

	my $patternCCandCDotCfor2 = "CC.*C.C.*C.C";
	my ($space1,$space2,$space3,$space4,$space5) ;
	my ($lenspace1,$lenspace2,$lenspace3,$lenspace4,$lenspace5) ;
    if($seq =~ /$patternCCandCDotCfor2/){
		my $tmpseq = $seq ;
		$tmpseq =~ s/CC/XXX/;
		$tmpseq =~ s/C.C.*C.C/C.C/;
		$tmpseq =~ s/XXX/CC/;
	    my @patt = ($tmpseq =~ /$patternCCandCDotC/g);
	    die "$i $tmpseq " if(@patt ne 1);
		my $xx = $patt[0];
		$space3 = $patt[0];
		$xx =~ s/CC//;
		$xx =~ s/C.C//;
	    my $lenpatternCCandCDotC = length($xx);
	    if($lenpatternCCandCDotC < 8 || $lenpatternCCandCDotC > $diffCCandCdotC){
		      printSeqInfo ($i,$ORIGSEQ,"ign.6.CCandCCdotCLenNotRight",$lenpatternCCandCDotC); $igncnt++;
		      print $ofhdebug "$i $lenpatternCCandCDotC = lenpatternCCandCDotC $patt[0]\n";
		      $ignore->{$i} = 1 ; 
			  next ;
	    }
	}
	else{
	    my @patt = ($seq =~ /$patternCCandCDotC/g);
	    die if(@patt ne 1);
		$space3 = $patt[0];
	    my $lenpatternCCandCDotC = length($patt[0]);
	    if($lenpatternCCandCDotC < 8 || $lenpatternCCandCDotC > $diffCCandCdotC){
		      printSeqInfo ($i,$ORIGSEQ,"ign.6.CCandCCdotCLenNotRight",$lenpatternCCandCDotC); $igncnt++;
		      print $ofhdebug "$i $lenpatternCCandCDotC = lenpatternCCandCDotC\n";
		      $ignore->{$i} = 1 ; 
			  next ;
	    }
	}

	if(1){
	   #$seq = $ORIGSEQ;
	   $seq =~ s/CC/ /;
	   my ($A,$B)  = split " ", $seq ;

	   if(!defined $B){
		  printSeqInfo ($i,$ORIGSEQ,"ign.7.CCatEnd",$seq); $igncnt++;
		  next ;
	   }

	   my @CinA = ($A =~ /C/g);
	   if(@CinA eq 1){
		  printSeqInfo ($i,$ORIGSEQ,"ign.7.CCOnlyOne",$seq); 
		  next ;
	   }

	   $A =~ s/C/ /g;;
	   my (@listA)  = split " ", $A ;
	   my $NNN = @listA ;

	   $space1 = $listA[$NNN-2];
	   $space2 = $listA[$NNN-1];
	   $lenspace1 = length($listA[$NNN-2]);
	   $lenspace2 = length($listA[$NNN-1]);
	   #print "$NNN lllllll $lenspace1, $lenspace2 \n";
	   die if(!@listA);

	   ## Divide B into two strings
	   $B =~ s/C.C/ /;
	   my (@listB)  = split " ", $B ;
	   

	   ## in between CC and C.C
	   $space3 = shift @listB;
	   $lenspace3 = length($space3);

	   if(@listB){
	       my $finalstr = shift @listB;
	       $finalstr =~ s/C/ /g;;
	       my (@listfinalstr)  = split " ", $finalstr ;
	       $lenspace4 = length($listfinalstr[0]);
	       $lenspace5 = length($listfinalstr[1]);
	   }
	   else{
	       $lenspace4 = 0 ;
	       $lenspace5 = 0 ;
       }
	   


	   my $anded = $lenspace1 && $lenspace2 && $lenspace3 && $lenspace4 && $lenspace5;
	   if(!defined $anded){
		  printSeqInfo ($i,$ORIGSEQ,"ign.9.allspacesnotdefined"); $igncnt++;
		  next ;
	   }



	   #if(!InRange($lenspace1,6,14) || !InRange($lenspace2,8,19) || !InRange($lenspace3,8,25) || !InRange($lenspace4,12,34) || !InRange($lenspace5,5,14)){
	   if(!InRange($lenspace1,6,16) || !InRange($lenspace2,8,19) || !InRange($lenspace3,8,25) || !InRange($lenspace4,12,34) || !InRange($lenspace5,3,20)){
		   printSeqInfo ($i,$ORIGSEQ,"ign.10.CCbutnotRight","$lenspace1 && $lenspace2 && $lenspace3 && $lenspace4 && $lenspace5 len= $ORIGLEN "); 
	   }
	   else{
		   print $ofhCCfa ">$i\n";
		   print $ofhCCfa "$ORIGSEQ\n";
	       print $ofhCC "$i $lenspace1 && $lenspace2 && $lenspace3 && $lenspace4 && $lenspace5 len= $ORIGLEN \n";
	   }


	 }
}

sub InRange{
	my ($val,$low,$high) = @_ ;
	my $retval = 1 ;
	if($val < $low  || $val > $high){
		$retval = 0 ;
	}
	#print "$retval === retval \n";
	return $retval ;
}

system("sortOnLast.pl -in $outfile.CC");
system("wc -l $outfile.*");

#my ($N,$first,$last,$range,$mean,$sd) = util_GetStatOfList(@l);
#my ($table) = util_make_list(\@list);
#my ($tablemerge,$newN) = util_table_merge($t1,$t2);
#my ($table) = util_mapFullLinetofirst($fname);
#my ($table,$N) = util_maketablefromfile($fname);
#my ($common,$inAbutnotinB,$inBbutnotinA) = util_table_diff($t1,$t2);
#my $junk = util_writelist2file($fname,@list);


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
