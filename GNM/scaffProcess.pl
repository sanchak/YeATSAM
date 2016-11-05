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
my ($fastafile,$idx,$infile,$postfix,$p1,$p2,$outfile,$cutoff,$which_tech,$groupfile,$protein);
my ($fullgenelist,@expressions);
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
            "groupfile=s"=>\$groupfile ,
            "fullgenelist=s"=>\$fullgenelist ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "idx=i"=>\$idx ,
            "verbose=i"=>\$verbose ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
my $ifh = util_read($infile);
usage( "Need to give a input file name => option -fullgenelist ") if(!defined $fullgenelist);
#my ($tablefasta,$infoSeq2PDB,$mapChainedName2Name) = util_parsePDBSEQRESNEW($fastafile,0,0);
#my ($str,$firstline) = util_readfasta($infile);


usage( "Need to give a groupfile -option -groupfile  ") if(!defined $groupfile);
#usage( "Need to give a protein pdb id -option -protein  ") if(!defined $protein);
my $CNT = 0 ; 
#my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
 my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP,$BLASTOUT,$BLASTDB) = util_SetEnvVars();
my $aaconfig = new AAConfig("$SRC/aa.config");
my  ($tableThree2One,$tableOne2Three,@sortedSingle) = Config_AACodes();


my @fullgenelist= util_read_list_sentences($fullgenelist);
my $fullgenetable = {};
map { s/\s*//g ; $fullgenetable->{$_} = 1 ; } @fullgenelist ;


my $infoRegion = {};
my $infoScaff = {};
my $infoDiff = {};
my $infoExon = {};
my $infoNNN = {};

my $scaff2trs = {};
while(<$ifh>){
     next if(/^\s*$/);
	 next if(/^\s*#/);
	 s/sS=//;
	 s/sE=//;
	 s/Ex=//;
	 s/NNinProximity=//;
	 s/NNinIntron=//;

	 my @l = split ;
	 my $trs = $l[0];
	 my $scaff = $l[1];
	 $scaff2trs->{$scaff} = {} if(! exists $scaff2trs->{$scaff}) ;
	 $scaff2trs->{$scaff}->{$trs} = 1 ;
	 my $info->{$trs} = {};
	 my $nExon = $l[3];
	 my $sS = $l[8];
	 my $sE = $l[9];
	 my $NNinProximity = $l[18];
	 my $NNinIntron = $l[19];
	 
	 $infoRegion->{$trs} = "$sS $sE";
	 $infoScaff->{$trs} = $scaff ;
	 $infoDiff->{$trs} = $sE - $sS ;
	 $infoExon->{$trs} = $nExon ;
	 $infoNNN->{$trs} = $NNinProximity * $NNinIntron ;
	 #print "$trs  $scaff $sS $sE\n";
}

my $ifhgroupfile = util_read($groupfile);
my $ALLGROUPEDTRS = {};
while(<$ifhgroupfile>){
	s/.*=N//;
	my @l = split ;
	my @newl ;
	my $ignore = 0 ;
	foreach my $i (@l){
		$i  = FixName($i);
		$ALLGROUPEDTRS->{$i} = 1 ;
		$ignore = 1  if(! exists $infoScaff->{$i});
		push @newl, $i ;
	}

	## one TRS needs to be ignored - see later TODO
	if(!$ignore){
		ProcessOneGroup($infoScaff,$infoRegion,@newl);
	}
}



## Now process the TRS not grouped
foreach my $trs (@fullgenelist){
	$trs = FixName($trs);
	next if(exists $ALLGROUPEDTRS->{$trs});
	my $scaff = $infoScaff->{$trs} ;
	if(!defined $scaff){
	   warn "$trs not there " ;
	   next ;
	}
	
	my $scaffwithtrs = $scaff2trs->{$scaff};

	my @trsL = (sort keys %{$scaffwithtrs});

	while(@trsL){
		my $trsX = shift @trsL;


		next if($trsX eq $trs);
		my $val = IsTRSMappedtoSameScaffold($trs,$trsX,$infoScaff,$infoRegion);
		if($val eq "samearea"){

		    my $X = $trs ;
		    my $Y = $trsX ;
		    $X =~ s/_.*//;
		    $Y =~ s/_.*//;
		    next if($X eq $Y);

			my $stra = $infoRegion->{$trs};
			my $strb = $infoRegion->{$trsX};
			my $diffA = $infoDiff->{$trs};
			my $diffB = $infoDiff->{$trsX};
			my $nExonA = $infoExon->{$trs};
			my $nExonB = $infoExon->{$trsX};


			my $nnnA = $infoNNN->{$trs};
			my $nnnB = $infoNNN->{$trsX};
			if(($nExonA < 5 && $nExonB < 5) && !$nnnA && !$nnnB ){
			    print $ofh  "$trs and $trsX map to $scaff in the samearea ($stra, $diffA  and $strb, $diffB)\n";
			}
			last ;
		}
	}
}

sub ProcessOneGroup{
	my ($infoScaff,$infoRegion,@trslist) = @_ ;
	my $firsttrs = $trslist[0];
	my $scaffs = {};
	foreach my $trs (@trslist){
		my $scaff = $infoScaff->{$trs} ;
		if(! exists $scaffs->{$scaff}){
			$scaffs->{$scaff} = {};
		}
		$scaffs->{$scaff}->{$trs} =1 ;
	}
	my @uniqscaffs = (sort keys %{$scaffs});
	my $Nscaff = @uniqscaffs ;

	## process TRS within each scaffold
	my $diffloci = 0 ;
	my $mappedignore = {};
	foreach my $scaff (@uniqscaffs){
		my @trsL = (sort keys %{$scaffs->{$scaff}});
		## pairwise
		while(@trsL){
			my $trsa = shift @trsL ;
			next if(exists $mappedignore->{$trsa});
		    $diffloci++;
			foreach my $trsb (@trsL){
				my $val = IsTRSMappedtoSameScaffold($trsa,$trsb,$infoScaff,$infoRegion);
				die "$trsa,$trsb" if($val eq 0); ## just check

				if($val eq "samearea"){
					$mappedignore->{$trsb} = 1;
				}
			}
		}
	}
	print $ofh "$firsttrs $Nscaff unique scaffolds, diffloci $diffloci \n";

}

sub IsTRSMappedtoSameScaffold{
	my ($trs1,$trs2,$infoScaff,$infoRegion) = @_ ;
    my $scaff1 = $infoScaff->{$trs1} ;
    my $scaff2 = $infoScaff->{$trs2} ;
	if($scaff1 ne $scaff2){
		return 0  ;
	}

	my $startend1 = $infoRegion->{$trs1};
	my $startend2 = $infoRegion->{$trs2};

	my ($s1,$e1) = split " " , $startend1 ;
	my ($s2,$e2) = split " " , $startend2 ;
	if(($s1 >= $s2 && $s1 <= $e2) || ($s2 >= $s1 && $s2 <= $e1)){
		return "samearea" ;
	}
	else{
		return "diffarea" ;
	}
}

sub FixName{
	my ($i) = @_ ;
	$i =~ s/.ORF.*//;
	$i =~ s/.MER.*//;
	$i =~ s/_(A|B|C)//;
	return $i ;
}

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
