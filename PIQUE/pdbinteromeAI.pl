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
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
if(-z $infile){
	print "Info: input file $infile is empty - possibly no chains. Exiting\n";
	exit(1);
}
my $ifh = util_read($infile);

usage( "Need to give a listfile -option -listfile  ") if(!defined $listfile);
usage( "Need to give a protein pdb id -option -protein  ") if(!defined $protein);

my $CNT = 0 ; 
#my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
 my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP,$BLASTOUT,$BLASTDB) = util_SetEnvVars();
my $aaconfig = new AAConfig("$SRC/aa.config");

# for reading fasta
#my ($fe,$x) =  util_readfasta($infile);
#chomp $x ;
my $PWD = cwd;

system ("mkdir -p WORKDIR/PYMOLCOMMANDS");
system ("mkdir -p WORKDIR/CLASPIN");

my @expected = qw ( HSP90_C1<->RAR1_C1
HSP90_C1<->SGT1A_C2
HSP90_C2<->RAR1_C2
HSP90_C2<->SGT1A_C1
HSP90_C1<->SGT1A_C1
HSP90_C2<->SGT1A_C1
HSP90_C2<->SGT1A_C2
HSP90_C3<->SGT1A_C2
HSP90_C3<->SGT1A_C3
HSP90_C1<->CDC37_C1
);

my $expectedtable = util_make_table(\@expected);

my $ignoretable ;
if(defined $ignorefile){
    ($ignoretable) = util_maketablefromfile($ignorefile);
}


my @list= util_read_list_sentences($listfile);
my $list = {};
map { s/\s*//g ; $list->{$_} = 1 ; } @list ;

foreach my $i (@list){
}


my $FULLIDS = {};
my $PREFIX_2_FULLIDS = {};
my $ID_2_EDGES = {};
## Read the contact file...
my $PDB ;
while(<$ifh>){
     next if(/^\s*$/);
     next if(/^\s*#/);
     chomp ;
	 if(/EDGE/){
	 	s/=//g;
	    my ($x1,$jjj1,$jjj2,$A,$B,$H1,$H2,$res1,$res2) = util_ParseContactFile($_);
	    my ($x2,$PDB1,$PDB2,$copyA,$copyB) = util_ParseContactFile($_);
		$PDB = $PDB1 ;
		$copyA =~ s/_.*//;
		$copyB =~ s/_.*//;
		$FULLIDS->{$A} = $copyA ;
		$FULLIDS->{$B} = $copyB ;

		$PREFIX_2_FULLIDS = Add2PrefixMapping($PREFIX_2_FULLIDS,$copyA,$A);
		$PREFIX_2_FULLIDS = Add2PrefixMapping($PREFIX_2_FULLIDS,$copyB,$B);

		$ID_2_EDGES->{$A} = {} if(! exists $ID_2_EDGES->{$A});
		$ID_2_EDGES->{$A}->{$_} = 1 ;

		$ID_2_EDGES->{$B} = {} if(! exists $ID_2_EDGES->{$B});
		$ID_2_EDGES->{$B}->{$_} = 1 ;
	 }
}

$PDB =~ s/.$//;
print "THE PDB is $PDB\n";

foreach my $k (sort keys %{$PREFIX_2_FULLIDS}){
	if($k eq $protein){
	    my $t = $PREFIX_2_FULLIDS->{$k};
	    my $N = (keys %{$t});
	    my @l = (sort keys %{$t});
		$, = " and ";
	    print "$k has $N number of chains, @l. \n";


		## process each chain
		my $ALLDNA = {};
		foreach my $chain (@l){
			print "\n+++++ processing $chain +++++++++\n";
			my @edges = (sort keys %{$ID_2_EDGES->{$chain}});

			## find the chains with which this chain makes contacts
			my $contactchains = {};
			foreach my $edge (@edges){
	            my ($x1,$PDB1,$PDB2,$A,$B,$H1,$H2,$res1,$res2) = util_ParseContactFile($edge);
				#print "$edge llllll\n";
				my $x = $A =~ $chain? $A : $B ;
				my $y = $A =~ $chain? $B : $A ;
				my $PDBA = $A =~ $chain? $PDB1 : $PDB2 ;
				my $PDBB = $A =~ $chain? $PDB2 : $PDB1 ;
				my $concatPDB = "$PDBA.$PDBB";
				$contactchains->{$y} = $concatPDB;
			}

			my @contactchains = (sort keys %{$contactchains});
	        #print  "$chain makes contact to @contactchains . \n";

			foreach my $otherchain (sort keys %{$contactchains}){
				 my $concatPDB = $contactchains->{$otherchain};
			     ProcessOneChainMakingContact($chain,$otherchain,$concatPDB,@edges);
			}

		}
	}
}

sub ProcessOneChainMakingContact{
	        my ($chain,$otherchain,$concatPDB,@edges) = @_ ;
			my $Chainresidues = "";
			my $otherChainresidues = "";
			my $ChainHELICES = "";
			my $otherChainHELICES = "";
			foreach my $edge (@edges){
	            my ($x1,$PDB1,$PDB2,$A,$B,$H1,$H2,$res1,$res2) = util_ParseContactFile($edge);
				my $x = $A =~ $chain? $A : $B ;
				my $y = $A =~ $chain? $B : $A ;

				## Choose the correct Helices and residues
				my $HEL = $A =~ $chain? $H1 : $H2 ;
				my $otherHEL = $A =~ $chain? $H2 : $H1 ;
				my $res = $A =~ $chain? $res1 : $res2 ;
				my $otherres = $A =~ $chain? $res2 : $res1 ;

				## Choose the otherchain
				if($y =~ /$otherchain/){
					## only DNA will come twice, else each chain will have one
					die "$Chainresidues " if($Chainresidues ne "" && $y ne "DNA");
					$Chainresidues = $Chainresidues . $res ;
					$otherChainresidues = $otherChainresidues . $otherres ;

				    $ChainHELICES = $ChainHELICES . ",".  $HEL ;
				    $otherChainHELICES = $otherChainHELICES . ",".  $otherHEL ;
				}
			}

			my ($tableChainresidues) = util_CreateTableAndListFromSymbolSeperated($Chainresidues,",");
			my @tableChainresidues = SortThreeLetterAminoAcid(sort keys %{$tableChainresidues});

			my ($tableOtherChainresidues) = util_CreateTableAndListFromSymbolSeperated($otherChainresidues,",");
			my @tableOtherChainresidues = SortThreeLetterAminoAcid(sort keys %{$tableOtherChainresidues});

			my ($tableChainHelices) = util_CreateTableAndListFromSymbolSeperated($ChainHELICES,",");
			my @tableChainHelices = (sort keys %{$tableChainHelices});

			my ($tableOthereChainHelices) = util_CreateTableAndListFromSymbolSeperated($otherChainHELICES,",");
			my @tableOthereChainHelices = (sort keys %{$tableOthereChainHelices});

			my $str = "$chain<->$otherchain";
			if(! exists $expectedtable->{$str}){
				my $ofherr = util_open_or_append("errfile");
				print $ofherr "$str does not exist in expected for $PDB\n";
			}

			my $commandfile = -e "WORKDIR/PYMOLCOMMANDS/$concatPDB.csh"? "WORKDIR/PYMOLCOMMANDS/$concatPDB.csh" :  "See reverese";

			my $strA = join "," , @tableChainresidues;
			my $strB = join "," , @tableOtherChainresidues;
			my $strHB_A = join "," , @tableChainHelices;
			my $strHB_B = join "," , @tableOthereChainHelices;

	        #print  "\t$chain makes contact to $otherchain with residues (A=@tableChainresidues, B=@tableOtherChainresidues)  through helices/b-sheets @tableChainHelices ($commandfile) \n";
	        print  "\t$chain -> $otherchain (A=$strA  B=$strB)  (H/B: A=$strHB_A B=$strHB_B ) ($commandfile) \n";

			my $CLASP_1 = "WORKDIR/CLASPIN/${concatPDB}_1.csh";
			my $CLASP_2 = "WORKDIR/CLASPIN/${concatPDB}_2.csh";
			my $ofhCLASP1 = util_write($CLASP_1);
			my $ofhCLASP2 = util_write($CLASP_2);
			print $ofhCLASP1 "@tableChainresidues\n";
			print $ofhCLASP2 "@tableOtherChainresidues\n";

			## else, the other pair option has worked
			if(-e $commandfile){
			   my $ofhPymolAppend = util_append("WORKDIR/PYMOLCOMMANDS/$concatPDB.csh");
			   my ($tmptable) = util_maketablefromfile("WORKDIR/PYMOLCOMMANDS/$concatPDB.csh");
			   if(! exists $tmptable->{"createCLASPinput.csh"}){

			       my $NNN1 = @tableChainresidues;
			       my $NNN2 = @tableOtherChainresidues;
    
			       ## Get the two pdbs, hacky...
			       my $TMP = $concatPDB ;
			       $TMP =~ s/\./ /;
			       my @l = split " ", $TMP;
			       my $P1 = $l[0];
			       my $P2 = $l[1];
    
			       print $ofhPymolAppend "createCLASPinput.csh $P1 $CLASP_1 $NNN1 1\n";
			       print $ofhPymolAppend "createCLASPinput.csh $P2 $CLASP_2 $NNN2 1\n";
			   }
			}

			#print "RRRR $chain<->$otherchain\n";

			my $strtableChainresidues = join "," , @tableChainresidues;
			my $toprint = "$PDB $str $strtableChainresidues";
			$expectedtable->{$str} = $toprint ;

			my $strtableOtherChainresidues = join "," , @tableOtherChainresidues;
}

foreach my $s (@expected){
    my $fname = "WORKDIR/contacts.$s";
    my $OFHRES = util_open_or_append($fname);
    my ($tableoutput) = util_maketablefromfile_firstentry($fname);
    my $v = $expectedtable->{$s} ;
	if($v eq 1){
		$v = "$PDB $s X";
	}
	if(! exists $tableoutput->{$PDB}){
		print $OFHRES "$v\n";
	}
}

sub SortThreeLetterAminoAcid{
    my @l = @_ ;
	my $table = {};
	foreach my $l (@l){
		my $num = $l ;
		$num =~ s/...//;
		my ($three) = ($l =~ /(...)/);

		## hack for HISTONEs
		#if($num > 400){
			#$num = $num -400 ;
		#}
		my ($tableThree2One,$tableOne2Three) = Config_AACodes();
		die "$three does not exist @#@#@@#@ " if(! exists $tableThree2One->{$three});
		my $oneletter = $tableThree2One->{$three} ;
		$l = $oneletter . $num ;
		$table->{$l} = $num ;
	}
	my @ret = (sort {$table->{$a} <=> $table->{$b}} keys %{$table});
	return @ret ;
}


sub Add2PrefixMapping{
	my ($table,$prefix,$full) = @_;
	if(! exists $table->{$prefix}){
		$table->{$prefix} = {};
	}
	$table->{$prefix}->{$full} = 1 ;
	return $table ;

}


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
