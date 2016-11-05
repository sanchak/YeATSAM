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
my ($fastafile,$idx,$infile,$postfix,$p1,$p2,$outfile,$cutoff,$which_tech,$mapfile,$protein);
my ($ignorefile,@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
my $INDEXforNonChaperones = 1 ;
GetOptions( "which_tech=s"=>\$which_tech , "protein=s"=>\$protein , "fastafile=s"=>\$fastafile , "postfix=s"=>\$postfix , "infile=s"=>\$infile , "p1=s"=>\$p1 , "p2=s"=>\$p2 , "mapfile=s"=>\$mapfile , "ignorefile=s"=>\$ignorefile , "outfile=s"=>\$outfile , "expr=s"=>\@expressions, "howmany=i"=>\$howmany , "idx=i"=>\$idx , "verbose=i"=>\$verbose , "cutoff=f"=>\$cutoff ,);
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
my $ifh = util_read($infile);

usage( "Need to give a mapfile -option -mapfile  ") if(!defined $mapfile);
my $CNT = 0 ; 
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP,$BLASTOUT,$BLASTDB) = util_SetEnvVars();
my $aaconfig = new AAConfig("$SRC/aa.config");

my $PWD = cwd;

## this is to identify the chaperones
#my @expected = qw (H3 H4 H2A H2B);
my ($expected) = util_maketablefromfile("expected");


my $originalmapping = {};
my $ifhmap = util_read($mapfile);
## this is the name map file - first time you have create with hand,
## next time you provide the one that has been created... (which might have an _ )
while(<$ifhmap>){
     next if(/^\s*$/);
	 next if(/^\s*#/);

	 my ($nm,$namemap) = split ; 
	 $namemap =~ s/_.*//;
	 $originalmapping->{$nm} = $namemap ;
}

my $tableName2Mappedname = {};

my $counter = {};

## This is the grouped file - 
while(<$ifh>){
     next if(/^\s*$/);
	 next if(/^\s*#/);
	 s/.*=N//;
	 s/#.*//;

	 my ($nm,@l) = split ; 

	 ## the first one needs to be mapped
	 die "Is not mapped $nm" if(! exists $originalmapping->{$nm});
	 my $mappedname = $originalmapping->{$nm};
	 print "Info: first mapped $nm to $mappedname\n";

	 $tableName2Mappedname->{$nm} = $mappedname ;
	 foreach my $x (@l){
	 	if(exists $originalmapping->{$x}){
			die "if($originalmapping->{$x} ne $mappedname) for $x"  if($originalmapping->{$x} ne $mappedname);
	 	    $tableName2Mappedname->{$x} = $mappedname ;
		}
		else{
			print "\tInfo: Adding new mapping for $x $mappedname\n";
	 	    $tableName2Mappedname->{$x} = $mappedname ;
		}
	 }

	 ## this is for printing later - multiple chains of the same type need a index
	 $counter->{$mappedname} = 1 ;
}


my $last ;
my $done = {};

my $CHAPERONE = {};
my $CHAPERONEPDBs = {};
foreach my $fiveLetterPDB (sort keys %{$tableName2Mappedname}){
	my $fourletterPDB = $fiveLetterPDB ;
	$fourletterPDB =~ s/.$//;

	## this is to blockify each PDB
	if(!defined $last){
		$last = $fourletterPDB ;
	}
	else{
		if($last ne $fourletterPDB){
			print $ofh "\n\n";
			$last = $fourletterPDB ;
			foreach my $i (keys %{$counter}){
				$counter->{$i} = 1;
			}
		}
	}

	## find the chaperones - but print later
	my $indexlessName = $tableName2Mappedname->{$fiveLetterPDB}; 
	my $globalindex ;

	## Chaperones
	if(! exists $expected->{$indexlessName}){
		## Each pdb can have multiple chaperones - so need to add both to make a key
		my $KET = "$fourletterPDB $indexlessName";

		if(!exists $done->{$KET}){
			$CHAPERONE->{$indexlessName} = [] if(!defined $CHAPERONE->{$indexlessName} );
			push @{$CHAPERONE->{$indexlessName}}, $fourletterPDB ;
			$CHAPERONEPDBs->{$fourletterPDB} = 1 ;
	        my @l = @{$CHAPERONE->{$indexlessName}};
	        $globalindex = @l ;
			$done->{$KET} = $globalindex ;
		}
		else{
	        $globalindex =  $done->{$KET};
		}
	}
	else{
		### these are the expected ones - like H3 H4...
		$globalindex = $INDEXforNonChaperones++;
	}

	## final printing of the mapping
	my $idx  = $counter->{$indexlessName} ;
	$counter->{$indexlessName} = $counter->{$indexlessName} + 1 ;

	#my $finalstr = "$indexlessName" . "_". $globalindex . "_C" . $idx ;
	my $finalstr = "$indexlessName" . "_C" . $idx ;
	print $ofh "$fiveLetterPDB $finalstr\n";

}

my $OFHNoChaperones = util_write("WORKDIR/list.nochaperones");
my $tableNoChaperones = {};
foreach my $fiveLetterPDB (sort keys %{$tableName2Mappedname}){
	my $fourletterPDB = $fiveLetterPDB ;
	$fourletterPDB =~ s/.$//;
	$tableNoChaperones->{$fourletterPDB} = 1  if(! exists $CHAPERONEPDBs->{$fourletterPDB});
}

my $fhcommandsnochaperones = util_write("commands.nochaperones.csh");
#print $fhcommandsnochaperones "cd WORKDIR\n";
my $idx_nochaperone = 0 ;
foreach my $fourletterPDB (sort keys %{$tableNoChaperones}){
    #print $OFHNoChaperones "$fourletterPDB\n";
	my $str = "NOCHAPERONE_" . $idx_nochaperone++;
    #print $fhcommandsnochaperones "echo $fourletterPDB > !  list.$str \n";
    #print $fhcommandsnochaperones "pdbgetlist.csh list.$str $str\n";

		print $fhcommandsnochaperones "pdbgetlist.csh $fourletterPDB\n";
		print $fhcommandsnochaperones "\\cp -f PDBINFO/$fourletterPDB.contactchains WORKDIR/$str.tmp\n";
		print $fhcommandsnochaperones "fixNameMap.pl -outf WORKDIR/$str.contactchains -inf WORKDIR/$str.tmp -map WORKDIR/namemap\n";
}
system("wc -l WORKDIR/list.nochaperones");


## Note that we are printing only those which have at least one chaperone
my $fhcommands = util_write("commands.chaperones.csh");
my $ofhextractinfo = util_write("commands.extractinfo.csh");
#print $fhcommands "cd WORKDIR\n";
foreach my $k (sort keys %{$CHAPERONE}){
	my @l = @{$CHAPERONE->{$k}};
	my $N = @l -1 ;
	foreach my $i (0..$N){
	    my $idx = $i + 1;
		my $str = $k.  "_". $idx ;

		my $PDBNM = $l[$i];
		#print $fhcommands "echo $l[$i] > !  list.$str \n";
		#print $fhcommands "pdbgetlist.csh list.$str $str\n";
		print $fhcommands "pdbgetlist.csh $PDBNM\n";
		print $fhcommands "\\cp -f PDBINFO/$PDBNM.contactchains WORKDIR/$str.tmp\n";
		print $fhcommands "fixNameMap.pl -outf WORKDIR/$str.contactchains -inf WORKDIR/$str.tmp -map WORKDIR/namemap\n";

		print $ofhextractinfo "pdbinteromeAI.pl -outf jjjjj -inf WORKDIR/$str.contactchains -list ll -pr HSP90\n";
	}
}
system(" wc -l commands*.csh ");

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
