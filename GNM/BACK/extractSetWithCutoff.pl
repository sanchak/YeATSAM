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
my ($infile,$p1,$p2,$outfile,$cutoff,$blastdir,$listfile,$protein);
my (@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
GetOptions(
            "blastdir=s"=>\$blastdir ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "verbose=i"=>\$verbose ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);

usage( "Need to give a input file name => option -infile ") if(!defined $infile);
usage( "Need to give a input file name => option -cutoff ") if(!defined $cutoff);
usage( "Need to give a input file name => option -blastdir ") if(!defined $blastdir);
my $ifh = util_read($infile);
my $CNT = 0 ; 
#my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
 my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP,$BLASTOUT,$BLASTDB) = util_SetEnvVars();

my $info = {};
my $ofhstrippedlist = util_write("list.withcutoff$cutoff");
while(<$ifh>){
	 my ($trs,@l) = split ; 
	 print $ofhstrippedlist "$trs ";
	 while(@l){
	 	my $nm = shift @l ;
	 	my $v = shift @l ;
		if($v > $cutoff){
			$info->{$nm} = [] if(!defined  $info->{$nm});
			push @{$info->{$nm}} ,$trs ;
	                print $ofhstrippedlist "$nm ";
		}
	 }
	 print $ofhstrippedlist "\n";
}


my $ofhexec = util_write("execPW.$cutoff.csh");
my $ofhAllScaff = util_write("allScaff.$cutoff.csh");
my $ignoredblast = 0 ;

foreach my $nm (keys %{$info}){
	my @l = @{$info->{$nm}};
	print $ofh "$nm @l \n";


	my $allfound = 1 ;

	my $INFOFILE = "INFO/$nm.info";

	## Commands for generaring the BLASTS
	foreach my $xxx (@l){
		my $FILE = "$blastdir/$xxx.$nm.blast";
		if (! -e  "$FILE"){
		    print $ofhexec "parseBlastLatestpairwiseStep1.csh $xxx $nm\n";
		    $allfound = 0 ;
		}
		else{
			$ignoredblast++;
		}
	}

	if(!$allfound){
		if(-e "$INFOFILE"){
		   system("unlink $INFOFILE");
		}
	}

	## Commands for parsing the BLAST - and write out the info
	if($allfound && ! -e "INFO/$nm.info"){
	    print $ofhAllScaff "parseBlastLatestpairwiseOneScaffold.pl -scaffoldfasta SCAFFOLDDIR/$nm.ALL.1.fasta -blastdir $blastdir -justchecking 0 -sname $nm  ";
	    foreach my $xxx (@l){
	    	print $ofhAllScaff " -trs $xxx ";
	    }
	    print $ofhAllScaff " \n";
	}
}


system("wc -l execPW.$cutoff.csh");
print "ignoredblast $ignoredblast\n";
system("wc -l allScaff.$cutoff.csh");
print "Will write only those scaffolds which have all BLASTS done\n";



sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
