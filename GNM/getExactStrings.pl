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
my $verbose = 0 ;
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
my $ofh = util_write("$outfile.listnew");
my $ofhsubset = util_write("$outfile.subset");
usage( "Need to give a input file name => option -fastafile ") if(!defined $fastafile);
my ($tablefasta,$infoSeq2PDB,$mapChainedName2Name) = util_parsePDBSEQRESNEW($fastafile,0,0);
 my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP,$BLASTOUT,$BLASTDB) = util_SetEnvVars();
my $aaconfig = new AAConfig("$SRC/aa.config");
my  ($tableThree2One,$tableOne2Three,@sortedSingle) = Config_AACodes();


my $strall = "";
my @listofinit ;
foreach my $k (keys %{$tablefasta}){
	push @listofinit, $k ;
	 my $v = $tablefasta->{$k} ; 
	 $strall = $strall . "ZZZ" . $v ;
}


my $doneremoving = {};

    my $remove = {};
    foreach my $k (keys %{$tablefasta}){
		 next if(exists $doneremoving->{$k});
	     my $v = $tablefasta->{$k} ; 
	     my @l = ($strall =~ /$v/g);
	     my $N = @l ;
	     next if($N eq 1);
		 
	     print "$v is in many places $N\n" if($verbose);
	     $remove->{$k} = $v ;
	     $doneremoving->{$k} = $v ;
    }


my @k = (keys %{$remove});
my $NN = @k ;
print STDERR "Info:There were $NN exact matches\n";

my $ofhDB = util_write("$fastafile.exact");
foreach my $i (@listofinit){
	if(exists $remove->{$i}){
		my $v = $remove->{$i} ;
		print $ofhsubset "$i $v\n";
	}
	else{
	     print $ofh "$i\n";
	     my $v = $tablefasta->{$i} ; 

		 print $ofhDB ">$i\n";
		 print $ofhDB "$v\n";
	}
}
system("wc -l $outfile.*");


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
