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
use AHB ;


use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($infile,$radii,$outfile,$values,$which_tech,$hetatm,$listfile,$protein);
my (@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
my $isCAOnly = 0 ;
GetOptions(
            "isCAOnly=i"=>\$isCAOnly ,
            "protein=s"=>\$protein ,
            "radii=f"=>\$radii ,
            "infile=s"=>\$infile ,
            "values=s"=>\$values ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a protein pdb id -option -values  ") if(!defined $values);
usage( "Need to give a protein pdb id -option -protein  ") if(!defined $protein);
usage( "Need to give a outfile pdb id -option -outfile  ") if(!defined $outfile);
my $ofh = util_write($outfile);
my $CNT = 0 ; 
#my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP,$BLASTOUT,$BLASTDB) = util_SetEnvVars();
my $PWD = cwd;

my $pdb = "$PDBDIR/$protein.pdb";
my $pdb1 = new PDB();
$pdb1->ReadPDB($pdb);



my ($RESTABLE,$AHinfo,$AHList,$Binfo,$BList,$DISUPLH,$HM) = ahb_ParseValues($protein,$pdb1,$values,$isCAOnly);

my @AHList = @{$AHList};
my @BList = @{$BList};

my $NAH = @AHList;
my $NB = @BList;



#die "Error: Did not find protein $protein in $values file" if(!@AHList || !@BList);


my $infoAHlist = ahb_printPairWiseDistance($protein,$pdb1,$RESTABLE,@AHList);
my $infoBlist = ahb_printPairWiseDistance($protein,$pdb1,$RESTABLE,@BList);
my $infoAHBLIST  = ahb_printCrossWiseDistance($protein,$pdb1,$RESTABLE,\@BList,\@AHList);

print $ofh "NAH $NAH \n";
print $ofh "NB $NB \n";
ahb_PrintData($infoAHlist,$ofh);
ahb_PrintData($infoBlist,$ofh);
ahb_PrintData($infoAHBLIST,$ofh);

## Note disulphide bonds are not exactly 2 
util_GetDisulphide ($protein,$pdb1,0);
my @disulphide = $pdb1->GetDisulphide();
my $DONEDISULPHIDE = {};
while(@disulphide){
	my $a = shift @disulphide ;
	my $b = shift @disulphide ;
	my @sort ;
	if(exists $DISUPLH->{$a} && exists $DISUPLH->{$b}){
		push @sort, $DISUPLH->{$a} ;
		push @sort, $DISUPLH->{$b} ;
		my @sorted = sort @sort ;
		my $str = $sorted[0] . " " .  $sorted[1];
		next if(exists $DONEDISULPHIDE->{$str});
		$DONEDISULPHIDE->{$str} = 1;
		print $ofh "DISUPLH $str \n";
	}
}








sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
