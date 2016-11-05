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


use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($ispolar,$infile,$outfile,$initial,$mapping,$hetatm,$listfile,$protein);
my (@expressions);
my $howmany = 4 ;
my $maxdist = 10 ;
my $verbose = 0 ;
my $INCR = 0.1 ;

my $IGNOREBACKBONE = 0 ;
my $doall = 1 ;
GetOptions(
            "mapping=s"=>\$mapping , "protein=s"=>\$protein , "infile=s"=>\$infile , "doall=i"=>\$doall , "listfile=s"=>\$listfile , "hetatm=s"=>\$hetatm , "outfile=s"=>\$outfile , "expr=s"=>\@expressions, "howmany=i"=>\$howmany , "maxdist=i"=>\$maxdist , "ispolar=i"=>\$ispolar , "initial=f"=>\$initial ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
## for debug 
usage( "Need to give a listfile -option -listfile  ") if(!defined $listfile);
#usage( "Need to give a protein pdb id -option -protein  ") if(!defined $protein);
#usage( "Need to give a ispolar pdb id -option -ispolar  ") if(!defined $ispolar);
my ($RESULTDIR,$PDBDIR) = util_SetEnvVars();
my @list= util_read_list_sentences($listfile);

my $ignoreHET = Config_HetIgnore();


my $pdbNOHET = {};
my $hetatmSIZES = {};
my $hetatm2PDBS = {};
my $PDB2hetatms = {};
my $pdbYESHET = {};
my @listREMOVEDNA ;
util_debug_printNList("Number of proteins being processed",@list);

my $ofhPDB2HET = util_write("$PDBDIR/PDB2HET");


foreach my $protein (@list){
    my $pdbinfo = "PDBINFO/$protein.hetinfo";
	if(! -e $pdbinfo){
		next;
	}
	my @list = util_read_list_words($pdbinfo);
	print $ofhPDB2HET "@list\n";
	shift @list ;
	foreach my $i (@list){
		my ($hetatm,$ninstances,$natoms) = split ":", $i ;

		if(!exists $ignoreHET->{$hetatm}){
		     $pdbYESHET->{$protein} = 1 ;
		}

		$hetatmSIZES->{$hetatm} = $natoms if(! defined $hetatmSIZES->{$hetatm} || $hetatmSIZES->{$hetatm} < $natoms);
		$hetatm2PDBS->{$hetatm} = {} if(!defined $hetatm2PDBS->{$hetatm});
		$hetatm2PDBS->{$hetatm}->{$protein} = 1 ;
	}
}

my $ofhHET2PDB = util_write("$PDBDIR/HET2PDB");

foreach my $protein (@list){
	$pdbNOHET->{$protein} = 1 if(! exists $pdbYESHET->{$protein});
}
my $ofhYESHET = util_write("$PDBDIR/YESHET");
my $ofhNOHET = util_write("$PDBDIR/NOHET");
foreach my $k (keys %{$pdbYESHET}){
	print $ofhYESHET "$k\n";
}
foreach my $k (keys %{$pdbNOHET}){
	print $ofhNOHET "$k\n";
}

print "Writing HET2PDB\n";
my @sortedSizehetatms = (sort {$hetatmSIZES->{$b} <=> $hetatmSIZES->{$a}} keys %{$hetatmSIZES});
foreach my $hetatm (@sortedSizehetatms){
   my $N = $hetatmSIZES->{$hetatm};
   my $tablePDBS = $hetatm2PDBS->{$hetatm} ;
   my @pdbs = (sort keys %{$tablePDBS});
   print $ofhHET2PDB "$hetatm $N @pdbs\n";

}
system("wc -l $PDBDIR/*HET* ");
    

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}

