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
            "mapping=s"=>\$mapping , "protein=s"=>\$protein , "infile=s"=>\$infile , "doall=i"=>\$doall , "listfile=s"=>\$listfile , "hetatm=s"=>\$hetatm , "outfile=s"=>\$outfile , "expr=s"=>\@expressions, "howmany=i"=>\$howmany , "maxdist=i"=>\$maxdist , "ispolar=i"=>\$ispolar , "initial=f"=>\$initial , "verbose=i"=>\$verbose,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
## for debug 
usage( "Need to give a protein -option -protein  ") if(!defined $protein);
#usage( "Need to give a protein pdb id -option -protein  ") if(!defined $protein);
#usage( "Need to give a ispolar pdb id -option -ispolar  ") if(!defined $ispolar);
my ($RESULTDIR,$PDBDIR) = util_SetEnvVars();

my $ignoreHET = Config_HetIgnore();

my $pdbNOHET = {};
my $hetatmSIZES = {};
my $hetatm2PDBS = {};
my $PDB2hetatms = {};
my $pdbYESHET = {};
my @listREMOVEDNA ;
my $ofhPDB2HET = util_write("PDBINFO/$protein.hetinfo");
print "Writing PDBINFO/$protein.hetinfo\n" if($verbose);

{
    my $pdbfile = "$PDBDIR/$protein.pdb";
    my $pdb1 = new PDB($pdbfile);
	if($pdb1->IsDNA()){
		print "Info: Nothing to be done for DNA $protein\n";
		exit ;
	}
	push @listREMOVEDNA,$protein ;

    my $infohetsTOP = $pdb1->FindHetAtms();
	foreach my $hetatm (keys %{$infohetsTOP}){

		 ## do not ignore anything except water 
		 #next if(exists $ignoreHET->{$hetatm});
		 next if($hetatm eq "HOH");

         my $infoHET = $infohetsTOP->{$hetatm};
		 my $N = keys %{$infoHET};
		 my $str  = "";
		 foreach my $hetNUM (keys %{$infoHET}){
			## number of atoms in this particular HETATM
			my $NN = 0 ;
			foreach my $i (@{$infoHET->{$hetNUM}}){
				my ($j1,$j2,$XX) = split " ", $i ;
				$NN++ if($XX !~ /^H/);
			}
		 	#my $NN = @{$infoHET->{$hetNUM}};
			$hetatmSIZES->{$hetatm} = $NN if(!defined $hetatmSIZES->{$hetatm} || $hetatmSIZES->{$hetatm} < $NN );
			$str = $str . " $NN " ;
		 }
		 
		 #print "$protein $hetatm $N ( $str )\n";
		 $pdbYESHET->{$protein} =1 ;
		 $PDB2hetatms->{$protein} = {}  if(! defined $PDB2hetatms->{$protein});
		 $PDB2hetatms->{$protein}->{$hetatm} = $N ;

		 $hetatm2PDBS->{$hetatm} = {} if(!defined $hetatm2PDBS->{$hetatm});
		 $hetatm2PDBS->{$hetatm}->{$protein} = 1;
	}
}


	   my $str2print = "";
   	   my $hetatmtable =  $PDB2hetatms->{$protein};
	   foreach my $hetatm (sort keys %{$hetatmtable}){
	   		 my $size = $hetatmSIZES->{$hetatm};
	   	     my $ninstances = $hetatmtable->{$hetatm} ;
			 $str2print = $str2print . " $hetatm:$ninstances:$size ";
	   }
       print $ofhPDB2HET "$protein $str2print \n";
       print  "$protein $str2print \n" if($verbose);

    

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}

