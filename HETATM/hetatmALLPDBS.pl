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

usage( "Need to give a listfile -option -listfile  ") if(!defined $listfile);


my ($RESULTDIR,$PDBDIR) = util_SetEnvVars();
my @list= util_read_list_sentences($listfile);

my $ignoreHET = Config_HetIgnore();

my $ifhmapping = util_read($mapping);
my $INFOMAPPING = {};
my $GRPMAPPING = {};
my $grpnumber = 0 ;
while(<$ifhmapping>){
	my @l = split ;
	$grpnumber++;
	foreach my $x (@l){
		$INFOMAPPING->{$x} = $grpnumber ;
	}
	$GRPMAPPING->{$grpnumber} = \@l ;
}

my $pdbNOHET = {};
my $hetatmSIZES = {};
my $hetatm2PDBS = {};
my $PDB2hetatms = {};
my $pdbYESHET = {};
my @listREMOVEDNA ;
util_printNList("Number of proteins being processed",@list);
foreach my $protein (@list){
    my $pdbfile = "$PDBDIR/$protein.pdb";
	if(! -e $pdbfile){
		my $ofhmissing = util_open_or_append("missingfiles");
		print $ofhmissing "$protein\n";
		next ;
	}
    my $pdb1 = new PDB($pdbfile);
	next if($pdb1->IsDNA());
	push @listREMOVEDNA,$protein ;

    my $infohetsTOP = $pdb1->FindHetAtms();
	foreach my $hetatm (keys %{$infohetsTOP}){
		 next if(exists $ignoreHET->{$hetatm});

         my $infoHET = $infohetsTOP->{$hetatm};
		 my $N = keys %{$infoHET};
		 my $str  = "";
		 foreach my $hetNUM (keys %{$infoHET}){
			## number of atoms in this particular HETATM
		 	my $NN = @{$infoHET->{$hetNUM}};
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

my $ofhPDB2HET = util_write("$PDBDIR/PDB2HET");
my $ofhHET2PDB = util_write("$PDBDIR/HET2PDB");

print "Writing PDB2HET\n";
foreach my $protein (@listREMOVEDNA){
   if(! exists $pdbYESHET->{$protein}){
       $pdbNOHET->{$protein} = 1 ;
   }
   else{
	   my $str2print = "";
   	   my $hetatmtable =  $PDB2hetatms->{$protein};
	   foreach my $hetatm (sort keys %{$hetatmtable}){
	   		 my $size = $hetatmSIZES->{$hetatm};
	   	     my $ninstances = $hetatmtable->{$hetatm} ;
			 $str2print = $str2print . " $hetatm:$ninstances:$size ";
	   }
       print $ofhPDB2HET "$protein $str2print \n";
   }
}

print "Writing HET2PDB\n";
my @sortedSizehetatms = (sort {$hetatmSIZES->{$b} <=> $hetatmSIZES->{$a}} keys %{$hetatmSIZES});
foreach my $hetatm (@sortedSizehetatms){
   my $N = $hetatmSIZES->{$hetatm};
   my $tablePDBS = $hetatm2PDBS->{$hetatm} ;
   my @pdbs = (sort keys %{$tablePDBS});
   print $ofhHET2PDB "$hetatm $N @pdbs\n";

   foreach my $x (@pdbs){
   	   next if(! exists $INFOMAPPING->{$x});
   	   my $grpnumber = $INFOMAPPING->{$x};
	   my @grp = @{$GRPMAPPING->{$grpnumber}};
       foreach my $y (@grp){
			if(exists  $pdbNOHET->{$y}){
			    print " YAHOO $hetatm $x $y \n";
			}
		}
   }
}
system("wc -l $PDBDIR/PDB2HET $PDBDIR/HET2PDB");
    

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}

