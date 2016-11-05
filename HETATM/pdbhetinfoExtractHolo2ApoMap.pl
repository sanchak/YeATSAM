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
usage( "Need to give a mapping -option -mapping  ") if(!defined $mapping);
usage( "Need to give a infile pdb id -option -infile  ") if(!defined $infile);
#usage( "Need to give a ispolar pdb id -option -ispolar  ") if(!defined $ispolar);
my ($RESULTDIR,$PDBDIR) = util_SetEnvVars();

my $ignoreHET = Config_HetIgnore();

my $INFOMAPPING = {};
my $GRPMAPPING = {};
my $grpnumber = 0 ;
if(defined $mapping){
  my $ifhmapping = util_read($mapping);
  while(<$ifhmapping>){
	  my @l = split ;
	  $grpnumber++;
	  foreach my $x (@l){
		  $INFOMAPPING->{$x} = $grpnumber ;
	  }
	  $GRPMAPPING->{$grpnumber} = \@l ;
  }
  close($ifhmapping);
}

my ($pdbNOHET) = util_maketablefromfile("$PDBDIR/NOHET");

my $IFH = util_read($infile);
while(<$IFH>){
   my ($hetatm,$N,@pdbs) = split ;
   ## do mapping analysis only for imp hetatms
   next if(exists $ignoreHET->{$hetatm});

   foreach my $x (@pdbs){
   	   next if(! exists $INFOMAPPING->{$x});
   	   my $grpnumber = $INFOMAPPING->{$x};
	   my @grp = @{$GRPMAPPING->{$grpnumber}};
       foreach my $y (@grp){
			if(exists  $pdbNOHET->{$y}){
			    print " YAHOO $hetatm $x $y \n";
				print "hetatmProcNew.pl -het $hetatm -pr $x -ispolar 1\n";
				print "hetatmHolo2ApoMatchSeqNum.pl -inf HETINFO/$x/$hetatm/$x.4.1.table.in -p1 $x -p2 $y\n";
			}
		}
   }
}
    

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}

