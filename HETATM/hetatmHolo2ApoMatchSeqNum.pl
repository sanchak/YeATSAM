#!/usr/bin/perl -w 
use strict ;
use PDB;
use FileHandle ;
use Algorithm::Combinatorics qw(combinations) ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use ConfigPDB;
use MyGeom;
use MyPymol;
use Math::Geometry ;
use Math::VectorReal qw(:all);  # Include O X Y Z axis constant vectors




use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($p1,$p2,$infile,$outfile,$which_tech,$listfile,$protein);
my (@expressions);
my $maxdist = 3 ;
my $verbose = 1 ;
my ($moveZ,$verify,$decaaf,$before1,$before2);
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "verify"=>\$verify ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "infile=s"=>\$infile ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "maxdist=i"=>\$maxdist ,
            "moveZ=i"=>\$moveZ ,
            "decaaf=i"=>\$decaaf ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
#usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofhclose = util_write("log.close");
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
usage( "Need to give a protein 1 id -option -p1  ") if(!defined $p1);
usage( "Need to give a protein 2 id -option -p2  ") if(!defined $p2);
 my ($RESULTDIR,$PDBDIR) = util_SetEnvVars();

my $file1 = "$PDBDIR/$p1.pdb";
my $file2 = "$PDBDIR/$p2.pdb";
{
my $pdb1 = new PDB();
$pdb1->ReadPDB($file1);
my $pdb2 = new PDB();
$pdb2->ReadPDB($file2);


my $IFH = util_read($infile);
while(<$IFH>){
   my @list= split ;

   my $ret = MatchListResidueIn2PDBS_string($pdb1,$pdb2,@list);
   print "$ret oooooooo\n";
   
   my @residuesAtoms1 = GetResiduesAtomsList_string($pdb1,@list);
   my @residuesAtoms2 = GetResiduesAtomsList_string($pdb2,@list);
   my @ld1 = GetPairwiseDiffFromResidueSet($pdb1,@residuesAtoms1);
   print "kkkkkkk\n";
   my @ld2 = GetPairwiseDiffFromResidueSet($pdb2,@residuesAtoms2);
   while (@ld1){
	   my $d1 = shift @ld1 ;
	   my $d2 = shift @ld2 ;
	   my $diff = $d1 - $d2;
	   print "$diff \n";
   }
}
}

sub GetResiduesAtomsList_string{
	my ($pdb,@l) = @_ ;
	my @objs ;
	foreach my $r (@l){
	      my @LL = split "/", $r ;
		  my $A = shift @LL ;
		  my $type = shift @LL ;
	      my ($name,$idx) = ($r =~ /([a-zA-Z]+)([0-9]+)/);
	      my $i = $pdb->GetResidueIdx($idx);
		  my $atom = $i->GetAtomType($type);
		  push @objs, $atom ;
	}
	return @objs ;
}


sub MatchListResidueIn2PDBS_string{
	   my ($pdb1,$pdb2,@l) = @_ ;
	   my $ret = 1 ;
	   foreach my $r (@l){
	         $ret = $ret * MatchSingleResidueIn2PDBS_string($pdb1,$pdb2,$r);
	   }
	   return $ret ;
}

sub MatchSingleResidueIn2PDBS_string{
	   my ($pdb1,$pdb2,$resNameNumber) = @_ ;
	   my ($name,$idx) = ($resNameNumber =~ /([a-zA-Z]+)([0-9]+)/);
	   return MatchSingleResidueIn2PDBS_idx($pdb1,$pdb2,$idx);
}

sub MatchSingleResidueIn2PDBS_idx{
	   my ($pdb1,$pdb2,$idx) = @_ ;
	   my $i = $pdb1->GetResidueIdx($idx);
	   my $j = $pdb2->GetResidueIdx($idx);
	   if(!defined $i || !defined $j){
	   	  print "Did not find atom for $idx\n";
	   	  return 0 ;
	   }
	   else{
	   	  my $nmA = $i->GetNameFull();
	   	  my $nmB = $j->GetNameFull();
		  if($nmA eq $nmB){
		     return 1 ;
		  }
		  else{
		     print "names $nmB and $nmA are not equal\n";
		     return 0 ;
		  }
	   }
}

sub GetPairwiseDiffFromResidueSet{
    my ($pdb,@residues) = @_ ; 
    my $iter1 = combinations(\@residues, 2);

    my @l ; 
    while (my $c1 = $iter1->next) {
        my @combo1 = @{$c1} ; 
		my ($a1) = $combo1[0];
		my ($a2) = $combo1[1];
 
	    my $d1 = $pdb->DistanceAtoms($a1,$a2);
		push @l, $d1 ;
		print "D = $d1 \n";
	}
	return @l ;
}


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
