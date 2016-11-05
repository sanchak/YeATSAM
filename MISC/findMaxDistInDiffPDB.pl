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

use AAConfig;

use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($infile,$radii,$p1,$p2,$outfile,$which_tech,$hetatm,$ahbsfile,$protein);
my (@expressions);
my $cutoff = 20 ;
my $cutoffCAAtoms = 10 ;
my $verbose = 1 ;
my $fordna = 0 ;
my $autodetect = 1 ;
GetOptions(
            "fordna=i"=>\$fordna ,
            "autodetect=i"=>\$autodetect ,
            "protein=s"=>\$protein ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "radii=f"=>\$radii ,
            "infile=s"=>\$infile ,
            "ahbsfile=s"=>\$ahbsfile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
#my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT) = util_SetEnvVars();
my $PWD = cwd;
my $aaconfig = new AAConfig("$SRC/aa.config");

usage( "Need to give a protein 1 id -option -p1  ") if(!defined $p1);
usage( "Need to give a protein 2 id -option -p2  ") if(!defined $p2);
my $CNT = 0 ; 
system ("mkdir -p MAXDIST");
$outfile = "MAXDIST/$p1.$p2.maxdist.out" if(!defined $outfile);
print STDERR "Writing to $outfile\n";
my $ofh = util_write($outfile);

die "Same PDB" if ($p1 eq $p2);

my $file1 = "$PDBDIR/$p1.pdb";
my $file2 = "$PDBDIR/$p2.pdb";
my $pdb1 = new PDB();
$pdb1->ReadPDB($file1);
my $pdb2 = new PDB();
$pdb2->ReadPDB($file2);

my @CA1 = GetCA($pdb1);
my @CA2 = GetCA($pdb2);



### autodetect DNA PDBs

my $bothdna = 0 ;
if($autodetect){
   my $NCA1 = @CA1 ;
   my $NCA2 = @CA2 ;
   if(!($NCA1 && $NCA2)){
	   $fordna = 1 ;
   }
   if($NCA1 eq 0 && $NCA2 eq 0){
   	 $bothdna = 1;
   }
}

if($bothdna){
    die "Info: Both are DNA - do nothing";	
}

my ($infoAHBS);
if(defined $ahbsfile){
   ($infoAHBS) = util_HELIXParseAHBS($ahbsfile);
}


### this just selects residues that have CA between 20 A first
### no point in doing all atoms
my @allatoms1 ;
my @allatoms2 ;
my $done1 = {};
my $done2 = {};

if(!$fordna){
    while(@CA1){
		    my $a1 = shift @CA1 ;
		    foreach my $a2 (@CA2){
		           my $d = $pdb1->DistanceAtoms($a1,$a2); 
				   if($d < $cutoffCAAtoms){
				   #if(1){
				   	   my $N1 = $a1->GetResNum();
				   	   my $N2 = $a2->GetResNum();
					   my $r1 = $pdb1->GetResidueIdx($N1);
					   my $r2 = $pdb2->GetResidueIdx($N2);

					   if(! exists $done1->{$N1}){
		                   my @atoms1 = $r1->GetAtoms();
		                   push @allatoms1, @atoms1;
						   $done1->{$N1} = 1 ;
					   }
					   if(! exists $done2->{$N2}){
		                   my @atoms2 = $r2->GetAtoms();
		                   push @allatoms2, @atoms2;
						   $done2->{$N2} = 1 ;
					  }

				   }
	        }

    }
}
else{
    @allatoms1 = $pdb1->GetAtoms();
    @allatoms2 = $pdb2->GetAtoms();
}

my $N1 = @allatoms1 ;
my $N2 = @allatoms2 ;
print "There are $N1 and $N2 atoms ..... fordna=$fordna cutoffCAAtoms=$cutoffCAAtoms A\n";
my ($sorted,$minvalue,$A,$B,$AB) = util_GetDistancesBetween2SetsOfAtoms($p1,$p2,$pdb1,$pdb2,\@allatoms1,\@allatoms2,$cutoff);


my $ofhclasp1 = util_write("clasp.in.1");
my $ofhclasp2 = util_write("clasp.in.2");
my $DONE1 = {};
my $DONE2 = {};

my $ELECTABLE={};
my $ELECBasiccnt = 0 ;
my $ELECAcidiccnt = 0 ;
foreach my $k (sort {$sorted->{$a} <=> $sorted->{$b}} keys %{$sorted}){
	my @l = split " ",$k ;
	my $pdb1 = shift @l ;
	my $pdb2 = shift @l ;
	my $atom1 = shift @l ;
	my $atom2 = shift @l ;
	my $dist = shift @l ;
	my ($res1,$num1,$type1) = split "/", $atom1 ;
	my ($res2,$num2,$type2) = split "/", $atom2 ;
	my $exists1 = IsAHBS($pdb1,$infoAHBS,$num1);
	my $exists2 = IsAHBS($pdb2,$infoAHBS,$num2);
	next if($res1 eq "HOH" || $res2 eq "HOH");

	my $IsElectrostatic = 0 ;
	if(!$fordna){
	$IsElectrostatic = $aaconfig->IsElectrostatic($res1,$res2);
	if($IsElectrostatic ne 0){
		my $str = $res1.$num1 . "<->" . $res2. $num2 ;
		if(! exists $ELECTABLE->{$str} && $dist < 4){
			$ELECTABLE->{$str} = 1;
			if($IsElectrostatic>0){
			   $ELECBasiccnt++;
			}
			else{
			   $ELECAcidiccnt++;
			}
		}
	}
	}

	my @lll = split " ", $k ;
	printf $ofh ("%6s%6s%18s%18s%10s%8s%8s%5s\n",@lll,$exists1,$exists2,$IsElectrostatic);
	#print $ofh "$k $exists1 $exists2 $IsElectrostatic\n"; 

}
my $Pstr = "# $ELECBasiccnt = ELECBasiccnt , $ELECAcidiccnt = ELECAcidiccnt  ";
foreach my $k (keys %{$ELECTABLE}){
	print "$k \n";
	$Pstr = $Pstr . " $k " ;
}
print "$Pstr \n";


sub IsAHBS{
	my ($pdb,$info,$resnum) = @_ ;
	if(exists $info->{$pdb}){
		my $tab = $info->{$pdb} ;
		if(exists $info->{$pdb}->{$resnum}){
			return $info->{$pdb}->{$resnum};
		}
	}
	return "X";
}


## turned off for now - cant see the use 
if(0){

print "For first: ";
foreach my $k (sort {$a <=> $b} keys %{$A}){
	my $v = $A->{$k};
	print "$v$k ";
}
print "\n";

print "For second: ";
foreach my $k (sort {$a <=> $b} keys %{$B}){
	my $v = $B->{$k};
	print "$v$k ";
}

print "\n";
print "For BOTH: ";
foreach my $k (keys %{$AB}){
	my $v = $AB->{$k};
	print "$k ";
}
print "\n";
}

close($ofh);


sub GetAllAtoms{
	my ($P) = @_ ;
    my @res = $P->GetResidues();
    my @allatoms ;
   foreach my $res (@res){
		   my @atoms1 = $res->GetAtoms();
		   push @allatoms, @atoms1;
   }
   return @allatoms ;
}
sub GetCA{
	my ($P) = @_ ;
    my @res = $P->GetResidues();
    my @allatoms ;
   foreach my $res (@res){
		   #my @atoms1 = $res->GetAtoms();
		   my $a = $res->GetAtomType("CA");
		   push @allatoms, $a if(defined $a);
   }
   return @allatoms ;
}

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
