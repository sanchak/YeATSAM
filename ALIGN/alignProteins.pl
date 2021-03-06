#!/usr/bin/perl -w 
use strict ;
use PDB;
use FileHandle ;
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
my $CNT = 0 ; 
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
my $PWD = cwd;


my $file1 = "$PDBDIR/$p1.pdb";
my $file2 = "$PDBDIR/$p2.pdb";
my $pdb1 = new PDB();
$pdb1->ReadPDB($file1);
my $pdb2 = new PDB();
$pdb2->ReadPDB($file2);

if(defined $verify){
	$before1 = $pdb1->DistanceMatrix();
	$before2 = $pdb2->DistanceMatrix();
}

my $mutafh ; 
if(defined $decaaf){
    $mutafh = util_write("mutate.directions");
}

my $info = {};


#print "Parsing results file\n";
my ($atoms1,$atoms2) = pymolin_getResultsLine($infile,$pdb1,$pdb2);


#print "Aligning geom_Align3PointsToXYPlane \n";
my ($done1,$remainingAtoms1) = geom_Align3PointsToXYPlane($pdb1,$atoms1,$verbose);
my ($done2,$remainingAtoms2) = geom_Align3PointsToXYPlane($pdb2,$atoms2,$verbose);


my ($absscore,$normalizedSum) = $pdb1->ScoreGivenSetOfAtoms($pdb2,$done1,$done2);
print STDERR " SCORES : $absscore,$normalizedSum .........\n";


if(defined $verify){
	my $after1 = $pdb1->DistanceMatrix();
	my $after2 = $pdb2->DistanceMatrix();
	$pdb1->VerifyDistanceMatices($after1,$before1);
	$pdb2->VerifyDistanceMatices($after2,$before2);
}

my $done2Res = {};
foreach my $i (@{$done2}){
	$done2Res->{$i->GetResNum()} = 1 ;
}


#$pdb2->PrintSeq();

print STDERR "Doing each remaining atom\n\n\n";
foreach my $i (@{$remainingAtoms1}){

	my @tmp1 = (@{$done1}, $i);
	print STDERR "Remaining Atom\n";
	print $ofhclose "Remaining Atom\n";
	print STDERR "=====================\n";
	$i->Print();
	$i->Print("",$ofhclose);
	my @atomlist ;
	push @atomlist, $i ;
	my ($results,$combined) = $pdb2->GetNeighbourHoodAtom(\@atomlist,$maxdist);
	print $ofhclose "Atoms close to this one\n";
	print STDERR  "Atoms close to this one\n";
	my $sort ;
    foreach my $j (@{$combined}){
		my $resnum = $j->GetResNum(); 
	    my $d = $i->Distance($j) ;
		$sort->{$j} = $d ;
	    #print $ofhclose "XXXX = $d $resnum  \n";
	}
	my @sorted = sort { $sort->{$a} <=> $sort->{$b} } (@{$combined});


	my $locallydone = {};
	foreach my $j (@sorted){
	    next if(exists $done2Res->{$j->GetResNum()});
	    next if(exists $locallydone->{$j->GetResNum()} && $verbose == 0);

		my $resnum = $j->GetResNum(); 
		## just CA - for decaaf
		if(defined $decaaf){
		    my $dist = $sort->{$j};
		    my $type = $j->GetType(); 
		    next if($type ne "CA");
		    next if($dist > $decaaf);
			my $atomnmj = $j->GetName();
			my $atomnmi = $i->GetName();
			print $mutafh "$atomnmi $atomnmj $dist \n";
		}


	    #print $ofhclose "YYY $resnum  \n";

		$locallydone->{$j->GetResNum()} = 1  ; 
	    my @tmp2 = (@{$done2}, $j);
		print STDERR " dist = $sort->{$j} \n";
		print $ofhclose " dist = $sort->{$j} \n";
	    $j->Print();
	    $j->Print("",$ofhclose);
        my ($absscore,$normalizedSum) = $pdb1->ScoreGivenSetOfAtoms($pdb2,\@tmp1,\@tmp2);
        print STDERR " RSMD  : $absscore,$normalizedSum \n";
		my $LEN = 3 ;
		my ($neighbouringresidues,$neighbouringresiduereplaced) = $pdb2->NeighbouringResiduesofAtom($j,$i,$LEN);
		my $str = "";
		foreach my $r (@{$neighbouringresidues}){
			my $x =  $r->PrintSingleLetter($pdb2);
			$str = $str . $x ;
		}
		my $len = length($str);
		my $s1 = substr $str , $LEN +1  ; 
		my $s2 = reverse(substr reverse($str) , $LEN +1)  ; 
		my $xxxx = $pdb1->{THREE2SINGLE}->{$i->GetResName()};
		my $str1 = $s2 . $xxxx . $s1 ;
		print STDERR "$str - existing\n";
		print STDERR "$str1 - after mutation\n";

	    print STDERR "=============\n\n\n";
	}
	print $ofhclose "=============\n\n\n";
	
}

if(defined $moveZ){
    $pdb1->MoveAtomsAlongAxis("Z",$moveZ);
}

my $outfile1 = "$p1.rotated.pdb";
my $outfile2 = "$p2.rotated.pdb";
$pdb1->WritePDB($outfile1,1);
$pdb2->WritePDB($outfile2,1);


my $outpymol = defined $outfile ? $outfile : "$p1.$p2.p1m";
my $outpymolfh = util_write($outpymol);
util_PrintPymolWith2Proteins($outpymolfh,$outfile1,$outfile2,$atoms1,$atoms2);

util_Banner("Wrote pymol file in $outpymol");


sub parseSingleLine{
	my ($line) = @_ ; 
	my ($num,$restype,$resnum,$atom,$x,$y,$z) = split " " , $line ; 
	return ($num,$restype,$resnum,$atom,$x,$y,$z);
}
print STDERR "Output written in $outfile\n";

chmod 0777, $outfile ;
sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
sub getResultsLine{
    my ($in,$p1,$p2) = @_ ; 
    my $ifh = util_read($in);
    my @l ; 
	my ($a1,$a2);
    while(<$ifh>){
	    next if(/RESULT/);
		print $_ ;
		if(!defined $a1){
		    $a1 = $p1->ParseResultLine($_) ;
		    next ;
		}
		if(!defined $a2){
		    $a2 = $p2->ParseResultLine($_) ;
		     next;
		}
    }
    return ($a1,$a2);
}

