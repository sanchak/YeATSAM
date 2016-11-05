#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use PDB;
use ConfigPDB;
use MyGeom;
# just test

use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($infile,$outfile,$which_tech,$listfile,$protein);
my (@expressions,$config);
my $howmany ;
my $verbose = 1 ;
my $LARGESTDISTALLOWED = 45 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "config=s"=>\$config,
            "howmany=i"=>\$howmany ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
usage( "Need to give a config file name => option -config ") if(!defined $config);
usage( "Need to give a listfile -option -listfile  ") if(!defined $listfile);
#usage( "Need to give a howmany -option -howmany  ") if(!defined $howmany);
usage( "Need to give a protein pdb id -option -protein  ") if(!defined $protein);
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
my $PWD = cwd;

my  ($seconds, $microseconds) = gettimeofday;

my $ofhresultstyle = util_write("$protein.resultstyle");

my $file1 = "$PDBDIR/$protein.pdb";
my $pdb1 = new PDB();
$pdb1->ReadPDB($file1);
ConfigPDB_Init($config);

my @list= util_read_list_sentences($listfile);
my $list = {};
my @reallist ; 
#map {my @j = split ; push @reallist, @j;  } @list ;
map { s/,/ /g ; my @j = split ; push @reallist, @j;  } @list ;
map { s/\s*//g ; $list->{$_} = 1 ; } @reallist ;

#my $size = @reallist ;
my @chain;
my @atoms ; 
my $done ;
my $original ;
my $CNT = 0 ; 
foreach my $i (@reallist){
	print "$i\n";
	$CNT++ ;
	$i =~ s/,//g;
	$i =~ s/-//g;

    my ($name,$number,$type1) ;
	if($i =~ /\//){
	   ($name,$number,$type1) = split "/", $i ;
	}
	else{
	   ($name,$number) = ($i =~ /([a-zA-Z]+)([0-9]+)/);
	}
	$name = uc($name);
	#print "$name,$number \n";
	my $len = length($name); 
	die "Wrong length $len" if($len != 1 && $len != 3);
	if($len == 1){
	      $name = $pdb1->GetThreeLetter($name);
	}
	#print "$name,$number \n";

	my ($res) = $pdb1->GetResidueIdx($number);
	my $type = ConfigPDB_GetAtom($res->GetName()) or die;
	print $ofhresultstyle  "$name/$number/$type ";
	print   "$name/$number/$type ";
	my ($atom1) = $pdb1->GetAtomFromResidueAndType($number,$type);
	$atom1->Print();
	push @atoms, $atom1 ;
	$done->{$number} = 1 ; 
	$original->{$number} = 1 ; 


	push @chain,$name;
	push @chain,$number;
	if(defined $howmany && $howmany && $CNT eq $howmany){
		last ; 
	}
}
print $ofhresultstyle  "\n";

my $DIFF = $howmany - $CNT ;


if($howmany && $DIFF > 0){
    my $dist = 1; 
	my $CCC = 0 ; 
    while($dist < 7 && $DIFF > 0){
	    print "Need $DIFF more atoms for dist $dist\n";
        foreach my $atom1 (@atoms){
	        my $list = util_make_list($atom1);
	        my ($junk,$neigh)  = $pdb1->GetNeighbourHoodAtom($list,$dist);
        
			my $addedcnt = 0 ;
	        foreach my $a (@{$neigh}){
	            my $number = $a->GetResNum();
				if(!exists $done->{$number} && $DIFF > 0 ){

	                my ($res) = $pdb1->GetResidueIdx($number);
				    my $name = $res->GetName();
	                my $type = ConfigPDB_GetAtom($res->GetName()) ;

    
					if(defined $type){
	                   my ($atom1) = $pdb1->GetAtomFromResidueAndType($number,$type);
					   $DIFF-- ;
	                   push @atoms, $atom1 ;
	                   $done->{$number} = 1 ; 
	                    push @chain,$name;
	                    push @chain,$number;
				       $addedcnt++;
					 }
				}
	        }
			print "Added $addedcnt for dist $dist\n";
        }
	    $dist++ ; 
    }
}

my @distlist = sort @{$pdb1->DistanceInGivenSetOfAtoms(\@atoms)};
my $N = @distlist - 1 ;
my $largestdist = $distlist[$N];
if($largestdist > $LARGESTDISTALLOWED){
	$, = " , " ;
	print "Distances are @distlist \n";
    die "Error: $protein: Largest distance = $largestdist  more than $LARGESTDISTALLOWED"  ;
}


my $size = @chain /2 ;
ConfigPDB_PrintOutConf($outfile,$protein,$size,\@chain,1);

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
