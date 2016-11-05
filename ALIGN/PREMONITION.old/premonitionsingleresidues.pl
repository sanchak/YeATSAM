#!/usr/bin/perl -w 
use strict ;
use PDB;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use Algorithm::Combinatorics qw(combinations) ;
use MyUtils;
use ConfigPDB;
use MyGeom;
use MyPymol;
use Math::Geometry ;
use Math::VectorReal qw(:all);  # Include O X Y Z axis constant vectors




use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($mlength,$p1,$p2,$infile,$outfile,$dist,$which_tech,$listfile,$protein);
my (@expressions);
my $maxdist = 3 ;
my $verbose = 1 ;
my ($moveZ,$verify,$decaaf,$before1,$resnum,$before2);
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
            "resnum=i"=>\$resnum ,
            "mlength=i"=>\$mlength ,
            "moveZ=i"=>\$moveZ ,
            "dist=f"=>\$dist ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_append($outfile);
my $ofhclose = util_write("log.close");
usage( "Need to give a protein 1 id -option -p1  ") if(!defined $p1);
usage( "Need to give a dist -option -dist  ") if(!defined $dist);
usage( "Need to give a mlength -option -mlength  ") if(!defined $mlength);
usage( "Need to give a resnum -option -resnum  ") if(!defined $resnum);
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
my $PWD = cwd;


my $file1 = "$PDBDIR/$p1.pdb";
my $pdb1 = new PDB();
$pdb1->ReadPDB($file1);

my @res = $pdb1->GetResidues();
my $N = @res;
my $sum = 0 ;
my $cntmatch = 0 ; 

my $DB = {};
my $done = {};
my ($res) = $pdb1->GetResidueIdx($resnum);
{
    return  if($res->GetAtomStr() ne "ATOM");
    my $CAatom = $pdb1->GetAtomFromResidueAndType($resnum,"CA");

	#### Find residues in a distance $dist from the CA atom
	    my @atomlist ;
	    push @atomlist, $CAatom ;
	    my ($results,$combined) = $pdb1->GetNeighbourHoodAtom(\@atomlist,$dist);
	    my $sort ;
		my @goodatoms ; 
		my $residues = {};
		my @Rlist ;
        foreach my $j (@{$combined}){
		    my $atomstr = $j->GetAtomStr();
			next if($atomstr eq "HETATM");
		    my $resnum = $j->GetResNum(); 
			if(!exists $residues->{$resnum}){
			    $residues->{$resnum} = 1 ;
			    my ($res) = $pdb1->GetResidueIdx($resnum);
				push @Rlist, $res ; 
			}
		}

		my $N = @Rlist ; 
		print "$resnum has $N resodues \n";
		print "Getting combined\n";
		my $iter = combinations(\@Rlist, $mlength);
		while (my $c = $iter->next) {
		    my @combo = @{$c} ;

			## sort based on single letter
		    my @resultssorted = sort { $a->PrintSingleLetter($pdb1) gt $b->PrintSingleLetter($pdb1) } @combo ;

			## create the index string for the amino acid string
		    my $nm = "";
		    my $index = "";
			foreach my $r (@resultssorted){
				my $n = $r->GetResNum();
				my $type = $r->PrintSingleLetter($pdb1);
				$nm = $nm . $type ;
				$index = $index . "." . $n ; 
			}

			## insert in the DB for that particular amino acid string
			if(!exists $done->{$index}){
				$done->{$index} = 1 ; 
				if(!exists $DB->{$nm}){
				    $DB->{$nm} = [];	
				}
				push @{$DB->{$nm}}, $index ; 
			}
		}
}

print $ofh "# dist=$dist mlength =$mlength \n";
foreach my $k (keys %{$DB}){
	my @l = @{$DB->{$k}};
    print $ofh "$k ";
	foreach my $i (@l){
		print $ofh "$i ";
	}
    print $ofh "\n";
}


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}

