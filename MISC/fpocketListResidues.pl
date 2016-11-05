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


use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($infile,$outfile,$which_tech,$fpocket,$listfile,$protein);
my ($castp,@expressions);
my $subtract = 0 ;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "fpocket=s"=>\$fpocket ,
            "castp=s"=>\$castp ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "subtract=i"=>\$subtract ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $fpocket);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
usage( "Need to give a output file name => option -protein ") if(!defined $protein);
my $ofh = util_write($outfile);

#usage( "Need to give a input file name => option -castp ") if(!defined $castp);
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT) = util_SetEnvVars();
my $PWD = cwd;
my $pdb = "$PDBDIR/$protein.pdb";
my $pdb1 = new PDB();
$pdb1->ReadPDB($pdb);

my ($infofpocket,$infoperatom,$volume) = parseFPocketFile($protein,$pdb1,0,$fpocket);

#my ($infocastp,$infobiggest) = parseCastp($castp,$subtract);

my $colorcoding = {};
$colorcoding->{"ARG"} = "blue";
$colorcoding->{"LYS"} = "blue";
$colorcoding->{"HIS"} = "blue";
$colorcoding->{"ASP"} = "red";
$colorcoding->{"GLU"} = "red";

my $CNTres = 0 ;
my $ID = "PROT1";
my $sortprint = {};
foreach my $k (keys %{$infofpocket}){
	my $v = $infofpocket->{$k};
	$v =~ s/\// /g;
	my ($name)  = split " " ,$v ;
    my $color = "orange" ;
    my $shape = "sticks" ;
	if(exists $colorcoding->{$name}){
		$color = $colorcoding->{$name} ;
		$shape = "spheres";
	}
	$sortprint->{$k} = $name ;
	$CNTres++ ;

	print $ofh  "select block_query$CNTres, /$ID//A/$k\n";
	print $ofh  "color $color, block_query$CNTres\n";
	print $ofh  "show $shape, block_query$CNTres\n";
}


foreach my $k (sort {$a <=> $b} keys %{$sortprint}){
	my $v = $sortprint->{$k};
	print " $v/$k ";
}
print "\n";

print "Volume fpocket = $volume\n";



sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
