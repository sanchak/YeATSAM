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
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);

#usage( "Need to give a listfile -option -listfile  ") if(!defined $listfile);
#usage( "Need to give a protein pdb id -option -protein  ") if(!defined $protein);
my $CNT = 0 ; 
#my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT) = util_SetEnvVars();
my $PWD = cwd;

my ($infofpocket,$infoperatom,$volume) = parseFPocketFile($fpocket);
util_table_print($infofpocket);

#usage( "Need to give a input file name => option -castp ") if(!defined $castp);
#my ($infocastp,$infobiggest) = parseCastp($castp,$subtract);
#util_table_print($infobiggest);

print "Vilume fpocket = $volume\n";


if(0){
### cnt number of matches - not working for now, since I have hashed out $CNT++
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
my $ifh = util_read($infile);
while(<$ifh>){
	chomp ;
	my $ORIG = $_ ;
	s/@.*//;
	my @l = split ;
	my $CNT = 0 ;
	foreach my $l (@l){
		my (@ll) = split "/", $l ;
		#my $k = $ll[0] . "/" . $ll[1];
		my $k = $ll[1];
		#$CNT++ if(exists $infobiggest->{$k});
	}
	if($CNT >= 4){
		print "$ORIG $CNT \n";
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
