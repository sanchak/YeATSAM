#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use MyGeom;
use MyGNM;
use MyConfigs;
use PDB;
use ConfigPDB;
use Math::Geometry ;
use Math::Geometry::Planar;
use Memory::Usage;
use Algorithm::Combinatorics qw(combinations) ;
my $mu = Memory::Usage->new();
$mu->record('');
use AAConfig;



use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($fastafile,$idx,$infile,$postfix,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$protein);
my ($ignorefile,@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "fastafile=s"=>\$fastafile ,
            "postfix=s"=>\$postfix ,
            "infile=s"=>\$infile ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "idx=i"=>\$idx ,
            "verbose=i"=>\$verbose ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
usage( "Need to give a input file name => option -fastafile ") if(!defined $fastafile);
my ($tablefasta,$infoSeq2PDB,$mapChainedName2Name) = util_parsePDBSEQRESNEW($fastafile,0,0);


my $tableallowed = {};
$tableallowed->{8} = 1 ;
$tableallowed->{9} = 1 ;
$tableallowed->{11} = 1 ;
$tableallowed->{12} = 1 ;
$tableallowed->{13} = 1 ;
$tableallowed->{14} = 1 ;
$tableallowed->{19} = 1 ;
$tableallowed->{20} = 1 ;


my $ofhCC = util_write("$outfile.CC");
my $ofhprobinmid = util_write("$outfile.probinmid");
my $ofhCdotCatend = util_write("$outfile.dotatend");
my $ofhCCatend = util_write("$outfile.CCatend");
my $ofhCdotCmorethanone = util_write("$outfile.dotmorethanone");
my $ofhNOCC = util_write("$outfile.NOCC");
my $ofhboth = util_write("$outfile.both");
my $ofhnone = util_write("$outfile.none");
my $ofhsubset = util_write("$outfile.subset.fa");
my $ofhfull = util_write("$outfile.full.fa");
my $ofhCCnotexactone = util_write("$outfile.CCnotexactone");
my $ofhPattern = util_write("$outfile.pattern");
#my $PATTERN = "CC.G.(K|R).*(D|E)(R|K)...C.C.(R|K).*CG";
my $PATTERN = "(T|S)(D|E).RY..Q";
#$PATTERN = "CC.G.(K|R).*C.C";
print "Checking for pattern =~ /$PATTERN/) \n";

my $ofhNOCatALLfa = util_write("$outfile.noCatALL.fa");
my $ofhNOCatALL = util_write("$outfile.noCatALL");

my $done = {};
foreach my $i (sort keys %{$tablefasta}){
    my $seq = $tablefasta->{$i} ;
	next if($done->{$seq});
	$done->{$seq} = 1 ;

	if($seq =~ /$PATTERN/){
		my @l = ($seq =~ /$PATTERN/g);
		print $ofhPattern ">$i\n";
		print $ofhPattern "$l[0]\n";

		print $ofhCC "$i\n";

	}
}

system("wc -l $outfile.*");

#my ($N,$first,$last,$range,$mean,$sd) = util_GetStatOfList(@l);
#my ($table) = util_make_list(\@list);
#my ($tablemerge,$newN) = util_table_merge($t1,$t2);
#my ($table) = util_mapFullLinetofirst($fname);
#my ($table,$N) = util_maketablefromfile($fname);
#my ($common,$inAbutnotinB,$inBbutnotinA) = util_table_diff($t1,$t2);
#my $junk = util_writelist2file($fname,@list);


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
