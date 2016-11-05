#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use MyGeom;
use MyGNM;
use PDB;
use ConfigPDB;
use Math::Geometry ;
use Math::Geometry::Planar;


use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($checkforsame,$infile,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$protein);
my ($ignorefile,@expressions,$trs);
my $howmany = 100000 ;
my $verbose = 0 ;

my $percentlength = 10;
my $percentmatched = 70;
my $percentidentity = 30;
$cutoff = 0.00000000001 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "trs=s"=>\$trs ,
            "checkforsame"=>\$checkforsame ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "verbose=i"=>\$verbose ,
            "percentlength=i"=>\$percentlength ,
            "percentmatched=i"=>\$percentmatched ,
            "percentidentity=i"=>\$percentidentity ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
$outfile = "TMP/out" if(!defined $outfile);
my $ofh = util_open_or_append($outfile);


usage( "Need to give a input file name => option -infile ") if(!defined $infile);


#my ($info,$querylength) = util_PARSEBLAST($infile);
my ($info,$querylength,$Subjectname,$queryname,$blastscore,$expect) = GNM_PARSEBLAST_BESTVAL($infile);
print "best $Subjectname,$queryname,$blastscore,$expect\n";


my @others ;
foreach my $k (@{$info}){
	my $org = $k ;
    my ($name,$description,$blastscore,$subjectlength,$iden_percent,$subjectmatched,$querystart,$queryend,$subjectstart,$subjectend,$expect) = split " ",$k ;

	my $BSintoIden = $blastscore * $iden_percent ;
	$description =~ s/ZZZ/ /g;
    print  "$queryname $name $description $blastscore $iden_percent $BSintoIden $blastscore $expect \n";
	die ;

	#my $querymatched = int(($subjectlength*$iden_percent)/100 );
	#print "$querylength $subjectlength $subjectmatched $querymatched \n";
}
die if (!defined $querylength && defined $checkforsame);

my $percentdiffofL = ($querylength * $percentlength) /100 ;
my $maxL = $querylength + $percentdiffofL ;
my $minL = $querylength - $percentdiffofL ;

my $percentfMatched = ($querylength * $percentmatched) /100 ;

my $grouping = {};
my $printed = 0 ;
foreach my $k (@{$info}){
	my $org = $k ;
    my ($name,$description,$blastscore,$subjectlength,$iden_percent,$subjectmatched,$querystart,$queryend,$subjectstart,$subjectend,$expect) = split " ",$k ;
	if($expect < $cutoff){ 
		$printed = 1 ;
		print $ofh "$name " if(! exists $grouping->{$name});
		$grouping->{$name} = 1 ;
	}
	else{
		print "not printing $name as less than $expect > $cutoff \n" if($verbose);
	}
}

my $ofhnotfound = util_open_or_append("notfoundingenome");
if($printed){
print $ofh "\n";
}
else{
}

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
