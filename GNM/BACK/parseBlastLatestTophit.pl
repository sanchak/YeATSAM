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
my ($checkforsame,$infile,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$protein);
my ($ignorefile,@expressions,$trs);
my $howmany = 100000 ;
my $verbose = 1 ;

my $percentlength = 10;
my $percentmatched = 70;
my $percentidentity = 30;
$cutoff = 1E-10;
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
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);


usage( "Need to give a input file name => option -infile ") if(!defined $infile);
usage( "Need to give a input file name => option -trs ") if(!defined $trs);


my ($info,$querylength) = util_PARSEBLAST($infile);





my $sorted = {};
foreach my $k (@{$info}){
    my ($name,$description,$blastscore,$subjectlength,$iden_percent,$subjectmatched,$querystart,$queryend,$subjectstart,$subjectend,$expect) = split " ",$k ;
    $sorted->{$k} = $blastscore ;
}

my @sorted = sort { $sorted->{$b} <=> $sorted->{$a} } (keys %{$sorted});
my $NNN = @sorted ;
if(!@sorted){
     print $ofh "$trs Warning:Did not find any characterized\n";
}

my $UNCHARSTRINGS = Config_GetUncharStrings();

my $first ; 
my $percent ; 
foreach my $k (@sorted){
    my $org = $k ;
    my ($name,$description,$blastscore,$subjectlength,$iden_percent,$subjectmatched,$querystart,$queryend,$subjectstart,$subjectend,$expect) = split " ",$k ;

    ## choose the first percent
    $description =~ s/ZZZ/ /g;
    my $str2print ;
    if(!defined $first){
    	if($expect > $cutoff){
           print $ofh "$trs Warning:Did not find any characterized\n";
	    }
        $percent = int(100 * ($querylength/$subjectlength));
        $str2print = "$trs\t$name\t$description $percent $blastscore $expect\n";
        $first = $str2print ;
    }

    $str2print = "$trs\t$name\t$description $percent $blastscore $expect\n";

    if($expect > $cutoff){
       print $ofh "$first\n";
    }
    next if($description =~ /$UNCHARSTRINGS/);

    print $ofh "$str2print";
    last ;
}

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
