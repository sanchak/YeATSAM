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
use Memory::Usage;
use Algorithm::Combinatorics qw(combinations) ;
my $mu = Memory::Usage->new();
$mu->record('');
use AAConfig;



use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($fastadir,$idx,$infile,$postfix,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$protein);
my ($ignorefile,@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "fastadir=s"=>\$fastadir ,
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
#usage( "Need to give a input file name => option -infile ") if(!defined $infile);
#my $ifh = util_read($infile);
usage( "Need to give a input file name => option -fastadir ") if(!defined $fastadir);

usage( "Need to give a listfile -option -listfile  ") if(!defined $listfile);
my @list= util_read_list_sentences($listfile);

## Simple program to get three ORFs from each trs
my @files;
foreach my $i (@list){
	$i =~ s/ //g;
	my $listf = "$fastadir/$i.list";
	die  " $listf is zero or not there " if(! -e $listf || -z $listf) ;
	my $ifh = util_read($listf);
	while(<$ifh>){
		my ($x) = split ;
	    my $ffile = "$fastadir/$x.ALL.1.fasta";
	    die " $ffile is zero or not there "  if(! -e $ffile || -z $ffile);
		push @files, $x ;
	}
	close($ifh);
}

my $ofh = util_write($outfile);
foreach my $f (@files){
	 print $ofh "$f \n";
}


chmod 0777, $outfile ;


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
