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
my ($f1,$idx,$infile,$f2,$postfix,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$blastdir);
my ($ignorefile,@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "blastdir=s"=>\$blastdir ,
            "f1=s"=>\$f1 ,
            "f2=s"=>\$f2 ,
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
my $ofh = util_write($outfile);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
my $ifh = util_read($infile);
usage( "Need to give a input file name => option -f1 ") if(!defined $f1);
usage( "Need to give a input file name => option -f2 ") if(!defined $f2);

usage( "Need to give a blastdir pdb id -option -blastdir  ") if(!defined $blastdir);
my $info = {};
system("mkdir -p $blastdir");
while(<$ifh>){
     next if(/^\s*$/);
	 next if(/^\s*#/);
	 my ($A,$B) = split ; 
	 print $ofh "myblastcompare2Fastafiles.csh $f1/$A.ALL.1.fasta $f2/$B.ALL.1.fasta $blastdir/$A.blast.nt P\n";
}
system ("wc -l $outfile");
close($ifh);


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
