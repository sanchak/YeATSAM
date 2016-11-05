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
my $ofh = util_write($outfile);
my $ofhlist = util_write("$outfile.list");
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
my $ifh = util_read($infile);

my $info = {};



my $done =  {};
while(<$ifh>){
     next if(/^\s*$/);
	 next if(/^\s*#/);
	 next if(/immobile/);
	 next if(/derived/);
	 next if(/verified/);

	 #my @l = split "YYYYYY", $_  ;
	 my @l = split "\t", $_  ;

	 s/.*genomes.all.//;
	 my ($a) = split ;
	 $a =~ s/YYYYYY//;

	 next if($a =~ /many\s*$/);
	 next if($a =~ /high\s*$/);
	 next if($a =~ /partial\s*$/);
	 next if($a =~ /John$/);
	 next if($a =~ /Seqman$/);
	 next if($a =~ /genome$/i);
	 next if($a =~ /missing/i);

	 #$a =~ s/_/AAA/;
	 #$a =~ s/_/BBB/;
	 #$a =~ s/AAA/_/;
	 #my @twoterms = split "BBB", $a ;



	 ## for Brevibacterium 
	 my @ll = split " ",$l[7] ;
	 $l[7] =~ s/\[/_/;
	 $l[7] =~ s/\]/_/;

	 ## remove _ in the real name
	 my $ORIG = $ll[0];
	 $ORIG =~ s/\[//;
	 $ORIG =~ s/\]//;
	 next if(exists $done->{$ORIG});



	 my $STR = "ftp://ftp.ncbi.nih.gov/genomes/refseq/bacteria/${ll[0]}_${ll[1]}/latest_assembly_versions/$a/${a}_genomic.fna.gz";
	 $STR =~ s/GCA/GCF/g;

	 my $FINAL = "wget  --spider -v $STR \n";

	 my $exit_status   = system("$FINAL");

	 next if($exit_status);

	 print $ofh "wget $STR\n";
	 $done->{$ORIG} = 1 ;
	 print $ofhlist "${a}_genomic.fna.gz\n";

}


