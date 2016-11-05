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
my ($fastafile,$idx,$infile,$postfix,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$protein);
my ($ignorefile,@expressions);
my $howmany = 100000 ;
my $verbose = 0 ;
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
#$cutoff = 1E-20 ;
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
usage( "Need to give a input file name => option -cutoff ") if(!defined $cutoff);
my $ifh = util_read($infile);

my $info = {};

##############################################################
### First get the last PDB - ie the first with blastscore > $cutoff
### We need to do this as we have to see each entry within the ">" lines that come afterwards
##############################################################
my $LASTNM ;
while(<$ifh>){
     next if(/^\s*$/);
	 next if(/^\s*#/);
	 last if(/^>/);

	 next if(!/^pdb/);
	 my (@l) = split ; 
	 my $N = @l ;
	 my $nm = $l[0];
	 my $blastscore = $l[$N - 2];
	 my $BREAK = 0 ;
	 if($blastscore < $cutoff){
		$LASTNM = $nm ;
		while(<$ifh>){
			if(/^\s*$/){
				$BREAK = 1 ;
				last ;
			}
		}

	 }
	 last if($BREAK);
}


($LASTNM =~ s/\|/ /g);
my @ll = split " ", $LASTNM ;
$LASTNM = $ll[1] ;

my $NMS = {};
my $NMSWITHCHAIN = {};
print "LASTNM = $LASTNM for infile $infile \n";
while(<$ifh>){
	last if(/$LASTNM/);
	s/>//;
	if(/^\s*pdb/){
        s/\|/ /g;
        my ($junk,$PDB,$CHAIN) = split ;
		my $PDBWITHCHAIN = $PDB . $CHAIN ;
		print " $PDB $CHAIN\n" if($verbose);
		$NMS->{$PDB} = 1 ;
		$NMSWITHCHAIN->{$PDBWITHCHAIN} = 1 ;
	}
}
close($ifh);


my $ofh = util_write($outfile);
foreach my $k (sort keys %{$NMS}){
	print $ofh "$k\n";
}

my $ofhwithchain = util_write("$outfile.withchain");
foreach my $k (sort keys %{$NMSWITHCHAIN}){
	print $ofhwithchain "$k\n";
}

system("wc -l $outfile");



sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
