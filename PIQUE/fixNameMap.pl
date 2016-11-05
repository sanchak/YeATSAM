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
my ($fastafile,$idx,$infile,$postfix,$p1,$p2,$outfile,$cutoff,$which_tech,$mapfile,$protein);
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
            "mapfile=s"=>\$mapfile ,
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
usage( "Need to give a mapfile -option -mapfile  ") if(!defined $mapfile);

my ($tableNames) = util_maketablefromfile($mapfile);

while(<$ifh>){
	if(/===/){
		my ($e,$p1,$p2,$NM1,$NM2,@l);
		if(/EDGE/){
		    ($e,$p1,$p2,$NM1,$NM2,@l) = split ;
		}
		else{
		    ($NM1,$NM2,@l) = split ;
		}
		if(exists $tableNames->{$NM1}){
             my $newnm1 = $tableNames->{$NM1} ;
			 $NM1 = $newnm1 ;
		}
		if(exists $tableNames->{$NM2}){
             my $newnm2 = $tableNames->{$NM2} ;
			 $NM2 = $newnm2 ;
		}
		if(/EDGE/){
		    print $ofh "$e $p1 $p2 $NM1 $NM2 @l\n";
		}
		else{
		    print $ofh " $NM1 $NM2 @l\n";
		}
	}
	else{
		print $ofh $_ ;
	}
}

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
