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
my ($mapfile,$idx,$infile,$postfix,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$protein);
my ($ignorefile,@expressions);
my $howmany = 100000 ;
my $verbose = 0 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "mapfile=s"=>\$mapfile ,
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
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
$outfile = "$infile.first.len.sort";
my $ofh = util_write($outfile);
my $ofhanoo = util_write("$outfile.anno");
my $ofhfirst = util_write("$outfile.first");
my $ifh = util_read($infile);
usage( "Need to give a input file name => option -mapfile ") if(!defined $mapfile);

#my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
 my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP,$BLASTOUT,$BLASTDB) = util_SetEnvVars();
my $aaconfig = new AAConfig("$SRC/aa.config");

# for reading fasta
#my ($fe,$x) =  util_readfasta($infile);
#chomp $x ;
my $PWD = cwd;


my ($maptable) = util_maketablefromfile($mapfile);

my $printsort = {};
my $LISTS = {};
while(<$ifh>){
     next if(/^\s*$/);
	 next if(/^\s*#/);
	 s/.*=N//;
	 my @l = split ;
	 my $N = @l ;
	 if($N eq 1 ){
	 	print $ofh $_ ;
	 }
	 else{
         my $info = {};
		 my @origlens ;
	     foreach my $i (@l){
	 	     die "$i does not exist in $mapfile" if(! exists $maptable->{$i});
		     my $len = $maptable->{$i};
			 push @origlens, $len;
			 $info->{$i} = $len ;
	     }
		 my @resultssorted = sort { $info->{$b} <=> $info->{$a} } @l ;
		 my @newlens ;
	     foreach my $i (@resultssorted){
		     my $len = $maptable->{$i};
			 push @newlens, $len;
		  }

		 if($verbose){
		    print "orig @l @origlens\n";
		    print "new @resultssorted @newlens\n";
		    print "\n\n";
		 }
		 print $ofh "@resultssorted\n";
		 my $NN = @resultssorted ;
		 my $first = $resultssorted[0];
		 $printsort->{$first} = $NN ;
		 $LISTS->{$first} = \@resultssorted ;
	}
}

my @SSS = sort { $printsort->{$b} <=> $printsort->{$a} } (keys %{$printsort});
foreach my $first (@SSS){
	my @l = @{$LISTS->{$first}};
	print $ofh "@l \n";
	print $ofhfirst "$first\n";
	foreach my $i (@l){
		 my $len = $maptable->{$i};
		 print $ofhanoo "$i $len   " ;
	}
	print $ofhanoo "\n";
}

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
