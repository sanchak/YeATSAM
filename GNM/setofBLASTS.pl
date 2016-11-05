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
my ($blastout,$fastafile,$idx,$infile,$dbdir,$postfix,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$dblist);
my ($ignorefile,@expressions);
my $howmany = 4 ;
my $verbose = 0 ;
my $direction = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "blastout=s"=>\$blastout ,
            "dblist=s"=>\$dblist ,
            "dbdir=s"=>\$dbdir ,
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
            "direction=i"=>\$direction ,
            "idx=i"=>\$idx ,
            "verbose=i"=>\$verbose ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_open_or_append($outfile);
usage( "Need to give a input file name => option -dbdir ") if(!defined $dbdir);
usage( "Need to give a listfile -option -listfile  ") if(!defined $listfile);
usage( "Need to give a dblist pdb id -option -dblist  ") if(!defined $dblist);
usage( "Need to give a blastout pdb id -option -blastout  ") if(!defined $blastout);

my (@lll) = util_read_list_firstitem($listfile);

my @dblist= util_read_list_sentences($dblist);
my ($tableTRS) = util_maketablefromfile($outfile);

foreach my $i (@lll){
	next if(exists $tableTRS->{$i}) ;

	my $infasta = "FASTADIR/$i.ALL.1.fasta";
	my $CNT = 0 ;
    foreach my $j (@dblist){
		my $OUT = "$blastout/$i.$j.blast.nt";
	    if(! -e $OUT){
		   #print "BLASTN $dbdir/$j $infasta  $OUT \n";
   	       system("BLASTN $dbdir/$j $infasta  $OUT");
		}
		my $ifh = util_read($OUT);
	    my $hit = 0 ;
		while(<$ifh>){
			if(/No hit/){
			   last ;
			}
			elsif(/Expect/){
				$hit = 1;
				last ;
			}
		}


		if($direction){
		    $CNT++ if($hit);
		    if($CNT >= $howmany){
			    print "$i has two matches : last was $j ...so last \n" if($verbose);
			    last ;
		    }
		}
		else{
		    $CNT++ if(!$hit);
		    if($CNT > $howmany){
			    last ;
		    }
		}
	}
	print $ofh "$i $CNT $howmany\n";
}


chmod 0777, $outfile ;

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
