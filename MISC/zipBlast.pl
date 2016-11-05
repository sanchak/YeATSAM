#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use MyGeom;
use PDB;
use ConfigPDB;
use Math::Geometry ;
use Math::Geometry::Planar;


use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($blastdir,$idx,$infile,$pairwise,$p2,$outfile,$cutoff,$which_tech,$listfile,$protein);
my ($ignorefile,@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "blastdir=s"=>\$blastdir ,
            "pairwise=i"=>\$pairwise ,
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
usage( "Need to give a input file name => option -blastdir ") if(!defined $blastdir);
usage( "Need to give a input file name => option -pairwise ") if(!defined $pairwise);

my $postfix = ".blast";
system ("cd $blastdir ; unlink LLL ; ls | grep $postfix > LLL ");
system ("echo there are ; wc -l $blastdir/LLL ");

my @list= util_read_list_sentences("$blastdir/LLL");
my $list = {};
map { s/\s*//g ; $list->{$_} = 1 ; } @list ;




my $outdir = "$blastdir". ".zblast.$pairwise" ;
system ("mkdir -p $outdir");

my $CNT = 0 ;
foreach my $infile (@list){
    next if (-e "$outdir/$infile");
    $CNT++;

    my $ifh = util_read("$blastdir/$infile");
    my $ofh = util_write("$outdir/$infile");


    my $scoreseen = 0 ;
    my $startprinting = 0 ;
    my $doneprinting = 0 ;
    my $cntinfo = 0 ;
    while(<$ifh>){
        next if(/^\s*$/);
	next if(/^\s*#/);
	next if(/^BLAST/);

	if($pairwise){

	     if(/Strand=Plus/){
	         $scoreseen = 1 ;
	         next ;
	     }

	     if($scoreseen && /\|/){
		     print $ofh "\n";
		     next ;
	     }
	     last if(/^Lambda/);

             print $ofh $_ ;
	}
	else{
		if(/^>/){
			$cntinfo++;
		}
		if($cntinfo>1){
			last ;
		}
		if(/^Query/){
		    $startprinting = 1;
		}
		if($startprinting){
                    print $ofh $_ ;
		}

	}
    }
    close($ofh);
    close($ifh);
    #die "$outdir/$infile";
}
print "Wrote $CNT number of files \n";

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
