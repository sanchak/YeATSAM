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
my ($append,$sleep,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$string2exec);
my ($force,$ignorefile,@expressions);
my $interval ;
my $verbose = 0 ;
my $dontsleep = 0 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "string2exec=s"=>\$string2exec ,
            "append"=>\$append ,
            "sleep=i"=>\$sleep ,
            "force"=>\$force ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "interval=i"=>\$interval ,
            "dontsleep=i"=>\$dontsleep ,
            "verbose=i"=>\$verbose ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a input file name => option -sleep ") if(!defined $sleep);
usage( "Need to give a listfile -option -listfile  ") if(!defined $listfile);
usage( "Need to give a string2exec pdb id -option -string2exec  ") if(!defined $string2exec);
usage( "Need to give a interval pdb id -option -interval  ") if(!defined $interval);
#my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
 my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP,$BLASTOUT,$BLASTDB) = util_SetEnvVars();
$outfile = "Sch.$listfile" if(!defined $outfile);
my $ofh ;
if(defined $append){
     $ofh = util_append($outfile);
	 print $ofh "sleep 300 # append sleep\n";
}
else{
     $ofh = util_write($outfile);
}

my $PWD = cwd;
die " no REPLACE, need that" if ($string2exec !~ /REPLACE/);

my $ignoretable = {};
if(defined $ignorefile){
   my @lll= util_read_list_sentences($ignorefile);
   map { s/\s*//g ; $ignoretable->{$_} = 1 ; } @lll ;
}

my @list= util_read_list_words($listfile);
my $list = {};
map { s/\s*//g ; $list->{$_} = 1 ; } @list ;


my $CNT = 0 ;
my $IGN = 0 ;
my $lastitem = 0 ;
while(@list){
	my $i = shift @list;
	if(!@list){
		$lastitem = 1 ;
	}
	my $JJ = "$string2exec";
	$JJ =~ s/REPLACE/$i/g;

    my @SS = split " ", $JJ ;
	my $NN = @SS -1 ;
	my $OUTFILE = $SS[$NN];


	if((! -e $OUTFILE || -z $OUTFILE) || defined $force){
	    if($lastitem){
	         print  $ofh "sleep 10 \n" if(!$dontsleep);; ## for quick running programs 
	         print  $ofh "$JJ \n";
		}
		else{
	         print  $ofh "$JJ > & ! /dev/null &\n";
		}
	    $CNT++;
	}
	else{
		print "ignoring $i\n" if($verbose);
		$IGN++;
	}
}
close ($ofh);
print "Addd $CNT commands, ignored $IGN \n";
system("scheduleprocessInsertingsleep.pl -interval $interval -infile $outfile -sleep $sleep ");
system("head -1 Sch.$outfile");



sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
