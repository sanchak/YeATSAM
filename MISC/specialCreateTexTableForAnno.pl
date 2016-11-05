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
my ($idx,$infile,$postfix,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$protein);
my ($ignorefile,@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
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
my $ifh = util_read($infile);
$outfile = "$infile.tex" if(!defined $outfile);
my $ofh = util_write($outfile);
util_printTablePre($ofh,"$infile");
print $ofh  "\\hline\n";
print $ofh  " TRID & Accid & Description & BBS\\\\ \n";
print $ofh  "\\hline\n";
while(<$ifh>){
     next if(/^\s*$/);
	 next if(/^\s*#/);

	 s/_/\\_/g;
	 my @l = split ;
	 
	 my $nm = shift @l ;
	 my $id = shift @l ;

	 my @RL = reverse @l ;
	 my $eval = shift @RL ;
	 my $bbs = shift @RL ;
	 my $percent = shift @RL ;
	 my $ISNT = shift @RL ;
	 my @anno = reverse @RL ;

	 my $A = lc(shift @anno) ;

	 my ($A1,$A2) = ($A =~ /(.)(.*)/);
	 $A1 = uc($A1);
	 my $B = lc(shift @anno) ;

	 my $ANNSTR = join " ", @anno;
	 $ANNSTR = lc($ANNSTR);


	 print $ofh  "$nm & $id & {\\it ${A1}$A2 $B} $ANNSTR & $bbs \\\\ \n";

}
print $ofh  "\\hline\n";
util_printTablePost($ofh,"$infile");
chmod 0777, $outfile ;

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
