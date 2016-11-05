#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use ConfigPDB;
use MyGeom;

use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($fix,$infile,$outfile,$which_tech,$listfile);
my (@expressions,$shift);
my $howmany = 100000 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "infile=s"=>\$listfile ,
            "shift"=>\$shift ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "fix=i"=>\$fix ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
if(!defined $outfile){
	$outfile = $listfile . ".tex";
}
my $ofh = util_write($outfile);
usage( "Need to give a listfile -option -listfile  ") if(!defined $listfile);
$ofh = util_write($outfile);
my $CNT = 0 ; 
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET) = util_SetEnvVars();
my $PWD = cwd;

    util_printTablePre($ofh,"caption");
my @diffs = ();

my @list= util_read_list_sentences($listfile);
my $list = {};


my $FLOATCUTTOFF = 1 ;
map { 
	s/_/\\_/g;
	my @diff = split ; 
	my $N = @diff -1 ;


	if(defined $fix){
	   foreach my $idx (1..$N){
		   my $v = $diff[$idx];
		   my $vfix = $v/$fix ;
		   $diff[$idx] = $vfix ;
	   }
	}

	shift @diff if(defined $shift);
     #print "$_ \n" ; $list->{$_} = 1 ; 
	util_printTableLine($ofh,\@diff,$FLOATCUTTOFF);
	push @diffs , \@diff ;
    } @list ;

#my ($means,$sds) = util_ProcessRowAndColumnsForMean(\@diffs,250,150);
#util_PrintMeanAndSD($ofh,$means,$sds);


util_printTablePost($ofh,"label");
print STDERR "Wrote to file $PWD/$outfile\n";


exit ; 


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
