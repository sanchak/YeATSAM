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

## FINFO
## Compares two files - with different blastscores....
## Prints the common ones 
## Not in each are not handled - todo???


use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($scores1,$idx,$scores2,$postfix,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$protein);
my ($ignorefile,@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "scores1=s"=>\$scores1 ,
            "postfix=s"=>\$postfix ,
            "scores2=s"=>\$scores2 ,
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
my $ofhfinal = util_write("$outfile.final");
usage( "Need to give a input file name => option -scores2 ") if(!defined $scores2);
my $ifh1 = util_read($scores2);
usage( "Need to give a input file name => option -scores1 ") if(!defined $scores1);
my $ifh2 = util_read($scores1);


my ($info1,$anno1) = ParseScoreFile($scores1);
my ($info2,$anno2) = ParseScoreFile($scores2);

my $common = {};
my $notinB = {};
my $both = {};
foreach my $k (sort keys %{$info1}){
	 my $s1 = $info1->{$k};
	 if(exists $info2->{$k}){
	     my $s1 = $info1->{$k};
	     my $s2 = $info2->{$k};
		 my $percent = util_format_float(100*(($s1 - $s2)/$s1));
		 next if(defined $cutoff && $percent < $cutoff);
		 my $diff = util_format_float(($s1 - $s2),1) ;
		 my $k = "$k $s1 $s2 $diff $percent";
		 $common->{$k} = $percent ;
		 $both->{$k} = $s1 ;
	 }
	 else{
	 	$notinB->{$k} = $s1 ;
	 	$both->{$k} = $s1 ;
	 }
}

my @sort = sort { $common->{$b} <=> $common->{$a}} (keys %{$common});
foreach my $k (@sort){
	print $ofh "$k\n";
}

my @sortboth = sort { $both->{$b} <=> $both->{$a}} (keys %{$both});
foreach my $k (@sortboth){
	my $ANNO = $anno1->{$k} ;
	print $ofhfinal "$ANNO";
}

close($ofh);
system("wc -l $outfile*");
system("extractindexfromfile.pl -in $outfile.final");

sub ParseScoreFile{
   my ($fname) = @_ ;
   my $ifh = util_read($fname);
   my $info = {};
   my $anno = {};
   while(<$ifh>){
     next if(/^\s*$/);
	 next if(/^\s*#/);
	 my $ORIG = $_ ;

     chomp ;
	 my (@l) = split ; 
	 my $N = @l; 
	 my $nm = $l[0];
	 my $score = $l[$N-2];
	$info->{$nm} = $score ;
	 $anno->{$nm} = $ORIG ;
   }
   close($ifh);
   return ($info,$anno) ;
}


sub parseSingleLine{
	my ($line) = @_ ; 
	my ($num,$restype,$resnum,$atom,$x,$y,$z) = split " " , $line ; 
	return ($num,$restype,$resnum,$atom,$x,$y,$z);
}
print STDERR "Output written in $outfile\n";

chmod 0777, $outfile ;

$mu->record('after something_memory_intensive()');
$mu->dump();


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
