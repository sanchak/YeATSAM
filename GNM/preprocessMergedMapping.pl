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
my ($infile,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$protein);
my ($ignorefile,@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "verbose=i"=>\$verbose ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
my $ifh = util_read($infile);

my $ignoretable = {};
if(defined $ignorefile){
   my @lll= util_read_list_sentences($ignorefile);
   map { s/\s*//g ; s/.ORF.*//; $ignoretable->{$_} = 1 ; } @lll ;
}


my $info = {};
while(<$ifh>){
chomp ;
	my @l = split ;


	my $getout = 0 ;
	foreach my $i (@l){
		if(exists $ignoretable->{$i}){
			$getout = 1;
			#print "Ignoring $i\n";
			last ;
		}
	}
	next if($getout);

	my $N = @l ;
	next if ($N eq 1);
	$info->{$_} = $N ;
}

my @sl = sort  { $info->{$b} <=> $info->{$a}} (keys %{$info});
my $CNT = 0 ;
my $ofhcnt = util_write("cntmapped");
foreach my $a (@sl){
	my $N = $info->{$a};
	$_ = $a ;
	my @l = split ;
	print $ofhcnt "$l[0] $N\n"  ;
	print $ofh "@l\n"  ;
	$CNT++;
}
print "Wrote $CNT\n";


chmod 0777, $outfile ;

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
