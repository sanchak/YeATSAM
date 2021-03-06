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
my ($reverse,$infile,$p1,$p2,$outfile,$ignoresingle,$cutoff,$which_tech,$listfile,$protein);
my ($ignorefile,@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "p1=s"=>\$p1 ,
            "ignoresingle"=>\$ignoresingle ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "reverse"=>\$reverse ,
            "verbose=i"=>\$verbose ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
my $ifh = util_read($infile);
$outfile = "$infile.mapped" if(!defined $outfile);
my $ofh = util_write($outfile);
my $ofhpw = util_write("$outfile.pw");
my $ofhsingle = util_write("$outfile.single");



my $info = {};
while(<$ifh>){
	my @l = split ;
	my $N = @l ;
	my $first = shift @l ;
	if($N eq 1){
		 print $ofhsingle "$first\n";
	}
	next if(defined $ignoresingle && $N eq 1);
	$info->{$_} = $N ;
	foreach my $i (@l){
		print $ofhpw "$first $i 0.1\n";
	}

}

my @sl ;
if(defined $reverse){
    @sl = sort  { $info->{$a} <=> $info->{$b}} (keys %{$info});
}
else{
    @sl = sort  { $info->{$b} <=> $info->{$a}} (keys %{$info});
}

my $ofhcount = util_write("$outfile.count");
foreach my $a (@sl){
	my $count = $info->{$a};
	print $ofh $a ;
	my @l = split " ",$a ;
	print $ofhcount "$l[0] $count \n";
}


chmod 0777, $outfile ;

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
