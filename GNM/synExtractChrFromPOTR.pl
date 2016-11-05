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
my ($fastafile,$idx,$infile,$postfix,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$protein);
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
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
my $ifh = util_read($infile);

`mkdir -p LISTS`;
my $info = {};
my $done = {};
my $namessorted = {};
while(<$ifh>){
     next if(/^\s*$/);
	 next if(/^\s*#/);
	 chomp ;


	 if(/^>/ && /chromosome/){
	 	s/>//;
	 	my ($CDS,$A,$name) = split ;
		$CDS = uc($CDS);
		my @l = split ":", $name ;
		$name = $l[2];
		$name = $name + 0;
		my $IDX = $name ;
		$name = "chr$name";

		if(! exists $info->{$name}){
			$namessorted->{$name} = $IDX ;
			$info->{$name} = util_write("LISTS/list.$name");
			$done->{$name} = {};
		}
		


		my $FH =  $info->{$name} ;
		print $FH "$CDS\n";
	 }


}
my $ofhlist = util_write("list");
foreach my $i (sort {$namessorted->{$a} <=> $namessorted->{$b}} keys %{$namessorted}){
	print $ofhlist "$i\n";
}


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
