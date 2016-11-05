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



use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($fastafile,$idx,$infile,$postfix,$p1,$p2,$outfile,$cutoff,$which_tech,$lengthfile,$protein);
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
            "lengthfile=s"=>\$lengthfile ,
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
usage( "Need to give a lengthfile -option -lengthfile  ") if(!defined $lengthfile);

my $ifhLen = util_read($lengthfile);
my $lenMap = {};
while(<$ifhLen>){
     next if(/^\s*$/);
	 next if(/^\s*#/);

     chomp ;
	 my ($nm,$len) = split ; 
	 $lenMap->{$nm} = $len ;
}

while(<$ifh>){
     next if(/^\s*$/);
	 next if(/^\s*#/);
	 chomp ;
	 my $orig = $_ ;
	 $orig =~ s/=N.*/=N/;
	 s/.*=N//;

	 my (@l) = split ; 
	 foreach my $i (@l){
	 	die if(! exists $lenMap->{$i});
	 }

	 my @sl = sort { $lenMap->{$b} <=> $lenMap->{$a}} @l ;
	 foreach my $i (@sl){
	 	my $len = $lenMap->{$i} ;
		#print "$i $len\n";
	 }
	 print $ofh "$orig @sl \n";
}
close($ifh);

sub parseSingleLine{
	my ($line) = @_ ; 
	my ($num,$restype,$resnum,$atom,$x,$y,$z) = split " " , $line ; 
	return ($num,$restype,$resnum,$atom,$x,$y,$z);
}
print STDERR "Output written in $outfile\n";

chmod 0777, $outfile ;

$mu->record('after something_memory_intensive()');
$mu->dump();

#my ($N,$first,$last,$range,$mean,$sd) = util_GetStatOfList(@l);
#my ($table) = util_make_list(\@list);
#my ($tablemerge,$newN) = util_table_merge($t1,$t2);
#my ($table) = util_mapFullLinetofirst($fname);
#my ($table,$N) = util_maketablefromfile($fname);
#my ($common,$inAbutnotinB,$inBbutnotinA) = util_table_diff($t1,$t2);
#my $junk = util_writelist2file($fname,@list);


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
