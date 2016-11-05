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
use MyConfigs;


use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($annotatedfile,$infile,$p1,$mapfile,$cutoff,$p2,$outfile,$which_tech,$ignorefile,$protein);
my (@expressions);
my $howmany = 100000 ;
my $verbose = 0 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "ignorefile=s"=>\$ignorefile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "annotatedfile"=>\$annotatedfile ,
            "cutoff=i"=>\$cutoff ,
            "verbose=i"=>\$verbose ,
            "mapfile=s"=>\$mapfile ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -infile ") if(!defined $infile);
usage( "Need to give a output file name => option -cutoff ") if(!defined $cutoff);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);


print "Writing to $outfile\n";
my $ofh = util_write($outfile);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
my $ifh = util_read($infile);

my $meaninfo = {};
my $sdindfo = {};

my $ignorelist = {};
if(defined $ignorefile){
    my @ignorelist= util_read_list_words($ignorefile);
    map { my @l = split ;; $ignorelist->{$l[0]} = 1 ; } @ignorelist ;
}
my $libstr = Config_getCountStr();
my @NAMES = split " " , $libstr ; shift @NAMES ;
my @NAMESSAVED = split " " , $libstr ; shift @NAMESSAVED ;
		

my $ofhpw = util_write("pairwisenames");

my @PWnames ;
my $NAMESTABLE = {};
while(@NAMES){
	my $a = shift @NAMES ;
	$NAMESTABLE->{$a} = 0;
	print "Adding $a to NAMESTABLE\n" if($verbose);
	foreach my $i (@NAMES){
		print $ofhpw "$a-$i\n";

		my @l;
		push @l, $a;
		push @l, $i;
		my @sorted = sort @l;

		push @PWnames, "$sorted[0]-$sorted[1]";
	}
}

my @NUMS ;
my $firstone = 1 ;
my $CCCC = 0 ;
while(<$ifh>){
	 my (@l) = split ; 
	 my ($name) = shift @l ;

	 if(defined $annotatedfile){
	 	shift @l ;
	 	shift @l ;
	 	my $mean = shift @l ;
		next if(defined $cutoff && $mean < $cutoff);
	 	shift @l ;
	 	shift @l ;
	 }
	 if(exists $ignorelist->{$name}){
	 	#print "Ignoring not existent $name\n";
		next ;
	 }
	 $CCCC++;
	 my $N = @l - 1;


	 foreach my $idx (0..$N){

		my $CNT = 0 ;
		## do the pairwise difference
		while(@l){
		    my $a = shift @l ;
			foreach my $i (@l){
				my $diff = abs ($a - $i);
				#print "$diff=diff $CNT\n" if($verbose);
				if($firstone){
					push @NUMS, $diff ;
				}
				else{
					$NUMS[$CNT] = $NUMS[$CNT] + $diff ;
				}
			    $CNT++;
			}
		}
	 }
	 $firstone = 0 ;
}
print "Processed $CCCC trs \n";


my $NN = @NUMS ;
die "Wrong index $NN " if($NN ne @PWnames);
$, = "\n";

## just sort first ...
my $VALUES = {};
my $max = 0 ;
foreach my $idx (1..$NN){
	my $realidx  = $idx -1 ;
    #print $ofh "$PWnames[$realidx] $NUMS[$realidx] \n"; 
	$VALUES->{$PWnames[$realidx]} = $NUMS[$realidx]  ;
	if($NUMS[$realidx] > $max){
		$max = $NUMS[$realidx] ;
	}
}

my @sorted = sort { $VALUES->{$a} <=> $VALUES->{$b}} (keys %{$VALUES});

my $GROUPCNT = 0 ;
my $GROUPS = {};

print "Max is $max ....\n";

my $NEWVALUES = {};
my $cnt = 0 ;
foreach my $k (@sorted){
	my $v = int(100*($VALUES->{$k}/$max));
	my ($a,$b) = split "-",$k ;
	print $ofh "$a $b $v\n";



}


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
