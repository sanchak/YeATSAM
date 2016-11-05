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
my ($lineidx,$mapfile,$idx,$infile,$postfix,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$protein);
my ($ignorefile,@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "mapfile=s"=>\$mapfile ,
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
            "lineidx=i"=>\$lineidx ,
            "verbose=i"=>\$verbose ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
my $ifh = util_read($infile);
usage( "Need to give a input file name => option -mapfile ") if(!defined $mapfile);
usage( "Need to give a input file name => option -lineidx ") if(!defined $lineidx);
my $ifhmap = util_read($mapfile);


print "WARNING ###### lineidx mattters \n";

## First create the maps
my $info = {};
my $infosingle = {};

while(<$ifhmap>){
     next if(/^\s*$/);
	 next if(/^\s*#/);
	 s/.*=N//; ## this may be there 

	 my @l = split ;


	 ## handle single seperately
	 if(@l eq 1){
	 	$infosingle->{$l[0]} = 1;
		next ;
	 }
	 my $main = shift @l ;

	 $info->{$main} = $main ;
	 foreach my $k (@l){
	      $info->{$k} = $main;
	 }
	 
}


my $mergetable = {};
my $CNT = 0 ;
my $CNTLINE = 0 ;
## find the ones that can be merged
while(<$ifh>){
     next if(/^\s*$/);
	 next if(/^\s*#/);
	 next if(!/^GSV/ && $lineidx eq 2);

	 $CNTLINE++;
	 my (@l) = split ; 
	 my $N = @l ;
	 my $nm = $l[0];
	 if(exists $infosingle->{$nm}){
	 	#print $ofh "XXX $nm @l \n";
		PrintIdx2Style($ofh,1,@l);
		next ;
	 }
	 
	 if(exists $info->{$nm}){
	 	 my $main = $info->{$nm};
	 	 $mergetable->{$main} = [] if(!defined $mergetable->{$main} );
		 push @{$mergetable->{$main}}, $_ ;
	 }
	 else{
	 	$CNT++;

		my @l = split ;
		if($lineidx eq 2){
		   PrintIdx2Style($ofh,1,@l);
		}
		else{

			my $A = shift @l;
			my $B = shift @l;
			my $C = shift @l;
			my @rl = reverse @l ;

			my $P = shift @rl;

			my $Q = shift @rl;
			my $R = shift @rl;
			my $S = shift @rl;
			my $T = shift @rl;
			my $U = shift @rl;
			my $V = shift @rl;

			my @rlrl = reverse @rl ;

	 	    print $ofh "1\t$A\t$B\t$C\t@rlrl\t$V\t$U\t$T\t$S\t$R\t$Q\t$P\n";
		}
	 }
}

my $NNN = (keys %{$mergetable});
print "Initially, there are $CNTLINE. There are $CNT with no maps, and $NNN with mapped\n";

foreach my $main (sort keys %{$mergetable}){
	my @l = @{$mergetable->{$main}};
	my $NMERGED = @l ;
	my $first = shift @l ;
	
	my @LL = split " ", $first ;
	if($LL[0] eq "GSVIVT01010591001"){
	    @LL = split " ", $l[0];
		#$LL[] = $l[1];
	}
	my $N = @LL ;

	    my $score = util_format_float($LL[$N-$lineidx],1);
	    #my $score = util_format_float($LL[1],1);
	    foreach my $others (@l){
	        my @LL = split " ", $others ;
	        my $N = @LL ;
	        my $otherscore = $LL[$N-$lineidx];
	        #my $otherscore = $LL[1];
		    $score = $score + $otherscore ;
	    }
	$LL[$N-$lineidx] = $score ;
	#$LL[1] = $score ;

		#my @l = split ;
		if($lineidx eq 2){
		   PrintIdx2Style($ofh,$NMERGED,@LL);
		#my @rl = reverse @LL ;
		#shift @rl ;
		#shift @rl ;
		#shift @rl ;
		#shift @rl ;
		#shift @rl ;
		#my @pl = reverse @rl ;
		#my $A = shift @pl ;
		#my $B = shift @pl ;
		#my $C = shift @pl ;
	    #$score = util_format_float($score,1);
	 	#print $ofh "$NMERGED\t$A\t$B\t@pl\t$score\n";
		}
		else{
			my @l = @LL ;
			my $A = shift @l;
			my $B = shift @l;
			my $C = shift @l;
			my @rl = reverse @l ;

			my $P = shift @rl;

			my $Q = shift @rl;
			my $R = shift @rl;
			my $S = shift @rl;
			my $T = shift @rl;
			my $U = shift @rl;
			my $V = shift @rl;

			my @rlrl = reverse @rl ;

	 	    print $ofh "$NMERGED\t$A\t$B\t$C\t@rlrl\t$V\t$U\t$T\t$S\t$R\t$Q\t$P\n";
		}
}


close($ifh);
close($ofh);


#my ($table) = util_make_list(\@list);
#my ($tablemerge,$newN) = util_table_merge($t1,$t2);
#my ($table) = util_mapFullLinetofirst($fname);
#my ($table,$N) = util_maketablefromfile($fname);
#my ($common,$inAbutnotinB,$inBbutnotinA) = util_table_diff($t1,$t2);
#my $junk = util_writelist2file($fname,@list);


sub PrintIdx2Style{
	my ($OFH,$NMERGED,@l) = @_ ;
	print $OFH "$NMERGED @l \n";
	return ;
		   my @rl = reverse @l ;
		   shift @rl ;
		   my $val = shift @rl;
		   $val = util_format_float($val,1);
		   shift @rl ;
		   shift @rl ;
		   shift @rl ;
		   my @pl = reverse @rl ;
		   my $A = shift @pl ;
		   my $B = shift @pl ;
		   my $C = shift @pl ;
		   $val = util_format_float($val,1);
	 	   print $OFH "$NMERGED\t$A\t$B\t@pl\t$val\n";
}

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
