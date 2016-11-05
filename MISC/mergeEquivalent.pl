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
my ($mergenames,$infile,$p1,$ignoreband,$cutoff,$p2,$mappingfile,$outfile,$which_tech,$ignorefile,$protein);
my (@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "mergenames"=>\$mergenames ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "ignorefile=s"=>\$ignorefile ,
            "outfile=s"=>\$outfile ,
            "mappingfile=s"=>\$mappingfile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "cutoff=i"=>\$cutoff ,
            "ignoreband=i"=>\$ignoreband ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -mappingfile ") if(!defined $mappingfile);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);


my $ofh = util_write($outfile);
my $ofhlog = util_write("merge.log");
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
my $ifh = util_read($infile);

my $ignoretable = {};
if(defined $ignorefile){
   my @lll= util_read_list_sentences($ignorefile);
   map { s/\s*//g ; $ignoretable->{$_} = 1 ; } @lll ;
}

my $ifhmap = util_read($mappingfile);
my $mappingtable = {};


while(<$ifhmap>){
	next if(/^\s*#/);
	next if(/^\s*$/);
	s/.*=N //;
	 my (@l) = split ; 
	 my $N = @l -1 ;
	 foreach my $i (0..$N){
	     my ($x) = $l[0];
		 if($x =~ /(_MERGED|_A|_B)/){
		 	my $y = shift @l ;
			push @l, $y ;
		 }
		 else{
		 	last ;
		 }
	 }
	 my ($first) = $l[0];
	 

	 $first = uc ($first);
	 #$first =~ s/.ORF.*//;

	 foreach my $i (@l){
	    #$i =~ s/.ORF.*//;
	    $mappingtable->{$i} = $first ;
	 }
}


my $NUM ; 
my $MAPPED = {};
my $CNT = 0 ;
my $DONE = {};
while(<$ifh>){
	$CNT++;
	 my (@l) = split ; 
	 my ($name) = $l[0];
	 $name = uc($name);
	 #$name =~ s/.ORF//;
	 next if(exists $ignoretable->{$name});
	 die "Done $name .......\n" if(exists $DONE->{$name});
	 $DONE->{$name} = 1 ;

	 
	 # ensure all lines have same nmber of fields
	 $NUM = @l  if(!defined $NUM);
	 die "$name not same" if(@l ne $NUM);
	

	 if(! exists $mappingtable->{$name}){
	 	## this does not have a mapping - so print as is
	 	$MAPPED->{$name} = [];
	 }
	 else{
	 	$name = $mappingtable->{$name} ;
	 	$MAPPED->{$name} = [] if(!exists $MAPPED->{$name});
	 }

	 push @{$MAPPED->{$name}}, \@l;
}


my $NN = (keys %{$MAPPED});
print "Started with $CNT, but ended with $NN, log in merge.log\n";

foreach my $k (keys %{$MAPPED}){

	my @list = @{$MAPPED->{$k}};
	my $first = shift @list ;
	my @first = @{$first};
	my $ZZZ = @first - 1 ;

	#if(@list eq 1){
	    #print $ofh "@first \n";
	    #next ;
	#}
	

	my $NEWSUM = 0 ;
	my $ORIGSUM = 0 ;
	foreach my $l (@list){
		my @ll = @{$l};
		my $N = @ll -1 ;
		die if($N ne $ZZZ);
		foreach my $idx (1..$N){
		        if(defined $mergenames){
			        #$ORIGSUM = $ORIGSUM + $first[$idx]  ;
				$first[$idx] = $first[$idx] + 0 ;
			        #$NEWSUM = $NEWSUM + $ll[$idx]  ;
			}
			else{
			$ORIGSUM = $ORIGSUM + $first[$idx]  ;
			$first[$idx] = $first[$idx] + $ll[$idx];
			$NEWSUM = $NEWSUM + $first[$idx]  ;
			}
		}
	}
	$first[0]= $k ;

	my $DIFFSUM = $NEWSUM - $ORIGSUM ;
	my $NAMEF = $first[0];
	print $ofhlog "$NAMEF $DIFFSUM\n";

	print $ofh "@first \n";
	
}

system("sort.pl -idx 1 -in merge.log -rev");


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
