#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use Carp ;
use MyUtils;
use MyGeom;
use PDB;
use ConfigPDB;
use Math::Geometry ;
use Math::Geometry::Planar;


use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($ignorefile,$direction,$infile,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$protein);
my (@expressions);
my $howmany = 100000 ;
my $verbose = 0 ;

$direction = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "p1=s"=>\$p1 ,
            "ignorefile=s"=>\$ignorefile ,
            "p2=s"=>\$p2 ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "direction=i"=>\$direction ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recog:nize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
usage( "Need to give a output file name => option -cutoff ") if(!defined $cutoff);
my $ofh = util_write($outfile);
my $ofhfirst = util_write("$outfile.first");
my $ofhallingroup = util_write("$outfile.allingroup");

usage( "Need to give a input file name => option -infile ") if(!defined $infile);
my $ifh = util_read($infile);

my $ignorelist = {};
if(defined $ignorefile){
my @ignorelist= util_read_list_sentences($ignorefile);
map { s/\s*//g ; $ignorelist->{($_)} = 1 ; } @ignorelist ;
}

my $info = {};
while(<$ifh>){
	next if(/^\s*#/);
    s/#.*//;
	my @l = split ;
	croak  "$_  in $infile was expecting  3 " if(@l ne 3);
	 my ($a,$b,$value) = split ; 
	 next if(exists $ignorelist->{$a} || exists $ignorelist->{$b});
	 next if ($a eq $b);
	 my $bool = $direction ? $value > $cutoff : $value < $cutoff ;
	 if($bool){
	     $info->{$a} = {} if(!defined $info->{$a});
	     $info->{$b} = {} if(!defined $info->{$b});
	     $info->{$a}->{$b} = $value ;
	     $info->{$b}->{$a} = $value  ;
	 }
}

my $SETS = {};
my $SETN = 0 ;
my $TOTALSETS = 0 ;

my @ALL = (sort keys %{$info});
my $N = @ALL ;
print "There were $N groups to start with\n";

my $element2Group = {};
foreach my $element (@ALL){
#print "$element \n";
	if(exists $element2Group->{$element}){
		#$SETN = $element2Group->{$element} ;
		next ;
	}
	else{
		$SETN++;
		$TOTALSETS++;
	    ## add element to SETN 
	    $element2Group->{$element} = $SETN ;
	    $SETS->{$SETN}->{$element} = 1 ;
	}


	## add each element of that element to the GROUP, 
	### and keep on tracking
	my $tab = $info->{$element};
    my @l = (sort keys %{$tab});
	my $N = @l ;
	print "element $element has $N connections\n" if($verbose);
	while(@l){
		my $l = shift @l ;
	    next if(exists $element2Group->{$l});	
		print "\t eeeee $l \n" if($verbose);
	    $element2Group->{$l} = 1 ;
	    $SETS->{$SETN}->{$l} = 1 ;
	 
	    my $tab = $info->{$l};
        my @ll = (sort keys %{$tab});
		push @l, @ll ;
	}
}



my $ofhlist = util_write("grouped.list");
my @sortN ;
my $sortN = {};
foreach my $SETN (1..$TOTALSETS){
    my (@l) = (sort keys %{$SETS->{$SETN}});

	if($verbose){
	foreach my $a (@l){
	   foreach my $b (@l){
	   		next if($a eq $b);
	   	    my $v = $info->{$a}->{$b};
			print "$v .......\n";
	   }
	}
	}
	my $NNNN = @l ;
    print $ofh "SET no: $SETN: $NNNN =N @l \n";
	foreach my $l (@l){
		print $ofhallingroup "$l\n";
	}
    print $ofhfirst "$l[0]\n";
}


print "There were $TOTALSETS groups for cutoff $cutoff \n";
close($ofh);
system ("sort.pl -idx 3 -inf $outfile");



sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
