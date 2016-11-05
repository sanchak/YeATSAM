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
my ($isNT,$samefile,$infile,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$protein);
my ($ignorefile,@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
my $removeN = 0 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "isNT=s"=>\$isNT ,
            "infile=s"=>\$infile ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            "samefile"=>\$samefile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "removeN=i"=>\$removeN ,
            "verbose=i"=>\$verbose ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
usage( "Need to give a input file name => option -isNT ") if(!defined $isNT);
my $ifh = util_read($infile);
$outfile = "$infile.fixed.fa" if(!defined $outfile);
my $ofh = util_write($outfile);




my $info = {};
my $totallen = 0 ;
my $totallenAfterNremoval = 0 ;
while(<$ifh>){
     if(/^\s*>/){
		my @l = split ;
		my $nm = shift @l ;

		my $orignm = $nm ;
		#$nm =~ s/_//g;
		$nm =~ s/\(//g;
		$nm =~ s/\)//g;
		#$nm =~ s/\.//g;
		$nm =~ s/\|//g;
		$nm =~ s/'//g;
		if($nm ne $orignm){
			print "Fixed $orignm to $nm\n";
        }


		$, = " ";
		print $ofh "$nm @l \n";
	      	
	 }
	 else{
	 	s/ //g;
		chomp ;
	 	$totallen = $totallen + length($_);
		if($isNT && $removeN){
		    $_ =~ s/NN*/N/g;
	 	    $totallenAfterNremoval = $totallenAfterNremoval + length($_);
		}
	 	print $ofh "$_\n" ;
	 }
}


my $diff = $totallen - $totallenAfterNremoval ;
my $percent = (100 * $diff)/$totallen;

print "$totallen = totallen, $totallenAfterNremoval = totallenAfterNremoval, diff = $diff, percent = $percent \n";

close($ifh);
close($ofh);
if(defined $samefile){
    util_printAndDo("cp -f $outfile $infile");
}


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
