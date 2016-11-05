#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use ConfigPDB;
use MyGeom;

use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($absolute,$reverse,$infile,$cutoff,$outfile,$which_tech,$listfile);
my (@binnumbers);
my $idx = 0;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "infile=s"=>\$infile ,
            "absolute"=>\$absolute ,
            "reverse"=>\$reverse ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "binnumbers=s"=>\@binnumbers,
            "idx=i"=>\$idx ,
            "cutoff=i"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
usage( "Need to give a input file name => option -binnumbers ") if(! @binnumbers);
my $ifh = util_read($infile);
my $CNT = 0 ; 
my ($RESULTDIR,$PDBDIR,$FASTADIR) = util_SetEnvVars();


my $numbers  = {};
my $warned = 0 ;
while(<$ifh>){
     next if(/^\s*$/);
	 my (@l) = split ;
	 my $N = @l -1 ;
	 next if($idx > $N);

	 my $NNN = $l[$idx];
	 next if(defined $cutoff && $NNN < $cutoff );
	 if($NNN eq "nan"){
		$warned ++;
		next ;
	 }
	 $NNN = abs($NNN) if(defined $absolute);
	  $numbers->{$NNN} = [] if(!defined  $numbers->{$NNN});
	 push @{$numbers->{$NNN}}, $_ ;
}
close($ifh);

if($warned){
	print "ignored $warned nan's\n";
}


my $bincnt = 0 ;
my $ofh = util_write("bin$bincnt");
my $prevvalue = 0 ;


my @sortedbin = sort {$a <=> $b} @binnumbers ;
my $NUMBEROFBINS = @sortedbin;
my $maxvalue = $sortedbin[$bincnt] ;
foreach my $i (sort {$a <=> $b} keys %{$numbers}){
	    if($i >= $prevvalue && $i < $maxvalue){
		}
		else{
			$prevvalue = $sortedbin[$bincnt];
			if($bincnt < $NUMBEROFBINS -1){
			     $bincnt++;
                 $maxvalue = $sortedbin[$bincnt] ;
			}
			else{
				$bincnt++;
				$maxvalue = 1000000 ;
			}
            $ofh = util_write("bin$bincnt");

		}

        print $ofh @{$numbers->{$i}};
}
print STDERR "Output written in $infile.bin..$bincnt\n";

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
