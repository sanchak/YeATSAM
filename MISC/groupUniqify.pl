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
my (@expressions);
my $howmany = 100000 ;
my $verbose = 0 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "verbose=i"=>\$verbose ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);

system ("touch logclose") if(! -e "logclose");
my $ofhlog = util_append("logclose");

usage( "Need to give a input file name => option -infile ") if(!defined $infile);
my $ifh = util_read($infile);
my $CNT = 0 ; 
#my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
 my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP,$BLASTOUT,$BLASTDB) = util_SetEnvVars();



my $info = {};
my $fullline = {};
while(<$ifh>){
     next if(/^\s*$/);
     chomp ;
     my ($nm,$junk,$val) = split ; 
     $info->{$nm} = {}  if(! defined $info->{$nm});
     $info->{$nm}->{$junk} = $val ;
}

my @nms = (sort keys %{$info});


my $done = {};
while(@nms){
        my $a = shift @nms ;
	my @A2B = (sort keys %{$info->{$a}});
        foreach my $b (@A2B){
		my $SA = $a.$b;
		my $SB = $b.$a;
		#next if(exists $done->{$SA} || exists $done->{$SB});
		$done->{$SA} = 1 ;
		$done->{$SB} = 1 ;

		next if(!(exists $info->{$a}->{$b} ||  exists $info->{$b}->{$a}));

		my $val ;
                if( exists $info->{$a}->{$b} &&  exists $info->{$b}->{$a}){
			$val =  $info->{$a}->{$b} > $info->{$b}->{$a} ? $info->{$a}->{$b} : $info->{$b}->{$a}  ;
			print "Choosing $val among  $info->{$a}->{$b} > $info->{$b}->{$a} \n" if($verbose);
			$CNT++;
		}
		else{
			if(exists $info->{$a}->{$b}){
			      $val = $info->{$a}->{$b};
			}
			elsif (exists $info->{$b}->{$a}){
			      $val = $info->{$b}->{$a};
			}
		}
		if(!defined $val){
			die " $a $b \n";
		}
		print $ofh "$a $b $val\n";
        }
}

print "ignored $CNT repeat values \n";

chmod 0777, $outfile ;

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
