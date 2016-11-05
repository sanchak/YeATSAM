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
my ($forR,$reverse,$infile,$outfile,$which_tech,$seperator,$listfile);
my (@idx);
$seperator = " ";
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "infile=s"=>\$infile ,
            "reverse"=>\$reverse ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "forR=s"=>\$forR ,
            "seperator=s"=>\$seperator ,
            "idx=i"=>\@idx,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);

usage( "Need to give a input file name => option -infile ") if(!defined $infile);
my $ifh = util_read($infile);
if(!@idx){
	push @idx, 0;
}

my $idx1 = $idx[0];

$outfile = "$infile.$idx1" if (! defined $outfile);
my $ofh = util_write($outfile);
my $CNT = 0 ; 
my ($RESULTDIR,$PDBDIR,$FASTADIR) = util_SetEnvVars();


my $first = 1 ;
my $info = {};


my $IDX = shift @idx ;
while(<$ifh>){
	 s/\)//g;
     next if(/^\s*$/);
	 my (@l) = split " " ,$_ ;

	 my $v = $l[$IDX];

	 if(! exists $info->{$v}){
	 	$info->{$v} = 1 ;
	 }
	 else{
	 	$info->{$v} = $info->{$v} +  1;
	 }
}

foreach my $k (sort { $info->{$a} <=> $info->{$b} } keys %{$info}){
	my $v = $info->{$k};
	print $ofh "$k $v \n";
}



print STDERR "Output written in $outfile\n";

chmod 0777, $outfile ;
sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
