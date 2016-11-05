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
my $info = {};

my $ofh = util_write("$infile.appended");

my @nms ;
while(<$ifh>){
	 my (@l) = split " " ,$_ ;
	 my $nm = shift @l ;
	 push @nms, $nm ;
	 if(!exists $info->{$nm}){
	 	 $info->{$nm} = [];
	 }
	 push @{$info->{$nm}}, @l ;
}

my $done= {};
foreach my $nm (@nms){
	next if(exists $done->{$nm});
	$done->{$nm} = 1;
	my @l = @{$info->{$nm}} ;
	print $ofh "$nm @l \n";
}


print "Wroting to $infile.appended\n";

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
