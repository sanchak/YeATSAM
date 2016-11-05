#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use Carp;
use MyUtils;
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($infile,$outfile,$which_tech,$listfile);
my (@expressions);
my $howmany = 100000 ;
my $infile1;
my $infile2; 
my ($idx1,$idx2);
$idx1 = $idx2 = 0 ; 
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "in1=s"=>\$infile1 ,
            "in2=s"=>\$infile2 ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "verbose=i"=>\$verbose ,
            "idx1=i"=>\$idx1 ,
            "idx2=i"=>\$idx2 ,
           );
die "Give 2 input files  " if(@ARGV != 2);
$infile1 = shift @ARGV ;
$infile2 = shift @ARGV ;


my $info1 = ReadInfile($infile1,$idx1);
my $info2 = ReadInfile($infile2,$idx2);

my $N1  = (keys %{$info1}) ;
my $N2  = (keys %{$info2}) ;
print "N1 = $N1 , N2 = $N2\n" if($verbose);

sub ReadInfile{
    my ($infile,$idx) = @_ ;
	croak " file $infile does not exist" if(! -e $infile);
    my $ifh = util_read($infile);
	my $info = {};
    while(<$ifh>){
         next if(/^\s*$/);
         chomp ;
	     my (@l) = split ; 
	     my $num = @l - 1 ;
	     my $str = $l[$idx];
	     #$str =~ s/.pdb.out//;
		 $info->{$str} = $_ ;
    }
	return $info ; 
}
my $ofhboth = util_write("ofhboth");
my $ofhinAbutnotinB = util_write("ofhinAbutnotinB");
my $ofhinBbutnotinA = util_write("ofhinBbutnotinA");

my $CNT = 0 ; 
foreach my $i (sort keys %{$info1}){
	if(!exists $info2->{$i}){
		print $ofhinAbutnotinB "$i\n";
	    $CNT++;
	}
}
#print  "$CNT exists in 1 but not 2 : in file ofhinAbutnotinB \n" if($verbose);

$CNT = 0 ;
foreach my $i (sort keys %{$info2}){
	if(!exists $info1->{$i}){
		print  $ofhinBbutnotinA "$i\n";
	    $CNT++;
	}
}
#print  "$CNT exists in 2 but not 1 in file ofhinBbutnotinA\n" if($verbose);

$CNT = 0 ;
foreach my $i (sort keys %{$info2}){
	if(exists $info1->{$i}){
	    my $l1 = $info1->{$i} ; 
	    my $l2 = $info2->{$i} ; 
		#print  $ofhboth "$i\n";
		print  $ofhboth "$l1\n";
	    $CNT++;
	}
}
#print  "$CNT exists in both :in files ofhboth.1 and ofhboth.2 \n" if($verbose);
system ("newfile.csh ofhConcat");
system ("cat ofhboth ofhinBbutnotinA ofhinAbutnotinB >>  ofhConcat");
system ("wc -l ofh*");

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
