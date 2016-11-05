#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use MyUtils;
use MyConfigs;
use Memory::Usage;
use MyGNM;
use KMER;
my $mu = Memory::Usage->new();
$mu->record('');


my $debugone = 0 ;
my $savetablebool = 0 ; ## applies only to seqeunce, and not genome

my ($genomefile,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$sequencefile);
my ($isNT,$fastadir,$ignorefile,@expressions);
my $ksize ;
my $verbose =0  ;
GetOptions(
            "sequencefile=s"=>\$sequencefile ,
            "fastadir=s"=>\$fastadir ,
            "genomefile=s"=>\$genomefile ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            "which_tech=s"=>\$which_tech ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "ksize=i"=>\$ksize ,
            "isNT=i"=>\$isNT ,
            "savetablebool=i"=>\$savetablebool ,
            "verbose=i"=>\$verbose ,
            "debugone=i"=>\$debugone ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);

usage( "Need to give a input file name => option -genomefile ") if(!defined $genomefile);
usage( "Need to give a input file name => option -outfile ") if(!defined $outfile);
usage( "Need to give a listfile -option -listfile  ") if(!defined $listfile);
usage( "Need to give a input file name => option -ksize ") if(!defined $ksize);
usage( "Need to give a input file name => option -isNT ") if(!defined $isNT);
usage( "Need to give a input file name => option -fastadir ") if(!defined $fastadir);


my $ofhcomm = util_write("$outfile/commands.csh");
my @K = qw (0 A C D E F G H I K L M N P Q R S T V W Y);
foreach my $i (@K){
	 #my $BACK = $i? " & " : " " ;
	 my $BACK  = " ";
     print $ofhcomm " testkmer.pl -genom $genomefile -out $outfile -lis $listfile -ksize $ksize -isNT $isNT -fasta $fastadir -mask $i $BACK \n";
}


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
