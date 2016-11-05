#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use MyUtils;
use MyConfigs;
use Memory::Usage;
use MyGNM;
use KMER;
use KMERGENOME;
use Algorithm::Combinatorics qw( combinations variations_with_repetition) ;
use Algorithm::Combinatorics qw( permutations ) ;
my $mu = Memory::Usage->new();
$mu->record('');


my $debugone = 0 ;
my $savetablebool = 0 ; ## applies only to seqeunce, and not genome

my ($infile,$genomefile,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$sequencefile);
my ($isNT,$fastadir,$ignorefile,@expressions);
my $ksize ;
my $verbose =0  ;
my $exitaftergenome =0  ;
my $maskl = 0  ;
my $nelemLimit = 3 ;
my $hopping ;
my $reverse ;
my $maxnumber = 3 ;
GetOptions(
            "sequencefile=s"=>\$sequencefile ,
            "fastadir=s"=>\$fastadir ,
            "genomefile=s"=>\$genomefile ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            "infile=s"=>\$infile ,
            "which_tech=s"=>\$which_tech ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "ksize=i"=>\$ksize ,
            "hopping=i"=>\$hopping ,
            "isNT=i"=>\$isNT ,
            "exitaftergenome=i"=>\$exitaftergenome ,
            "maskl=s"=>\$maskl ,
            "savetablebool=i"=>\$savetablebool ,
            "verbose=i"=>\$verbose ,
            "nelemLimit=i"=>\$nelemLimit ,
            "debugone=i"=>\$debugone ,
            "maxnumber=i"=>\$maxnumber ,
            "reverse=i"=>\$reverse ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);

usage( "Need to give a input file name => option -infile ") if(!defined $infile);
$outfile = "$infile.orf" if(!defined $outfile);
my $ofh = util_write($outfile);

my $ifh = util_read($infile);
	my ($NAME,$sequence);
	my $seenNAMES = {};
    while(<$ifh>){
	     next if(/^\s*$/);
         if(/^\s*>/){
	 	    if(defined $NAME){
			   if(! exists $seenNAMES->{$NAME} ){
	               $_ = util_WriteORFSingleFasta($NAME,$sequence,$ofh,$_,$maxnumber);
			       $seenNAMES->{$NAME} = 1 ;
			   }
			   else{
			        warn "Warning: repeat name $NAME\n" ;
			   }
               $sequence = "";
		    }
    
		    s/>//;
		    ($NAME) = split ;
	     }
	     else{
	        die if(/>/);
		    ## remove spaces 
		    s/ //g;
		    chomp;
		    $sequence = $sequence . $_ ;
	     }
    }
	if(! exists $seenNAMES->{$NAME} ){
	    $_ = util_WriteORFSingleFasta($NAME,$sequence,$ofh,$_,$maxnumber);
	}
	close($ifh);

