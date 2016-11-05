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
my ($ksize,$isNT,$fastadir,$ignorefile,@expressions);
my $verbose =0  ;
my $exitaftergenome =0  ;
my $maskl = 0  ;
my $nelemLimit = 3 ;
my $hopping ;
my $reverse ;
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
            "ksize=f"=>\$ksize ,
            "hopping=i"=>\$hopping ,
            "isNT=i"=>\$isNT ,
            "exitaftergenome=i"=>\$exitaftergenome ,
            "maskl=s"=>\$maskl ,
            "savetablebool=i"=>\$savetablebool ,
            "verbose=i"=>\$verbose ,
            "nelemLimit=i"=>\$nelemLimit ,
            "debugone=i"=>\$debugone ,
            "reverse=i"=>\$reverse ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);

usage( "Need to give a input file name => option -infile ") if(!defined $infile);
usage( "Need to give a input file name => option -outfile ") if(!defined $outfile);
usage( "Need to give a input file name => option -ksize ") if(!defined $ksize);

my $ksizetab = {};
$ksizetab->{14} = 132269611 ;
$ksizetab->{8} = 102269611 ;
$ksizetab->{7} = 82269611 ;

die if(!exists $ksizetab->{$ksize});
my $filesize = $ksizetab->{$ksize};
my $ifh = util_read($infile);
my $map2name = {};
$map2name = util_mapID2fullStringInFastaFile($infile);
my ($NAME,$sequence);
my $seenNAMES = {};
my $currentCount = 0 ;
my $allCount = 0 ;
my $idx = 0 ;
my $ofhlist = util_write("$outfile.$ksize.list");
my $OFH = util_write("$outfile.$ksize.$idx");
while(<$ifh>){
	     next if(/^\s*$/);
         if(/^\s*>/){
	 	    if(defined $NAME){
			   if(! exists $seenNAMES->{$NAME} ){
	               ProcessSingleSequence($NAME,$sequence);
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
	    ProcessSingleSequence($NAME,$sequence);
	}
	close($ifh);


## Wrapper to have it run on the complimentary	
sub ProcessSingleSequence{
	my ($i,$sequence) = @_ ;
	my $len = length($sequence);
	$allCount = $allCount + $len ;
	$currentCount = $currentCount + $len ;
	
	## open a new file 
	if($currentCount > $filesize){
		$idx++;
		close($OFH);
		$OFH = util_write("$outfile.$ksize.$idx");
		print $ofhlist "$outfile.$ksize.$idx\n";
		$currentCount = 0 ;

	}
	my $nm = $map2name->{$NAME} ;
	print $OFH "$nm\n";
	print $OFH "$sequence\n";
}




sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
