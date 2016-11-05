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
my $mu = Memory::Usage->new();
$mu->record('');


my $debugone = 0 ;
my $savetablebool = 0 ; ## applies only to seqeunce, and not genome

my ($genomefile,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$sequencefile);
my ($isNT,$fastadir,$ignorefile,@expressions);
my $ksize ;
my $verbose =0  ;
my $maskl = 0  ;
my $hopping = 1 ;
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
            "hopping=i"=>\$hopping ,
            "isNT=i"=>\$isNT ,
            "maskl=s"=>\$maskl ,
            "savetablebool=i"=>\$savetablebool ,
            "verbose=i"=>\$verbose ,
            "debugone=i"=>\$debugone ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);

usage( "Need to give a input file name => option -genomefile ") if(!defined $genomefile);
usage( "Need to give a input file name => option -outfile ") if(!defined $outfile);
usage( "Need to give a input file name => option -ksize ") if(!defined $ksize);
usage( "Need to give a input file name => option -isNT ") if(!defined $isNT);
my $info = {};
system("mkdir -p $outfile");
my $outfileanno = "$outfile/$outfile.$maskl.$ksize.anno";
my $outfilelist = "$outfile/$outfile.$maskl.$ksize.list";
my $outfilemap = "$outfile/$outfile.$maskl.$ksize.mapone";
my $outfilemapall = "$outfile/$outfile.$maskl.$ksize.mapall";
#my $outfiledone = "$outfile/$outfile.$maskl.$ksize.done";
my $outfilenotmapped = "$outfile/$outfile.$maskl.$ksize.notmapped";

my $ofhanno = util_write($outfileanno);
my $ofhlist = util_write($outfilelist);
my $ofhmapping = util_write($outfilemap);
my $ofhmappingall = util_write($outfilemapall);
my $ofhnotmapped = util_write($outfilenotmapped);

#print STDERR "Info: Loggin in tmp2monitor\n";

my $STR = "";
my $NAME;



print "Reading DB $genomefile\n";

my $GENOME = new KMERGENOME($genomefile,$ksize,$hopping,$isNT,1);
my $GENOMETABLE = $GENOME->SPLITSEQ();
$GENOME->Print();


my $allsequences = $GENOME->GetAllSequences();

my @list = (sort keys %{$allsequences});

my $Nfiles = @list ;

my $CNT = 0 ;
my $CNTFOUND = 0 ;
print "Searching for fasta $Nfiles \n";
my $mapped = {};
foreach my $i (@list){
    my ($sequence) = $allsequences->{$i} ;
	my $LEN = length($sequence);

    my $kmerSingle = new KMER($i,$sequence,$ksize,$isNT);

	my $retval = $kmerSingle->SanityCheck();
	if($retval eq -1){
		print "Sequence smalller than $ksize\n" if($verbose);
		next ;
	}

	## hack for plants
	if($maskl){
	    $retval = $kmerSingle->SetMASKL($maskl);##
	}

    my ($seqtable) = $kmerSingle->SplitSliding();
    #$kmerSingle->Print();
     

	$CNT++ ;
    my $N1 = (keys %{$seqtable}) ;

	my @matches ;
	my @strings ;
    foreach my $s (keys %{$seqtable}){
	   if(exists $GENOMETABLE->{$s}){
		  my @l = @{$GENOMETABLE->{$s}} ;
		  push @matches,@l;

		  $CNTFOUND++;
		  push @strings,$s ;
	   }
    }
	my $NSTRINGS = @strings ;
	#print "$N1 $NSTRINGS ===========\n";
	next if(!$NSTRINGS);


	my $str = $strings[0];
    my $NElem = util_FindAlphabetStrings($str);

	## this counts the number too
	my ($tableofMatches)  = util_make_table(\@matches);

	my $nNumberofMatchesinGenome = keys %{$tableofMatches};
	my @uniqmatches = (sort {$tableofMatches->{$b} <=> $tableofMatches->{$a}} keys %{$tableofMatches});
	my $v = $uniqmatches[0];
	my $NUM = $tableofMatches->{$v} ;

    print $ofhanno "$i example string $str with NElem $NElem example $v (total= $nNumberofMatchesinGenome ) BLASTP $genomefile TMP.$i\n";
    print $ofhlist "$i\n";
    print $ofhmapping "$i $v $NUM\n";
    print $ofhmappingall "$i @uniqmatches\n";

	## This is for checking mapped when doing self - TODO
	my $TAB = util_make_table(\@uniqmatches);
	my $NNNN = keys %{$TAB};
	if($NNNN > 1){
	        $mapped->{$i} = 1 ;
	}
}



#my $ofhdone = util_write($outfiledone);

foreach my $i (@list){
	if(! exists $mapped->{$i}){
		print $ofhnotmapped "$i\n";
	}
}
system("wc -l $outfile/$outfile.*.$ksize*");

$mu->record('after something_memory_intensive()');
$mu->dump();


exit ;



####################



sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
