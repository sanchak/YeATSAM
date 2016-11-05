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

my ($queryfile,$splitfile,$genomefile,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$sequencefile);
my ($isNT,$fastadir,$ignorefile,@expressions);
my $ksize ;
my $verbose =0  ;
my $maskl = 0  ;
my $exitaftergenome =0  ;
my $nelemLimit = 3 ;
my $hopping ;
my $reverse = 0 ;
my $postprocess = 0 ;
GetOptions(
            "sequencefile=s"=>\$sequencefile ,
            "fastadir=s"=>\$fastadir ,
            "genomefile=s"=>\$genomefile ,
            "splitfile=s"=>\$splitfile ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            "queryfile=s"=>\$queryfile ,
            "which_tech=s"=>\$which_tech ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "ksize=i"=>\$ksize ,
            "hopping=i"=>\$hopping ,
            "isNT=i"=>\$isNT ,
            "maskl=s"=>\$maskl ,
            "savetablebool=i"=>\$savetablebool ,
            "verbose=i"=>\$verbose ,
            "nelemLimit=i"=>\$nelemLimit ,
            "debugone=i"=>\$debugone ,
            "exitaftergenome=i"=>\$exitaftergenome ,
            "reverse=i"=>\$reverse ,
            "postprocess=i"=>\$postprocess ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);

usage( "Need to give a input file name => option -genomefile ") if(!defined $genomefile);
usage( "Need to give a input file name => option -outfile ") if(!defined $outfile);
usage( "Need to give a input file name => option -ksize ") if(!defined $ksize);
usage( "Need to give a input file name => option -isNT ") if(!defined $isNT);
usage( "Need to give a input file name => option -hopping ") if(!defined $hopping);


if($postprocess){
	PrintPostProcess();
	exit ;
}

system("mkdir -p $outfile");
my $outfileanno = "$outfile/$outfile.$maskl.$ksize.anno";
my $outfiledone = "$outfile/$outfile.$maskl.$ksize.done";
my $outfilelist = "$outfile/$outfile.$maskl.$ksize.list";
my $outfilemap = "$outfile/$outfile.$maskl.$ksize.mapone";
my $outfilemapall = "$outfile/$outfile.$maskl.$ksize.mapall";
my $outfilenotmapped = "$outfile/$outfile.$maskl.$ksize.notmapped";
my $outfileOrigCommand = "$outfile/$outfile.$maskl.$ksize.origcommand.csh";

if(-e $outfileanno){
	PrintPostProcess();
	die "File $outfileanno exists, are you sure ?";
}

my $ofhanno = util_write($outfileanno);
my $ofhlist = util_write($outfilelist);
my $ofhmapping = util_write($outfilemap);
my $ofhmappingall = util_write($outfilemapall);
my $ofhnotmapped = util_write($outfilenotmapped);
my $ofhdone = util_write($outfiledone);
my $ofhorigcommand = util_write($outfileOrigCommand);
print $ofhorigcommand "kmer.pl -quer $queryfile -genome $genomefile -outf $outfile -ks $ksize -hop $hopping -isNT=$isNT \n";


my $ARGS = "-isNT=$isNT -outf $outfile -ks $ksize -hop $hopping -reverse 0 -verb $verbose -exitaftergenome $exitaftergenome";
if(!defined $splitfile){
    system("testkmer_lower.pl -quer $queryfile -geno $genomefile $ARGS");
}
else{
	 my $splitfilelist = "$splitfile";
	 if(! -e $splitfilelist){
	 	die ;
	 }
	 else{
         ##system("splitFastaInMultipleFiles.pl -in $genomefile -out $splitfile -ksize $ksize");
	     my @list= util_read_list_sentences($splitfile);
		 foreach my $i (@list){
               die "File $i does not exist" if(!-e $i);
			   system("nfasta $i");
		 }
		 sleep(1);
		 foreach my $i (@list){
               system("testkmer_lower.pl -quer $queryfile -geno $i $ARGS");
		 }
	 }
}

my ($doneTable) = util_maketablefromfile($outfiledone);
my ($mapped) = util_maketablefromfile($outfileanno);

foreach my $i (sort keys %{$doneTable}){
	if(! exists $mapped->{$i}){
		print $ofhnotmapped "$i\n";
	}
}

system("printInColumnUniquely.pl -in $outfile/$outfile.$maskl.$ksize.mapall");
system("printInPairwiseFromMapfile.pl -inf $outfile/$outfile.$maskl.$ksize.mapall");
system("groupBasedonCutoff.pl -in $outfile/$outfile.$maskl.$ksize.mapall.pw -out $outfile/$outfile.$maskl.$ksize.mapall.pw.group -cutoff 0");

system("wc -l $outfile/$outfile.*.$ksize*");
system("extractindexfromfile.pl -in $outfile/$outfile.$maskl.$ksize.anno");



exit ;

sub PrintPostProcess{
    system("processPWmatchforOne2Many.pl -in $outfile/$outfile.0.$ksize.mapall -outf $outfile/$outfile.0.$ksize.runPW.csh -f1 FASTADIR_NT/ -f2 FASTADIR_NT/ -outd TMP -get 1 -isNT $isNT");
    #print "processPWmatchforOne2Many.pl -in $outfile/$outfile.0.$ksize.mapall -outf $outfile/$outfile.0.$ksize.final.anno -f1 FASTADIR_NT/ -f2 FASTADIR_NT/ -outd TMP -get 0 -genomefile $genomefile -isNT $isNT\n";
	print "processPWmatchALLinOneGo.pl -in $outfile/$outfile.0.$ksize.mapall -out qqq -f1 FASTADIR_NT/ -f2 FASTADIR_NT/ -outd TMP -isNT 0 -geno DB/plantpep.fasta\n";
}


####################



sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
