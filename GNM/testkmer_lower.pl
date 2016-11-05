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

my ($queryfile,$genomefile,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$sequencefile);
my ($isNT,$fastadir,$ignorefile,@expressions);
my $ksize ;
my $verbose =0  ;
my $exitaftergenome =0  ;
my $maskl = 0  ;
my $nelemLimit = 3 ;
my $hopping ;
my $reverse = 0 ;
GetOptions(
            "sequencefile=s"=>\$sequencefile ,
            "fastadir=s"=>\$fastadir ,
            "genomefile=s"=>\$genomefile ,
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

usage( "Need to give a input file name => option -genomefile ") if(!defined $genomefile);
usage( "Need to give a input file name => option -outfile ") if(!defined $outfile);
usage( "Need to give a input file name => option -ksize ") if(!defined $ksize);
usage( "Need to give a input file name => option -isNT ") if(!defined $isNT);
usage( "Need to give a input file name => option -hopping ") if(!defined $hopping);

my ($tmp) = Config_getCodonTable();
my @tmp = (values %{$tmp});
my $AA = util_make_table(\@tmp);
my @AA = (keys %{$AA});
my $NN = @AA ;
my $ignorestrings = {};
foreach my $aa (@AA){
	my $str = "";
	foreach my $idx (1..$ksize){
		$str = $str . $aa ;
	}
	$ignorestrings->{$str} = 1 ;
}

my @IGNORE ;
my $UNKNOWN ;
if($isNT){
   my @NTs = qw (A T G C);
   push @IGNORE,  FindRepititive(1,12,@NTs);
   push @IGNORE, FindRepititive(2,5,@NTs);
   push @IGNORE, FindRepititive(3,3,@NTs);
   push @IGNORE, FindRepititive(4,3,@NTs);
   $UNKNOWN = "N";
}
else{
   if(1){
   push @IGNORE,  FindRepititive(1,8,@AA);
   push @IGNORE, FindRepititive(2,3,@AA);
   #push @IGNORE, FindRepititive(3,3,@AA);
   #push @IGNORE, FindRepititive(4,2,@AA);
   }
   $UNKNOWN = "X";
}

my $igntable = util_maketablefromlist(\@IGNORE);
my $expr2ignore = "(";
foreach my $k (sort keys %{$igntable}){
	$expr2ignore = $expr2ignore . "$k|";
}
$expr2ignore = $expr2ignore . "$UNKNOWN)";

sub FindRepititive{
	my ($N,$nREP,@NTs) = @_ ;
    my $expr2ignore = "";
    my $iter = combinations(\@NTs, $N);
	my @RET;
    while (my $c = $iter->next) {
            my @combo = @{$c} ; 
            #my $iter1 = permutations(\@combo,$N);
			my $iter1 = variations_with_repetition(\@combo, $N);
            while (my $d = $iter1->next) {
                my @combo1 = @{$d} ; 
		        my $x = join "",@combo1;
			    my $multx = util_ConcatStringNTimes($x,$nREP);
				push @RET, $multx ;
		    }
    }
	return @RET ;
}


print "Info:expr2ignore = $expr2ignore \n";
print "Info:nelemLimit = $nelemLimit \n";
print "Info:reverse = $reverse \n";




die "Cant have both listfile and queryfile defined" if (defined $listfile && $queryfile);

my ($tablefasta,$infoSeq2PDB,$mapChainedName2Name) ;
if(defined $listfile){
   usage( "Need to give a input file name => option -fastadir ") if(!defined $fastadir);
}
else{
   usage( "Need to give a input file name => option -queryfile ") if(!defined $queryfile);
   #print "Reading $queryfile\n";
   #($tablefasta,$infoSeq2PDB,$mapChainedName2Name) = util_parsePDBSEQRESNEW($queryfile,0,0);
}

my $info = {};
system("mkdir -p $outfile");
my $outfileanno = "$outfile/$outfile.$maskl.$ksize.anno";
my $outfiledone = "$outfile/$outfile.$maskl.$ksize.done";
my $outfilelist = "$outfile/$outfile.$maskl.$ksize.list";
my $outfilemap = "$outfile/$outfile.$maskl.$ksize.mapone";
my $outfilemapall = "$outfile/$outfile.$maskl.$ksize.mapall";
my $outfilenotmapped = "$outfile/$outfile.$maskl.$ksize.notmapped";

my $ofhanno = util_open_or_append($outfileanno);
my $ofhlist = util_open_or_append($outfilelist);
my $ofhmapping = util_open_or_append($outfilemap);
my $ofhmappingall = util_open_or_append($outfilemapall);
my $ofhnotmapped = util_open_or_append($outfilenotmapped);
my $ofhdone = util_open_or_append($outfiledone);

#print STDERR "Info: Loggin in tmp2monitor\n";

my $STR = "";
my $NAME;


my $name2seq ;

if(defined $listfile){
    my @list= util_read_list_words($listfile);
    my $allfilesthere = 1 ;
    foreach my $i (@list){
	    my $sequencefile = "$fastadir/$i.ALL.1.fasta";
	    if (! -e $sequencefile){
	        print "$sequencefile does not exist\n" ; 
		    $allfilesthere = 0; 
	    }
	    else{
             my ($sequence,$firstline) = util_readfasta($sequencefile);
		     $name2seq->{$i} = $sequence ;
	    }
    }
    die "Files missing" if(!$allfilesthere);
    my $NFiles = @list ;
    print "There are $NFiles to match, all good\n";
}




print "Reading DB $genomefile\n";

my $delfactor = 2 ;
my $GENOME = new KMERGENOME($genomefile,$ksize,$hopping,$isNT,0,$reverse,$delfactor);
my $GENOMETABLE = $GENOME->SPLITSEQ();
$GENOME->Print($verbose);

$mu->record('after something_memory_intensive()');
$mu->dump();
if($exitaftergenome){
   exit ;
}


my $mapped = {};
my $ofhreadAlreadyAnnotated = util_read($outfileanno);
while(<$ofhreadAlreadyAnnotated>){
	my ($x) = split ;
	$mapped->{$x} = 1 ;
}



print "Searching for fasta\n";
if(defined $listfile){
    foreach my $i (sort keys %{$name2seq}){
	    #my $sequencefile = "$fastadir/$i.ALL.1.fasta";
	    my $sequence = $name2seq->{$i} ;
	    ProcessSingleSequence($i,$sequence);
    }
}
else{
	my $ifh = util_read($queryfile);
	my ($NAME,$sequence);
	my $seenNAMES = {};
    my $docheck = 1 ;
    while(<$ifh>){
	     next if(/^\s*$/);
         if(/^\s*>/){
	 	    if(defined $NAME){
			   if(! exists $seenNAMES->{$NAME} ){
	               print $ofhdone "$NAME\n";
	               $docheck = ProcessSingleSequence($NAME,$sequence,$docheck);
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
	    print $ofhdone "$NAME\n";
	    $docheck = ProcessSingleSequence($NAME,$sequence,$docheck);
	}
	close($ifh);
}


## Wrapper to have it run on the complimentary	
sub ProcessSingleSequence{
	my ($i,$sequence,$docheck) = @_ ;
    my $len = length($sequence);
	if($docheck && $len > 60 ){
		$docheck = 0 ;
		util_CheckSequence($i,$sequence,$isNT);
	}

	return $docheck if(exists $mapped->{$i});
	ProcessSingleSequenceLowerLevel($i,$sequence);

	return $docheck if(exists $mapped->{$i});
	print "Runnin COMPLIMENTARY\n" if($verbose);
	my $complimentarysequence = util_getComplimentaryString($sequence) ;
	ProcessSingleSequenceLowerLevel($i,$complimentarysequence);
	return $docheck ;
}

sub ProcessSingleSequenceLowerLevel{
	my ($i,$sequence) = @_ ;

    my $CNT = 0 ;
    my $CNTFOUND = 0 ;
    my $kmerSingle = new KMER($i,$sequence,$ksize,$isNT);

	my $retval = $kmerSingle->SanityCheck();
	if($retval eq -1){
		print "Sequence smalller than $ksize\n" if($verbose);
		return ;
	}

	## hack for plants
	if($maskl){
	    $retval = $kmerSingle->SetMASKL($maskl);##
	}

    my ($seqtable) = $kmerSingle->SplitSliding();
    #$kmerSingle->Print();

	$CNT++ ;
    my $N1 = (keys %{$seqtable}) ;

	my @idMatches ;
	my @strings ;
    foreach my $s (keys %{$seqtable}){
	   if($verbose){
	         my $complimentarysequence = util_getComplimentaryString($s) ;
	         print "File $i, string $s, compl = $complimentarysequence \n" ;
	   }
	   next if(exists $ignorestrings->{$s});

	   ## Only for nt - imp to get nt value right, else you will ignore all asn
	   if($s =~ /$expr2ignore/){
		      next ;
	   }
	   if(exists $GENOMETABLE->{$s}){
		     my @l = @{$GENOMETABLE->{$s}} ;
			 if(@l eq 1){
			     if($i eq $l[0]){
				 	next ;
				 }
			 }
		     $CNTFOUND++;

		     push @idMatches,@l;
		     push @strings,$s ;
	   }
    }
	return if(!@strings);
	my $str = $strings[0];
    my $NElem = util_FindAlphabetStrings($str);
	return if($isNT && $NElem < $nelemLimit);
	
	$mapped->{$i} = 1 ;




	## this counts the number too
	my ($tableofMatches)  = util_make_table(\@idMatches);
	my $nNumberofMatchesinGenome = keys %{$tableofMatches};

	my @uniqmatches = (sort {$tableofMatches->{$b} <=> $tableofMatches->{$a}} keys %{$tableofMatches});
	my $v = $uniqmatches[0];
	my $NUM = $tableofMatches->{$v} ;

    print $ofhanno "$i example string $str with NElem $NElem example $v (total= $nNumberofMatchesinGenome ) BLASTP $genomefile  TMP.$i\n";
    print $ofhanno "@strings" if($verbose);
    print $ofhlist "$i\n";
    print $ofhmapping "$i $v $NUM\n";
    print $ofhmappingall "$i @uniqmatches\n";

	## This is for checking mapped when doing self - TODO
}





foreach my $i (sort keys %{$name2seq}){
	if(! exists $mapped->{$i}){
		print $ofhnotmapped "$i\n";
	}
}
system("printInColumnUniquely.pl -in $outfile/$outfile.$maskl.$ksize.mapall");
system("printInPairwiseFromMapfile.pl -inf $outfile/$outfile.$maskl.$ksize.mapall");
system("groupBasedonCutoff.pl -in $outfile/$outfile.$maskl.$ksize.mapall.pw -out $outfile/$outfile.$maskl.$ksize.mapall.pw.group -cutoff 0");

system("wc -l $outfile/$outfile.*.$ksize*");



exit ;



####################



sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
