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
use MyConfigs;


use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($duplicate,$blastout,$infile,$p1,$p2,$outfile,$trs,$cutoff,$listfile,$protein);
my ($chooseone,$ignorefile,$donone,@expressions);
my $howmany = 100000 ;
my $verbose = 0 ;
GetOptions(
            "trs=s"=>\$trs ,
            "blastout=s"=>\$blastout ,
            "infile=s"=>\$infile ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            "outfile=s"=>\$outfile ,
            "donone"=>\$donone ,
            "duplicate"=>\$duplicate ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "verbose=i"=>\$verbose ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
#usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
#my $ofh = util_write($outfile);

my $orfdir = "FASTADIR_ORF";
usage( "Need to give a trs -option -trs  ") if(!defined $trs);
usage( "Need to give a blastout -option -blastout  ") if(!defined $blastout);
usage( "Need to give a cutoff -option -cutoff  ") if(!defined $cutoff);


## When checking for different annotations within a same file, have a cutoff of this 
my $EVALFORSECOND = 1E-15 ;
my $EVALFORCOMMON = 1E-10 ;

my $BLASTEXT = ".blast.nt";

system("mkdir -p ANN/ERR");
system("mkdir -p ANN/WARN");
system("mkdir -p ANN/ANNOTATE");
system("mkdir -p ANN/COMMANDS");
system("mkdir -p ANN/LISTS");
system("mkdir -p FASTADIR_ORFBEST");


my $ofhwarn = util_write("ANN/WARN/$trs.w");
my $ofhannotate = util_write("ANN/ANNOTATE/$trs.a");
my $ofhcommands = util_write("ANN/COMMANDS/$trs.c");
my $ofhlist = util_write("ANN/LISTS/$trs.list");


my $trsfasta = "$trs.ALL.1.fasta";
my $trsnt = "$trs.blast.nt";

my @list= util_read_list_sentences("$orfdir/$trs.list");

my @BLASTFILES ;
my @orffiles ;
foreach my $i (@list){
	my $inblast = "$blastout/$i" ."$BLASTEXT";
	die "$inblast not there " if(! -e $inblast);
	my $infasta = "$orfdir/$i.ALL.1.fasta";
	die if(! -e $infasta);
	push @BLASTFILES, $inblast ;
	push @orffiles, $infasta ;
}


print "Files are @BLASTFILES \n" if($verbose);


my $sort = {};
my $sortEVALS2SCORES = {};
my $sortEVALS2ALLSTRS = {};
my $sortEVALS2INFILE = {};
my $MAXBLASTBITSCORE = 0 ;
my $ALLSTRS = "";

my $sortTRSLen = {};

foreach my $idx (0..2){
   my $infile = $BLASTFILES[$idx];
   my $OOO = $orffiles[$idx];
   my $orfile = $list[$idx];
   my ($fulllength,$TRSnm,$STRS,$LLLs,$SCORES,$IDENTITIES,$EVALUES) = util_ParseWebBlast($infile);
   if(!defined $fulllength ) {
      my $ofhERR = util_write("ANN/ERR/$trs.e");
   	  print $ofhERR "File not complete, rerun $trs\n";
   	  exit ;
   }


	my ($str,$firstline) = util_readfasta($OOO);
	my $len = length($str);
	$sortTRSLen->{$orfile} = $len;
	#print "Mapping $orfile to $len thru sortTRSLen \n" if($verbose);

   ## No hits here ... it has to be here after the above few lines
   next if(!defined $TRSnm);

   my @TRSnm = @{$TRSnm};
   my @STRS = @{$STRS};
   my @LLLs = @{$LLLs};
   my @SCORES = @{$SCORES};
   my @EVALUES = @{$EVALUES};
   my @IDENTITIES = @{$IDENTITIES};

   $infile =~ s/^$blastout.//;
   $ALLSTRS = $ALLSTRS . "\n" . "\t$infile ". $STRS[0];



   my $EEE = $EVALUES[0] ;
   print "EVALUES = $EEE\n" if($verbose);
   while(exists $sort->{$EEE}){
	  print "incrementing $EEE $trs \n"  if($verbose);
   	  $EEE  = $EEE + $EEE/10000 + 1E-100 ;
   }

   $sort->{$EEE} = $STRS ;
   $sortEVALS2SCORES->{$EEE} = $SCORES[0] ;
   $sortEVALS2ALLSTRS->{$EEE} = $STRS ;
   $sortEVALS2INFILE->{$EEE} = $infile ;
   if($MAXBLASTBITSCORE < $SCORES[0]){
   	  $MAXBLASTBITSCORE = $SCORES[0];
   }
   
}

print "Max score = $MAXBLASTBITSCORE   \n" if($verbose);

my $UNCHARSTRINGS = Config_GetUncharStrings();

my $cnt = 0 ;
my $trsrealname  ;
my $firsteval = {};

foreach my $k (sort { $a <=> $b} keys %{$sort}){
	print "Score are $k \n" if($verbose);
}
foreach my $k (sort { $a <=> $b} keys %{$sort}){

	my $SCORE = $sortEVALS2SCORES->{$k};
	if(!defined $SCORE){
        my $ofhERR = util_write("ANN/ERR/$trs.e");
	    print $ofhERR "$k $trs has not SCORE " ;
		die ;
	}
	my $infile = $sortEVALS2INFILE->{$k};
	$infile =~ s/$BLASTEXT//;


	my $len = $sortTRSLen->{$infile};
	print "$infile blastscore=$SCORE len=$len eval=$k\n" if($verbose);


	print "Doing $k < $cutoff \n" if($verbose);
	if(($k < $cutoff)){
		## only for the first 
		if(!defined $trsrealname){
	        $trsrealname = $infile ;
			$firsteval = $k ;
		}


		$cnt++;
	    my @sortedLIST = @{$sort->{$k}};

		## try to get an charaterization - only if the evalue does not fall off
		foreach my $i (@sortedLIST){
			my $EEE = getEvalue($i);
			next if($EEE > $cutoff);

			if(!($i =~ /($UNCHARSTRINGS)/i)){
				last ;
             }
		}

	}
}

print "cnt is $cnt \n" if($verbose);
if($cnt eq 0){
		## choose the longest ....
		my @sorted = sort {$sortTRSLen->{$b} <=> $sortTRSLen->{$a}} (keys %{$sortTRSLen});
		my $NNN = @sorted ;

		my $trsrealname = $sorted[0];
		die "$trsrealname is not found in $trs with $NNN" if(! $trsrealname);

		print "trsrealname = $trsrealname \n" if($verbose);


		my $len = $sortTRSLen->{$trsrealname};
	    #print $ofhcommands "\\cp -f  $orfdir/$trsrealname.ALL.1.fasta FASTADIR_ORFBEST/$trsfasta # none \n";
	    print $ofhcommands "fixFastaNames.pl -trs $trs -inf $orfdir/$trsrealname.ALL.1.fasta -out FASTADIR_ORFBEST/$trsfasta \n";
	    print $ofhcommands "\\cp -f BLASTOUT_AA//$trsrealname.blast.nt  BLASTOUT_AA/$trsnt \n";
	    print $ofhlist "$trs\n";
	    print $ofhannotate "$trs #none none \n";
	    print "None - choosing $trsrealname with len = $len\n" if($verbose);
        if($verbose){
            foreach my $k (keys %{$sortTRSLen}){
                my $v = $sortTRSLen->{$k} ;
                print "$k $v\n";
            }
        }

        exit ;
}

my @sorted = (sort {$a <=> $b} keys %{$sortEVALS2INFILE});
if(!@sorted){
    my $ofhERR = util_write("ANN/ERR/$trs.e");
	print $ofhERR "Expect a value here \n";
	die ;
}

## Evalue of best match
my $X = shift @sorted;
if($cnt eq 1){		

	my $SSS = "uniqhigh";
	if($X > 1E-10){
		$SSS = "uniqlow";
	}
	my @SCORESA = @{$sortEVALS2ALLSTRS->{$X}};
	my $STRA  = $SCORESA[0];

	print $ofhannotate "$trs #$SSS $STRA \n";
	#print $ofhcommands "\\cp -f  $orfdir/$trsrealname.ALL.1.fasta FASTADIR_ORFBEST/$trsfasta # uniq \n";
    print $ofhcommands "fixFastaNames.pl -trs $trs -inf $orfdir/$trsrealname.ALL.1.fasta -out FASTADIR_ORFBEST/$trsfasta \n";
	print $ofhcommands "\\cp -f BLASTOUT_AA//$trsrealname.blast.nt  BLASTOUT_AA/$trsnt \n";
	print $ofhlist "$trs\n";
	print "$SSS $trsrealname\n" if($verbose);

}
elsif($cnt > 1 ){

		my $Y = shift @sorted;

		my $Z ;

		my @SCORESA = @{$sortEVALS2ALLSTRS->{$X}};
		my @SCORESB = @{$sortEVALS2ALLSTRS->{$Y}};
		my @SCORESC ;

		my $STRA  = $SCORESA[0];
		my $STRB  = $SCORESB[0];
		my $STRC  ;


		my ($Atable,@Atable) = getTableofUniprotIds(\@SCORESA,$EVALFORCOMMON);
		my ($Btable,@Btable) = getTableofUniprotIds(\@SCORESB,$EVALFORCOMMON);

		my $AEqualB = 0 ;
		my $AEqualC = 0 ;
		my $BEqualC = 0 ;
		foreach my $uniprotID (@Atable){
			if(exists $Btable->{$uniprotID}){
				$AEqualB = $uniprotID;
				last ;
			}
		}

		my $Ctable = {};
		my @Ctable = {};
		if($cnt eq 3){
		    $Z = shift @sorted;
		    @SCORESC = @{$sortEVALS2ALLSTRS->{$Z}};
		    $STRC = $SCORESC[0];
		    ($Ctable,@Ctable) = getTableofUniprotIds(\@SCORESC,$EVALFORCOMMON);

		    foreach my $uniprotID (@Atable){
			   if(exists $Ctable->{$uniprotID}){
				   $AEqualC = $uniprotID;
				   last ;
			   }
			}
		    foreach my $uniprotID (@Btable){
			   if(exists $Ctable->{$uniprotID}){
				   $BEqualC = $uniprotID;
				   last ;
			   }
			}
		}

		my $foundcommon = 0 ;

		### All (either 2 or 3) point to the same ID 
		if($AEqualB && ($AEqualC || $cnt eq 2)){
				my $OFHMER = util_open_or_append("ANN/mergeUniqFixIssue.csh");
				my $uniprotID = $AEqualB ;
				my $EXTRA = $AEqualC ? $BLASTFILES[2] : " ";
				my $str2print = $AEqualC ? "allthree" : "alltwo";
				print $OFHMER "mergeUniqFixIssue.pl $uniprotID FASTADIR_ORF FASTADIR_PLANTPEP $BLASTFILES[0] $BLASTFILES[1] $EXTRA  #$str2print \n";

	            print $ofhannotate "$trs #uniqfix $STRA  \n";
	            print $ofhcommands "fixFastaNames.pl -trs $trs -inf $orfdir/$trsrealname.ALL.1.fasta -out FASTADIR_ORFBEST/$trsfasta \n";
	            print $ofhcommands "\\cp -f BLASTOUT_AA//$trsrealname.blast.nt  BLASTOUT_AA/$trsnt \n";
	            print $ofhlist "$trs\n";

				exit ;
				
		}


	    my $trsorig = $trs.  ".ALL.1.fasta";
		my @printAnnotate ;
		if(1){
		    my $a = $sortEVALS2INFILE->{$X};
		    my $b = $sortEVALS2INFILE->{$Y};
	        my $trsaonlyname = $trs. "_A" ;
	        my $trsbonlyname = $trs. "_B" ;
	        my $trsa = $trsaonlyname.  ".ALL.1.fasta";
	        my $trsb = $trsbonlyname.  ".ALL.1.fasta";

            $a =~ s/$BLASTEXT/\.ALL.1.fasta/;
            $b =~ s/$BLASTEXT/\.ALL.1.fasta/;

		    my $BLASTA = $sortEVALS2INFILE->{$X};
		    my $BLASTB = $sortEVALS2INFILE->{$Y};
	        print $ofhcommands "fixFastaNames.pl -trs $trsaonlyname -inf $orfdir/$a -out FASTADIR_ORFBEST/$trsa \n";
	        print $ofhcommands "\\cp -f BLASTOUT_AA//$BLASTA  BLASTOUT_AA/$trsaonlyname.blast.nt \n";
		    print $ofhcommands "\\cp -f  FASTADIR_NT/$trsorig FASTADIR_NT/$trsa \n";
	        print $ofhlist "$trsaonlyname  \n";
			push @printAnnotate, "$trsaonlyname #dupli $STRA  \n";

			if(!$AEqualB){

	           print $ofhcommands "fixFastaNames.pl -trs $trsbonlyname -inf $orfdir/$b -out FASTADIR_ORFBEST/$trsb \n";
	           print $ofhcommands "\\cp -f BLASTOUT_AA//$BLASTB  BLASTOUT_AA/$trsbonlyname.blast.nt \n";
		       print $ofhcommands "\\cp -f  FASTADIR_NT/$trsorig FASTADIR_NT/$trsb \n";
	           print $ofhlist "$trsbonlyname  \n";
	           #print $ofhannotate "$trsbonlyname #dupli $STRB  \n";
			   push @printAnnotate, "$trsbonlyname #dupli $STRB  \n";
			}
			else{
				### Do not merge if there is duplicate
				if(0){
				   my $OFHMER = util_open_or_append("ANN/mergeUniqFixIssue.csh");
				   my $uniprotID = $AEqualB ;
				   print $OFHMER "mergeUniqFixIssue.pl $uniprotID FASTADIR_ORF FASTADIR_PLANTPEP $BLASTFILES[0] $BLASTFILES[1] #AEqualB\n";
				}
			}

            if($cnt eq 3)
		    {
		        my $c = $sortEVALS2INFILE->{$Z};
                $c =~ s/$BLASTEXT/\.ALL.1.fasta/;
	            my $trsconlyname = $trs . "_C" ;
	            my $trsc = $trsconlyname . ".ALL.1.fasta";
		        my $BLASTC = $sortEVALS2INFILE->{$Z};
     
				if(!$AEqualC && !$BEqualC){
		            print $ofhcommands "\\cp -f  FASTADIR_NT/$trsorig FASTADIR_NT/$trsc \n";
	                print $ofhcommands "fixFastaNames.pl -trs $trsconlyname -inf $orfdir/$c -out FASTADIR_ORFBEST/$trsc \n";
     
	                print $ofhcommands "\\cp -f BLASTOUT_AA//$BLASTC  BLASTOUT_AA/$trsconlyname.blast.nt \n";
	                print $ofhlist "$trsconlyname  \n";
			        push @printAnnotate, "$trsconlyname #dupli $STRC  \n";
				}
				else{
					die if($AEqualC && $BEqualC);

				    ### Do not merge if there is duplicate
					if(0){
			           $printAnnotate[0] = "$trsaonlyname #dupli_fixB $STRA  \n";
				       my $OFHMER = util_open_or_append("ANN/mergeUniqFixIssue.csh");
				       my $uniprotID = $AEqualC? $AEqualC : $AEqualB ;
				       my $file = $AEqualC? $BLASTFILES[0] : $BLASTFILES[1];
				       my $str2print = $AEqualC? "AEqualC" : "AEqualB";
				       print $OFHMER "mergeUniqFixIssue.pl $uniprotID FASTADIR_ORF FASTADIR_PLANTPEP $file $BLASTFILES[2] #$str2print \n";
					}
				}

		    }
			foreach my $s (@printAnnotate){
				print $ofhannotate $s ;
			}
		    print $ofhcommands "\n";
		}
	}

sub getUniprotID{
	my ($str) = @_ ;
	#$str =~ s/\|/ /g;
	my @l = split " ", $str ;
	return $l[0];
}
sub getEvalue{
	my ($str) = @_ ;
	$str =~ s/\|/ /g;
	my @l = split " ", $str ;
	my $N = @l - 1;
	my $EEE = $l[$N];
	return $EEE ;
}
   
sub getTableofUniprotIds{
		    ## make table ofA loAwer uniprot IDs 
		    my ($l,$CUTOFF) = @_ ;
		    my @SCORESB = @{$l};
		    my $Atable = {};
		    my @Atable ;
			#print " HHHH \n";
		    foreach my $str (@SCORESB){
			    my $uniprotID = getUniprotID($str);
			    my $EEE = getEvalue($str);
                next if($EEE > $CUTOFF);

			    $Atable->{$uniprotID} = $EEE ;
				push @Atable,$uniprotID ;
				#print "$uniprotID $EEE\n";
		    }
			return ($Atable,@Atable) ;
}

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
