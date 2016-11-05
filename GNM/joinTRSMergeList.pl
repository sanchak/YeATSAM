#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use MyGeom;
use MyConfigs;
use PDB;
use MyGNM;
use ConfigPDB;
use Math::Geometry ;
use Math::Geometry::Planar;
use Memory::Usage;
use Algorithm::Combinatorics qw(combinations) ;
my $mu = Memory::Usage->new();
$mu->record('');
use AAConfig;



use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($fastadir,$idx,$DBfile,$annfile,$p1,$p2,$outfile,$cutoff,$which_tech,$mergedir,$blastdir);
my ($ignorefile,@expressions);
my $ignoresamescaff = 1 ;
my $verbose = 0 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "blastdir=s"=>\$blastdir ,
            "fastadir=s"=>\$fastadir ,
            "annfile=s"=>\$annfile ,
            "DBfile=s"=>\$DBfile ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "mergedir=s"=>\$mergedir ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "ignoresamescaff=i"=>\$ignoresamescaff ,
            "idx=i"=>\$idx ,
            "verbose=i"=>\$verbose ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);
usage( "Need to give a input file name => option -DBfile ") if(!defined $DBfile);
usage( "Need to give a input file name => option -fastadir ") if(!defined $fastadir);
usage( "Need to give a input file name => option -annfile ") if(!defined $annfile);

my $MINEVALE =  1E-10 ;
my $PERCENTCUTOFF_EVALUE =  0.9 ;

my $joinlog = "$mergedir/joinlog.$ignoresamescaff";
my $joingood = "$mergedir/good.$ignoresamescaff";
my $joinbad = "$mergedir/bad.$ignoresamescaff";
my $newannfile = "$mergedir/newannfile.$ignoresamescaff";

my $ofhlog = util_write($joinlog);
my $ofhgood = util_write($joingood);
my $ofhbad = util_write($joinbad);
my $ofhann = util_write($newannfile);
my $ofhnewlist = util_write("$mergedir/newlist.$ignoresamescaff");

usage( "Need to give a mergedir -option -mergedir  ") if(!defined $mergedir);
usage( "Need to give a blastdir pdb id -option -blastdir  ") if(!defined $blastdir);

my $listfile = "$mergedir/MERGEALL";

my $anninfo = {};
my $anninfofullline = {};
my $annfh = util_read($annfile);
## Read annotations - for evlue, and later writing the new annotations
while(<$annfh>){
     next if(/^\s*$/);
	 chomp ;
	 next if(/^\s*#/);

	 my (@l) = split ; 
	 my $lastidx = @l-1;
	 my $nm = $l[0];
	 my $last = $l[$lastidx];
	 $anninfo->{$nm} = $last ;
	 $anninfofullline->{$nm} = $_ ;
}


##Read mergeable TRSs
## Ignore TRS with possible repetitions
my $ifh = util_read($listfile);
my $runblasts  = 0 ;
while(<$ifh>){
     next if(/^\s*$/);
	 next if(/^\s*#/);

	 my ($str,$f1,$f2) = split ; 
	 my $newnm = "$f2.MER.$f1";
	 my $blastresult = "$blastdir/$newnm.blast.nt";
	 my $len = length($str);
	 my $nElem = util_FindAlphabetStrings($str);
	 my $allowed = $len eq 5 ? 2 : 3 ;
	 if($nElem <$allowed){
	 	print $ofhlog "$f1 $f2 $str is repetitive so next\n";
		next ;
	 }
	 if(! -e $blastresult){
	 	$runblasts = 1 ;
		last ;
	 }
}

print "runblasts = $runblasts\n";

my ($scafftable) = util_maketablefromfile("map.full.TRS2Scaffold");

$ifh = util_read($listfile);
my $finallymerged = {};
my $newMERGED = {};
while(<$ifh>){
     next if(/^\s*$/);
	 next if(/^\s*#/);

	 my ($str,$f1,$f2) = split ; 
	 my $newnm = "$f2.MER.$f1";
	 my $blastresult = "$blastdir/$newnm.blast.nt";
	 my $len = length($str);
	 my $nElem = util_FindAlphabetStrings($str);
	 my $allowed = $len eq 5 ? 2 : 3 ;
	 if($nElem <$allowed){
		next ;
	 }

	 my ($str1) = util_readfasta("$fastadir/$f1.ALL.1.fasta");
	 my ($str2) = util_readfasta("$fastadir/$f2.ALL.1.fasta");
	 die if(!($str1 =~ /^$str/));
	 die if(!($str2 =~ /$str$/));
     
     $str1 =~ s/^$str//;
	 my $finalstr = $str2 . $str1 ;
	 my $L1 = length($str1);
	 my $L2 = length($str2);


	 my $e1 = exists $anninfo->{$f1} ?  $anninfo->{$f1} : -1 ;
	 my $e2 = exists $anninfo->{$f2} ?  $anninfo->{$f2} : -1 ;

     ## write fasta each time - otherwise gets complicated
	      #my ($str1) = util_readfasta("$fastadir/$f1.ALL.1.fasta");
	      #my  ($str2) = util_readfasta("$fastadir/$f2.ALL.1.fasta");
	      #die if(!($str1 =~ /^$str/));
	      #die if(!($str2 =~ /$str$/));
     
          #$str1 =~ s/^$str//;
	      #my $finalstr = $str2 . $str1 ;
     
	      my $newfasta = "$fastadir/$newnm.ALL.1.fasta";
	      print "Writing $newnm\n" if($verbose);
	      my $ofhnewfast = util_write($newfasta);
	      print $ofhnewfast ">$newnm MERGE=$str \n";
	      print $ofhnewfast "$finalstr\n";
     
          
	      print " $e1 $e2 \n" if($verbose);
	      close($ofhnewfast);

	 my $newevalue = 1000  ;
	 my $DESC ;
	 if($runblasts){
	    if(! -e $blastresult){
	      print $ofh "BLASTP $DBfile $newfasta $blastresult\n";
	    }
	 }
	 else{
	 	 ## BLASTs are done, now process the data
	     my ($info,$querylength,$Subjectname,$queryname) = GNM_PARSEBLAST($blastresult);
		 foreach my $k (@{$info}){
		     my $org = $k ;
			 my ($name,$description,$blastscore,$subjectlength,$iden_percent,$subjectmatched,$querystart,$queryend,$subjectstart,$subjectend,$expect) = split " ",$k ;

			if($newevalue > $expect){
			    $newevalue = $expect  ;
				$DESC = $description ;
			    $DESC =~ s/ZZZ/ /g;
			}
			
		}
		$DESC = "XXX" if(!defined $DESC);

		my $plsplrint = 0 ;
		my $IsSameScaffolds = $ignoresamescaff ? GNM_IsSameScaffolds($scafftable,$f1,$f2) : 1 ;

		## Ignore any TRS that has been merged

		if(exists $finallymerged->{$f1} || exists $finallymerged->{$f2}){
			$plsplrint = 0 ;
		}
		else{
			 #print "$newevalue $blastresult\n";
			 #print " $IsSameScaffolds && $newevalue < $MINEVALE && $newevalue < $e1 && $newevalue  < $e2 \n";
			 
		    if($IsSameScaffolds && $newevalue < $MINEVALE && $newevalue < $e1 && $newevalue  < $e2){

			  my $E1 = $e1 ;;
			  my $E2 = $e2 ;;
			  $E1 =~ s/.E-//i;
			  $E2 =~ s/.E-//i;
			  my $TOTAL = $E1 + $E2 ;
		  	 if($newevalue eq "0.0"){
			 	$plsplrint = 1 ;
			 }
			 else{
			 	my $NEWVALUE = $newevalue ;
				$NEWVALUE =~ s/.E-//i;
				if($NEWVALUE > $PERCENTCUTOFF_EVALUE*$TOTAL){
			         $plsplrint = 1 ;
				}
				#print "jjj $TOTAL $NEWVALUE \n";

			 }

		}
		}
			if($plsplrint){
		        print $ofhgood "$newnm $L1 $L2 $str e1=$e1 e2=$e2 newevalue=$newevalue $DESC\n";
				$finallymerged->{$f1} = $str ;
				$finallymerged->{$f2} = $str ;
				$newMERGED->{$newnm} = "$DESC $newevalue";
			}
			else{
		        print $ofhbad "$newnm $L1 $L2 $str e1=$e1 e2=$e2 newevalue=$newevalue $DESC\n";
			}

	 }
	if($e1 eq -1 && $e2 eq -1 ){
	 	print $ofhlog "$str $f1 and $f2 are not anno ,so next\n";
	 }
}


## Remove both TRSs being merged
foreach my $k (sort keys %{$anninfo}){
	if(! exists $finallymerged->{$k}){
		my  $v = $anninfofullline->{$k};
		print $ofhann  "$v\n";
		print $ofhnewlist  "$k\n";
	}
}

## write the new merged TRS id
foreach my $k (sort keys %{$newMERGED}){
		my  $v = $newMERGED->{$k};
		print $ofhann  "$k $v\n";
		print $ofhnewlist  "$k\n";
}

system("wc -l $outfile");
system("wc -l $joinlog");
system("wc -l $joingood");

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
