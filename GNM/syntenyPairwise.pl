#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use MyGeom;
use MyConfigs;
use PDB;
use ConfigPDB;
use Math::Geometry ;
use Math::Geometry::Planar;
use Memory::Usage;
use Algorithm::Combinatorics qw(combinations) ;
my $mu = Memory::Usage->new();
$mu->record('');
use AAConfig;


my $diffAllowed = 2 ;
my $PREFIX1stOrganism = "MDP" ;
my $NMATCHES = 6 ;
my $MIN_LEN_OF_RUN = 10 ;

print "Info: diffAllowed = $diffAllowed, PREFIX1stOrganism = $PREFIX1stOrganism , NMATCHES =$NMATCHES \n";

use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($annofile,$fastafile,$idx,$infile,$postfix,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$protein);
my ($ignorefile,@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "fastafile=s"=>\$fastafile ,
            "postfix=s"=>\$postfix ,
            "infile=s"=>\$infile ,
            "annofile=s"=>\$annofile ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "diffAllowed=i"=>\$diffAllowed ,
            "idx=i"=>\$idx ,
            "verbose=i"=>\$verbose ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
my $ifh = util_read($infile);

usage( "Need to give a listfile -option -listfile  ") if(!defined $listfile);
usage( "Need to give a p1 pdb id -option -p1  ") if(!defined $p1);
usage( "Need to give a p2 pdb id -option -p2  ") if(!defined $p2);
usage( "Need to give a annofile pdb id -option -annofile  ") if(!defined $annofile);
my $CNT = 0 ; 
#my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
 my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP,$BLASTOUT,$BLASTDB) = util_SetEnvVars();


print "Info: diffAllowed between two syntentic matches is $diffAllowed, anytning more than NMATCHES=$NMATCHES considered as moveable elements, PREFIX1stOrganism=$PREFIX1stOrganism \n";



## Get the mapping info for genes for both organisms... Note "$which/LISTS/list.$i is supposed to be ordered!!
my ($infoMapChr1,$infoMapCnt1,$reverseMap1,@listchrNames1) = GetNumberedInfo($p1);
my ($infoMapChr2,$infoMapCnt2,$reverseMap2,@listchrNames2) = GetNumberedInfo($p2);


## Annoation of the first organism genes (currently through Arabidopsis only)
my $annoMapArab = util_mapFullLinetofirst($annofile);
my $annoMapScaff = util_mapFullLinetofirst("MALUS/MALUS2WAL.100.appended.removevalues");

## Now parse the grouping file
my $mappedCounts = {};
my $p1Numbered = {};
my $p2Numbered = {};

my $mapGeneA2GenesinB = {};
### first find out which gene(s) map to which in the two organisms...
while(<$ifh>){
     next if(/^\s*$/);
	 next if(/^\s*#/);

	 my (@l) = split ;
	 my $nmatches = @l -1 ; ## first one is the MDP malus gene

	 next if($nmatches > $NMATCHES);

	 my $cntP1 = 0;
	 my $cntP2 = 0;
	 ## Only one 2 one match - not really doing much with it...
	 if($nmatches eq 2){
	 	my ($g1,$g2) = @l ;
		my $chr1 = $infoMapChr1->{$g1} ;
		my $index1 = $infoMapCnt1->{$g1} ;
		my $chr2 = $infoMapChr2->{$g2} ;
		my $index2 = $infoMapCnt2->{$g2} ;

		## not all genes are mapped to chromosomes
		next if(!defined $chr1 || !defined $chr2);

		$mappedCounts->{$chr1} = {} if(!defined $mappedCounts->{$chr1});
		$mappedCounts->{$chr1}->{$chr2} = 0 if(!defined $mappedCounts->{$chr1}->{$chr2});
		$mappedCounts->{$chr1}->{$chr2} = $mappedCounts->{$chr1}->{$chr2} + 1 ;

	 }


	 my @A ;
	 my @B ;
	 foreach my $gene (@l){
		if($gene =~ /$PREFIX1stOrganism/){
		   my $g1 = $gene ;
		   my $chr1 = $infoMapChr1->{$g1} ;
		   my $index1 = $infoMapCnt1->{$g1} ;
		   ## not all genes are mapped to chromosomes
		   next if(!defined $chr1);
		   $p1Numbered->{$chr1} = {} if(!defined $p1Numbered->{$chr1});
		   $p1Numbered->{$chr1}->{$index1} = 1 ;
		   push @A, $gene ;
		}
		else{
		    my $g2 = $gene ;
		    my $chr2 = $infoMapChr2->{$g2} ;
		    my $index2 = $infoMapCnt2->{$g2} ;
		   ## not all genes are mapped to chromosomes
		   next if(!defined $chr2);
		   $p2Numbered->{$chr2} = {} if(!defined $p2Numbered->{$chr2});
		   $p2Numbered->{$chr2}->{$index2} = 1 ;
		   push @B, $gene ;
		}

	 }
	 foreach my $a (@A){
		  die if(exists $mapGeneA2GenesinB->{$a});
	 	  $mapGeneA2GenesinB->{$a} = {};
	      foreach my $b (@B){
		  	 $mapGeneA2GenesinB->{$a}->{$b} = 1;
	      }
	 }



}
close($ifh);


## Now find synteny for each chromosome
foreach my $k (@listchrNames1){
    my $ofhRunsdebug = util_write("FINAL/runs.$k");
	my $tab = $p1Numbered->{$k} ;
	print $ofh "$k ";
	my @orderedlistofmappedgenes ;
    foreach my $k1 (sort {$a <=> $b} keys %{$tab}){
		print $ofh " $k1 ";
		push @orderedlistofmappedgenes,$k1 ;
	}
    print $ofh " \n";
	ProcessOneChromosome($k,$ofhRunsdebug,@orderedlistofmappedgenes);
}


sub ProcessOneChromosome{
	my ($chromosome,$ofhRunsdebug,@orderedlistofmappedgenes) = @_ ;
	system ("mkdir -p SCAFF_per_CHR");
    my $OFHTMP = util_write("SCAFF_per_CHR/$chromosome");
	my $doneprintofhtmp = {};

	my $prev ; 
	my @runs ;

	my $sorting = {};
	## club ordered numbers based on diffAllowed:
	foreach my $i (@orderedlistofmappedgenes){
		if(!defined $prev){
			$prev = $i ;
			next ;
		}

		my $diff = abs($prev - $i);
		if($diff > $diffAllowed){


			## break run...
			if(@runs> $MIN_LEN_OF_RUN){
			   my @finalruns ;
			   push @finalruns, @runs ;
			   my @nameruns = MapBack($chromosome,$reverseMap1,@finalruns);
			   my $NUM = @nameruns ;
			   my $str = join " ", @nameruns ;
			   $sorting->{$str} = $NUM ;
			}
			## reset ...
			@runs = ();
			push @runs,$i ;
		}
		else{
			push @runs,$i ;
		}
		$prev = $i 
	}


	my $NUM ;
	## now look at the other organism
    foreach my $fullrun (sort {$sorting->{$b} <=> $sorting->{$a}} keys %{$sorting}){
	    my $secondOrganism = {};
	    my $juglansOrganism = {};
	    my $doneOnce = {};
	    my $first ;
	    my $last ;

		## there may be multiple ... just do the max 
		if(!defined $NUM){
	        $NUM = $sorting->{$fullrun};
		}
		else{
			last if($NUM > $sorting->{$fullrun});
		}

		my @genes1 = split " ", $fullrun ;

		my $TOT = @genes1 ;
		foreach my $gene1 (@genes1){
			my $annArabidopsisStr = "$gene1 Not anno to arab";
			if(exists $annoMapArab->{$gene1}){
				$annArabidopsisStr = $annoMapArab->{$gene1} ;
				$annArabidopsisStr =~ s/isNT.*//i;
			}



		    my $chr1 = $infoMapChr1->{$gene1} ;
		    my $index1 = $infoMapCnt1->{$gene1} ;

			## save the first and last to find the range
			$first = $index1 if(!defined $first);
			$last = $index1 ;

			print $ofhRunsdebug "$chr1-$index1 $annArabidopsisStr \n";

			my ($_gene,$ArabidopisAnnoscaff) = split " ",$annArabidopsisStr;

			my @genesin2 = (keys %{$mapGeneA2GenesinB->{$gene1}});


			## search for corresponding genes in the secondOrganism;
			print $ofhRunsdebug "--> @genesin2 ";
			foreach my $g2 (@genesin2){
				my $chr = $infoMapChr2->{$g2};
				my $cnt = $infoMapCnt2->{$g2};
				print $ofhRunsdebug  "$chr $cnt ";
				$secondOrganism->{$chr} = [] if(!defined $secondOrganism->{$chr});

				## just take an unique
				if(!exists $doneOnce->{$ArabidopisAnnoscaff}){
				    push @{$secondOrganism->{$chr}}, $cnt ;
			        print $OFHTMP "$gene1\n" if(!exists $doneprintofhtmp->{$gene1});
					$doneprintofhtmp->{$gene1} = 1;
				}
			}
			print $ofhRunsdebug "\n";


			## not the right place currently?
			if(exists $annoMapScaff->{$gene1}){
				 my $str = $annoMapScaff->{$gene1} ;
				 my @l = split " ", $str ;
				 shift @l ;
				 if(@l){
				     print $ofhRunsdebug "JUG---->  @l\n";

					 my $CNT = 0 ;
					 foreach my $chr (@l){
				         $juglansOrganism->{$chr} = [] if(!defined $juglansOrganism->{$chr});
				         if(!exists $doneOnce->{$ArabidopisAnnoscaff}){
				             push @{$juglansOrganism->{$chr}}, $CNT ;
							 $CNT++;
				         }
					 }
				 }
			}
			
			$doneOnce->{$ArabidopisAnnoscaff} = 1;
		}
		my $UNIQ = (keys %{$doneOnce});

		my $diffinseq = $last - $first ;
		print $ofhRunsdebug "---------- \n";
		print $ofhRunsdebug "$chromosome rang = $diffinseq, total = $TOT, unique = $UNIQ \n";
		print $ofhRunsdebug "---------- \n";

	    foreach my $chrOrganism2 (sort keys %{$secondOrganism}){
	        my @l = @{$secondOrganism->{$chrOrganism2}} ;
		    next if(@l < $MIN_LEN_OF_RUN);
	        my ($N,$first,$last,$range,$mean,$sd) = util_GetStatOfList(@l);
		    print $ofhRunsdebug "secondOrganism $chrOrganism2 @l \t N=$N,range=$range \n";


			## fix, by remobving outliers
			my @s = sort { $a <=> $b} @l;
			my @newl ;
			foreach my $x (@s){
				if(abs($x - $mean) < 2*$sd){
					push @newl, $x ;
				}
			}

			if(@newl){
	            ($N,$first,$last,$range,$mean,$sd) = util_GetStatOfList(@newl);
		        print $ofhRunsdebug "fixedsecondOrganism $chrOrganism2 @newl \t N=$N,range=$range \n\n";

				## MAPPING ...
				if(abs($diffinseq -$range) < 50 && $N*2 > $UNIQ){
					print "a$chromosome b$chrOrganism2\n";
				}
			}
			else{
		        print $ofhRunsdebug "fixedsecondOrganism nothing here \n\n";
			}
	    }

	    foreach my $chr (sort keys %{$juglansOrganism}){
	        my @l = @{$juglansOrganism->{$chr}} ;
		    next if(@l < $MIN_LEN_OF_RUN);
	        my ($N,$first,$last,$range,$mean,$sd) = util_GetStatOfList(@l);
		    print $ofhRunsdebug "juglansOrganism $chr @l \t N=$N,range=$range \n";
	    }
	    print $ofhRunsdebug "\n\n";

	} ## foreach synteny


} ## Foreach chromosome




## from a chromosome name and (a set of) number(s), get the names
sub MapBack{
    my ($chrnm,$reverseMap,@finalruns) = @_ ;
    my @ret ;
    my $tab = $reverseMap->{$chrnm};
    foreach my $i (@finalruns){
  	   my $scaff = $tab->{$i} ;
	   push @ret, $scaff ;
    }
    return @ret 
}

sub GetNumberedInfo{
	my ($which) = @_ ;
    my $L1 = "$which/list";
    my @list= util_read_list_sentences($L1);
	my $infoMapChr = {};
	my $infoMapCnt = {};
	my $reverseMap = {};
    foreach my $i (@list){
	    my $infasta = "$which/LISTS/list.$i";
	    my $ifh = util_read($infasta);
		my $cnt = 1 ;
	    $reverseMap->{$i} = {};
	    while(<$ifh>){
			my ($l) = split ;
			$infoMapChr->{$l} = $i;
			$infoMapCnt->{$l} = $cnt;
	        $reverseMap->{$i}->{$cnt} = $l;
			$cnt++;
	    }
    }
	return ($infoMapChr,$infoMapCnt,$reverseMap,@list) ;
}

sub util_CreateMapFromGroupFile{
  my ($inf,$prefix) = @_ ;
  my $ifh = util_read($inf);
  my $mapSet12Set2 = {};
  while(<$ifh>){
     next if(/^\s*$/);
	 next if(/^\s*#/);

	 #SET no: 56: 2 =N MDP0000026761 POPTR_0014S17270.1 
	 my ($j1,$j2,$j3,$nmatches,$j4,@l) = split ; 

	 #next if($nmatches > $NMATCHES);

	 my @s1 ;
	 my $s2 ;
	 foreach my $gene (@l){
		if($gene =~ /$prefix/){
			push @s1, $gene;
		}
		else{
			$s2->{$gene} = 1 ;
		}
	 }
	 foreach my $g (@s1){
	 	$mapSet12Set2->{$g} = $s2 ;
	 }


   }
   close($ifh);
   return $mapSet12Set2 ;
}

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
