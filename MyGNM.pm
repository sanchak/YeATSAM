package MyGNM;
use Carp ;
use POSIX ;
require Exporter;
use Algorithm::Combinatorics qw(combinations) ;
use Math::NumberCruncher;
use Math::MatrixReal;  # Not required for pure vector math as above
use Math::Geometry ; 
use Math::VectorReal qw(:all);  # Include O X Y Z axis constant vectors
#use Math::Trig;
#use Math::Trig ':radial';
no warnings 'redefine';
my $EPSILON = 0.01;
use MyPymol;
use MyConfigs;
use MyUtils;

my $verbose = 0 ;
my $debugone = 0 ;
no warnings 'recursion'; 

  local $SIG{__WARN__} = sub {};

@ISA = qw(Exporter);
@EXPORT = qw( 
GNM_PARSEBLAST
GNM_DirectionBlast
GNM_IsSameScaffolds
GNM_MakeGroupOfNR_OrAnnotate
GNM_MapOneTRS2Scaffold
GNM_parseCountSummed
GNM_parseANNFile
GNM_KMERIZE_sliding
GNM_KMERIZE_hopping
GNM_GenomeBreakHoppingOrSliding
GNM_PARSEBLAST_BESTVAL
        );

use strict ;
use FileHandle ;
use Getopt::Long;




## this functions writes in a group file format - (:$trs $name $blastscore}  (forWGS = 0)
## or write the annotation 

sub GNM_MakeGroupOfNR_OrAnnotate{
    my ($info,$querylength,$Subjectname,$isNT,$findcharstring,$forWGS,$trs,$infile,$ofh,$ofhanno,$blastcutoff,$MAPSscafftoTRS,$checkforsame) = @_ ;
    
    my $sorted = {};
    my $DONEBSCORE = {};


	#print STDERR "Sorting based on BLAST SCORE\n";
    foreach my $k (@{$info}){
        my $org = $k ;
        my ($name,$description,$blastscore,$subjectlength,$iden_percent,$subjectmatched,$querystart,$queryend,$subjectstart,$subjectend,$expect) = split " ",$k ;
        while(exists $DONEBSCORE->{$blastscore}){
            $blastscore = $blastscore - 0.1;
        }
        $DONEBSCORE->{$blastscore} = 1 ;
        $sorted->{$k} = $blastscore ;
    }
     my @sorted = sort { $sorted->{$b} <=> $sorted->{$a} } (keys %{$sorted});
     my $NNN = @sorted ;
    
    my $done = {};
    my @others ;
    my $found100percent = 0;
    foreach my $k (@sorted){
        my $org = $k ;
        my ($name,$description,$blastscore,$subjectlength,$iden_percent,$subjectmatched,$querystart,$queryend,$subjectstart,$subjectend,$expect) = split " ",$k ;
        print "$name,$description,$blastscore,$subjectlength,$iden_percent,$subjectmatched,$querystart,$queryend,$subjectstart,$subjectend,$expect\n" if($verbose);
        $name =~ s/;//g;
        if($name eq $trs && $iden_percent eq 100){
            $found100percent  = 1 ;
            $querylength = $subjectlength;
            $done->{$name} = 1 ;
        }
        else{
            push @others, $k ;
        }

        #while(exists $DONEBSCORE->{$blastscore}){
            #$blastscore = $blastscore - 0.1;
        #}
        #$DONEBSCORE->{$blastscore} = 1 ;
        #$sorted->{$k} = $blastscore ;
    }
    die if (!$found100percent && defined $checkforsame);
    die if (!defined $querylength && defined $checkforsame);

    

     if(!@sorted){
          print $ofhanno "$trs Warning:Did not find any characterized\n";
		  return -1 ;
     }



	### Annotate ###
    my $UNCHARSTRINGS = Config_GetUncharStrings();

    my $first ; 
    my $percent ; 
    foreach my $k (@sorted){
        my $org = $k ;
        my ($name,$description,$blastscore,$subjectlength,$iden_percent,$subjectmatched,$querystart,$queryend,$subjectstart,$subjectend,$expect) = split " ",$k ;
    
        ## choose the first percent
        $description =~ s/ZZZ/ /g;
        my $str2print ;
        if(!defined $first){
    	    if($blastscore < $blastcutoff){
               print $ofhanno "$trs Warning:Did not find any characterized\n";
			   last ;
	        }

			## note when matching to scaffolds this percent might be very low
            $percent = int(100 * ($querylength/$subjectlength));
            $str2print = "$trs\t$name\t$description $percent $blastscore $expect\n";
            $first = $str2print ;
        }
    
        $str2print = "$trs\t$name\t$description isNT=$isNT $percent $blastscore $expect\n";
    
        next if($findcharstring && $description =~ /$UNCHARSTRINGS/i);

		if($blastscore < $blastcutoff){
           print $ofhanno "$first";
        }
		else{
             print $ofhanno "$str2print";
		}
        last ;
    }
        
	########## Make groups ############
    my $grouping = {};
    my $foundone = 0 ;
    my $cntfound = 0 ;
    if($forWGS){
       print $ofh "$trs ";    
     }
	my $ofhCommands = util_open_or_append("pairwiseTRS2scaffold.csh");
	my $last ; 

	my $FOUNDSOMETHING ;
    my $DONENAMES = {};
    foreach my $k (@sorted){
        my $org = $k ;
        my ($name,$description,$blastscore,$subjectlength,$iden_percent,$subjectmatched,$querystart,$queryend,$subjectstart,$subjectend,$expect) = split " ",$k ;
    
        $name =~ s/;//g;
        #$name =~ s/.ORF.*//;
        next if($trs eq $name);
    
        print STDERR "This one has $name $blastscore $subjectlength $subjectmatched $iden_percent.  querylength = $querylength\n" if($verbose);
        if($forWGS){
          
          if(exists $DONENAMES->{$name}){
		  		## this is happening as we are doing collison correction
                #warn "$trs has $blastscore > $DONENAMES->{$name}" if($blastscore > $DONENAMES->{$name});
          }
		  $FOUNDSOMETHING = $blastscore  if(!defined $FOUNDSOMETHING); ## keep the first 
          if($blastscore > $blastcutoff){
              print $ofh " $name $blastscore " ;

			  ## HARDCODED
			  if(! -e "BLASTOUT_SCAFFOLD/$trs.$name.blast"){
                 print $ofhCommands " parseBlastLatestpairwiseStep1.csh $trs $name\n";
			  }
			  $MAPSscafftoTRS->{$name} = {} if(!defined $MAPSscafftoTRS->{$name});
			  $MAPSscafftoTRS->{$name}->{$trs} = 1 ;
              $DONENAMES->{$name} = $blastscore ;
          }
        }
        elsif($blastscore > $blastcutoff){
            $foundone = 1 ;
            if(! exists $done->{$name}){
                print $ofh "$trs $name $blastscore\n" ;
                $done->{$name} = 1;
                $cntfound++;
             }
    
        }
    }
    if($forWGS){
       print $ofh "\n";    
     }
     else{
		 ## No matches
         if(!$foundone){
             #print $ofh "$trs $trs 0 # only one \n";
         }
    }
	return $FOUNDSOMETHING ;
}


sub GNM_PARSEBLAST_BESTVAL{
	my ($infile) = @_ ;
    my ($info,$querylength,$Subjectname,$queryname) = GNM_PARSEBLAST($infile);
    my $maxscore = 0 ;
    my $maxStr = 0 ;
    foreach my $k (@{$info}){
        my $org = $k ;
        my ($name,$description,$blastscore,$subjectlength,$iden_percent,$subjectmatched,$querystart,$queryend,$subjectstart,$subjectend,$expect) = split " ",$k ;
        if($blastscore > $maxscore){
            $maxscore = $blastscore ;
            $maxStr = $k ;
        }
    }
    my ($name,$description,$blastscore,$subjectlength,$iden_percent,$subjectmatched,$querystart,$queryend,$subjectstart,$subjectend,$expect) = split " ",$maxStr ;
    return ($info,$querylength,$Subjectname,$queryname,$blastscore,$expect);
}



## Note that this parses not the first lines (after "sequences") - but the complete file.
sub GNM_PARSEBLAST{
    my ($infile,$verbose) =@_ ;
	$verbose = 0 if(!defined $verbose);
	print "GNM_PARSEBLAST for $infile \n" if($verbose);
    
    my $ifh = util_read($infile);
    my $junk ;
    my $Query ;
    my $Querylength ;
    my $Subjectlength ;
    my @SAVED ;
    my $found = 0 ;
	my $Subjectname ; 

	## Sequences producing significant alignments:  decides pairwise or otherwise
	## This block reads the input file till "Sequences producing significant alignment";
    while(<$ifh>){
        chomp ;
        next if(/^BLAST/);
        next if(/^\s*$/);
        push @SAVED, $_ ;
        if(/Subject=/){
            chomp;
            s/Subject=//;
            $Subjectname = $_ ;
            next ;
        }
        if(/Query= /){
            ($junk,$Query) = split ;
			$Query =~ s/;.*//;
            next ;
        }
        if(/Length=/){
            s/Length=//;
            if(!defined $Querylength){
                $Querylength =$_ ;
            }
            else{
                ## for pw only
                $Subjectlength =$_ ;
            }
            next ;
        }
        if(/Sequences producing significant alignments/){
            $found = 1 ;
            last ;
        }

    }   


    ## All lines are in @SAVED
	## this is pairwise PW
    my @alllines ;
    if(!$found){
        print "Processing pairwise\n" if($verbose);
        @alllines = @SAVED ;
        my $CNT = 0 ;
        my @newalllines ;
        my $subject ;
        my $subjectlen ;
        my $annotation ;
		my @lll ;
        while(@alllines){
            $_ = shift @alllines ;
            if(/^\s*Subject/){
			   s/Subject=//;
               ($subject, @lll) = split ;
			   $annotation = join " ", @lll ;
			   #print "$subject $annotation kkkkkkkkkkkkkkkkkkkkk\n";
			   my $dowhile =  1;
			   while($dowhile){
                   $_ = shift @alllines ;
				   if(/Length=/){
                        s/\s*//g;
                       ($subjectlen) = (/Length=(.*)/);
					   $dowhile = 0 ;
			        }
				}
            }
            if(/Score/){
                push @newalllines, ">$subject $annotation";
                push @newalllines, "Length=$subjectlen";
                $CNT++;
            }
            push @newalllines, $_ ;
        }

        @alllines = @newalllines ;
    }
    else{
    
        while(<$ifh>){
           chomp ;
           push @alllines, $_ ;
        }
    }

    
    my $info = [];
    while(@alllines){
        $_ = shift @alllines ;

		## for pairwise, we have inserted the ">" above
        if(/^\s*>/){
            ParseSingleMatch($ifh,$_,\@alllines,$info);
            $found = 1 ;
            last ;
        }
    }


    return ($info,$Querylength,$Subjectname,$Query);
}


sub ParseSingleMatch{
    my ($IFH,$line,$alllines,$ALLINFO) = @_ ;
    my $N = @{$alllines};
    #print "$alllines $N kkkkkkkkk\n";
    $line =~ s/>//;


	



    my ($name,@l) = split " " ,$line ;
    my $description = join "ZZZ", @l;


    my ($junk,$subjectlength,$blastscore,$expect,$iden_percent,$ratio,$subjectmatched) ;
    my ($querystart,$queryend) ;
    my ($subjectstart,$subjectend) ;

    while(@{$alllines}){
        ## recurse
        $_ = shift @{$alllines};

        if(/^\s*>/){
            ParseSingleMatch($IFH,$_,$alllines,$ALLINFO);
        }
        else{
           if(/Length=/){
               s/Length=//;
               $subjectlength =$_ ;
           }
           elsif(/Score/){
                   if(defined $querystart){
                     $description = "UUU" if($description eq "");
                   my $str = "$name $description $blastscore $subjectlength $iden_percent $subjectmatched  $querystart $queryend $subjectstart $subjectend $expect";
                   #print "XXXXXXX $str\n";
                   push @{$ALLINFO}, $str ;
                }
                   undef $querystart;
                   undef $queryend;
                   undef $subjectstart;
                   undef $subjectend;
                   ($junk,$junk,$blastscore,$junk,$junk,$junk,$junk,$expect)= split ;
                $expect =~ s/,//;
            }
           elsif(/Identities/){
                   ($junk,$junk,$ratio,$iden_percent) = split ;
				   die "$_ " if(!defined $iden_percent);
                   ($junk,$subjectmatched) = split "/", $ratio;
                   $iden_percent =~ s/\%//;
                   $iden_percent =~ s/,//;
                   $iden_percent =~ s/\(//;
                   $iden_percent =~ s/\)//;
           }
           elsif(/Query/){
                 my (@l) = split ;
              my $N= @l -1;
              $querystart = $l[1] if(!defined $querystart);;
              $queryend = $l[$N] ;
           }
           elsif(/Sbjct/){
                 my (@l) = split ;
              my $N= @l -1;
              $subjectstart = $l[1] if(!defined $subjectstart);;
              $subjectend = $l[$N] ;
              #print "$subjectstart $subjectend \n";
           }


        }
    }

    $description = "UUU" if($description eq "");
    my $str = "$name $description $blastscore $subjectlength $iden_percent $subjectmatched  $querystart $queryend $subjectstart $subjectend $expect";
    #print "XXXXXXX $str\n";
    push @{$ALLINFO}, $str ;

}


sub GNM_DirectionBlast{
	my ($infile) = @_ ;
    my $ifh = util_read($infile);
    while(<$ifh>){
		if(/^\s*Strand/){
			if(!/Strand=Plus\/Plus/){
				return 0 ;
			}
			else{
		       return 1 ;
			}
		}
    }
	close($ifh);
	## no hits
	return -1 ;
}

sub GNM_parseCountSummed{
	my ($annofile,$cutoff) = @_ ;
	my $ifh = util_read($annofile);
	my $table = {};
	my $allinfocounts = {};
    while(<$ifh>){
         next if(/^\s*$/);
	     next if(/^\s*#/);
		 my @l = split ;
		 my $N = @l -1 ; 
		 my $nm = $l[0] ;
		 my $mean = $l[3] ;
		 if($mean > $cutoff){
		     $table->{$nm} = $mean ;
		     $allinfocounts->{$nm} =  {};
		     $allinfocounts->{$nm}->{TOTAL} = $l[1];
		     $allinfocounts->{$nm}->{NUMBER} = $l[2];
		     $allinfocounts->{$nm}->{MEAN} = $l[3];
		     $allinfocounts->{$nm}->{SD} = $l[4];
		     $allinfocounts->{$nm}->{RATION} = $l[5];
		 }
		 
	}
	return ($table,$allinfocounts) ;
    
}
sub GNM_parseANNFile{
	my ($annofile,$tableRRNA) = @_ ;
	my $ifh = util_read($annofile);
	my $table = {};
    while(<$ifh>){
         next if(/^\s*$/);
	     next if(/^\s*#/);
		 my @l = split ;
		 my $N = @l -1 ; 
		 my $nm = $l[0];
		 next if(exists $tableRRNA->{$nm});
		 my $eval = $l[$N ] ;
		 my $blastscore = $l[$N - 1] ;
		 my $percent = $l[$N - 2] ;
		 $table->{$nm} = {};
		 $table->{$nm}->{EVAL} = $eval ;
		 $table->{$nm}->{BLASTSCORE} = $blastscore ;
		 $table->{$nm}->{PERCENT} = $percent ;
		 $table->{$nm}->{FULLLINE} = $_ ;
		 if($nm =~ /_A/){
		 	$nm =~ s/_A//;
		    $table->{$nm} = {};
		    $table->{$nm}->{EVAL} = $eval ;
		    $table->{$nm}->{BLASTSCORE} = $blastscore ;
		    $table->{$nm}->{PERCENT} = $percent ;
		    $table->{$nm}->{FULLLINE} = $_ ;
		 }
	}
	return $table ;
    
}

sub GNM_KMERIZE_sliding{
	my ($str,$table,$LEN,$name,$slidingunit,$genome) = @_ ;
	
    #HACK
	if(0){
	    $str =~ s/L//g;
	    $str =~ s/S//g;
	}
	my $len = length($str);
	print "Processing GNM_KMERIZE_sliding with ($len < $LEN) $name \n" if($verbose);

	return if($len < $LEN);
	

	my $start = 0 ;
	my $end = $LEN ;
	my $mapped = 0 ;
	while(1){
	    my $s1 = substr ($str , $start, $LEN)   ;
		die "$name substr ( , $start, $LEN)  \n" if(!defined $s1);
		my $lll = length($s1);

		if(!defined $genome ){ ## save if genome is not defined
			$table->{$s1} = [] if(! exists $table->{$s1});
		    push @{$table->{$s1}}, $name ;
		}
		else{
			if(exists $genome->{$s1}){
				my $genomemap = $genome->{$s1};
				$mapped = 1;
			}
		}

		if($debugone){
		    my $lll = length($s1);
		    print "$lll $start $end \n";
		    print "$s1\n";
		}
		last if($mapped);

		last if($end eq $len );
		$start = $start + $slidingunit ;
		$end = $start + $LEN ;
	}
	die "debugone = $debugone" if($debugone);
	return $mapped ;

}

sub GNM_KMERIZE_hopping{
	my ($maskL,$str,$kmertable,$LEN,$name) = @_ ;

	## HACK
	if($maskL){
	    $str =~ s/$maskL//g;
	}
	#print "Genome string = $str\n";

	my $len = length($str);
	print "Processing GNM_KMERIZE_hopping with ($len < $LEN) $name \n" if($verbose);
	return if($len < $LEN);


	my $start = 0 ;
	my $end = $LEN ;
	while(1){

	    my  $s1 = $end > $len ? substr ($str , $start) : substr ($str , $start, $LEN)   ;

		

		### Add data
		    if(! exists $kmertable->{$s1}){
			    $kmertable->{$s1} = [];
		    }
		    push @{$kmertable->{$s1}}, $name ;



		if($debugone){
		    my $lll = length($s1);
		    print "$lll $start $end \n";
		    print "$s1\n";
		}

		last if($end > $len );

		$start = $start + $LEN ;
		$end = $start + $LEN ;
	}
	die "debugone = $debugone" if($debugone);
}

## genometable might be undefined
sub GNM_GenomeBreakHoppingOrSliding{
	my ($maskL,$INFILE,$ksize,$hopping,$isNT,$genometable) = @_; 
    my $kmertable = {};
    my $NAMES = {};
	die "file $INFILE does ot exist" if(! -e $INFILE);
    my $ifh = util_read($INFILE);
	my $NAME ;
	my $STR ;
    while(<$ifh>){
	     next if(/^\s*$/);
         if(/^\s*>/){
	 	    ## process one fasta
	 	    if(defined $NAME){
	           _ProcessSingleFasta($maskL,$STR,$kmertable,$ksize,$hopping,$isNT,$genometable,$NAME) ;
			   $NAMES->{$NAME} = 1 ;

                $STR = "";
		    }
    
		    s/>//;
		    ($NAME) = split ;
	     }
	     else{
	        die if(/>/);
		    ## remove spaces 
		    s/ //g;
		    chomp;
		    $STR = $STR . $_ ;
	     }
    }
	 _ProcessSingleFasta($maskL,$STR,$kmertable,$ksize,$hopping,$isNT,$genometable,$NAME) ;
	$NAMES->{$NAME} = 1 ;
    close($ifh);
    return ($kmertable,$NAMES);
}

sub _ProcessSingleFasta{
	my ($maskL,$STR,$kmertable,$ksize,$hopping,$isNT,$genometable,$NAME) = @_ ;
	my $len = length($STR);
	if($hopping){
       GNM_KMERIZE_hopping($maskL,$STR,$kmertable,$ksize,$NAME);
	}
	else{
		my $slidingunit = 1 ;
        my $mapped = GNM_KMERIZE_sliding($STR,$kmertable,$ksize,$NAME,$slidingunit,$genometable);
		if(!$mapped && $isNT){
		die ;
            my $rev = util_getComplimentaryString($STR) ;
            GNM_KMERIZE_sliding($rev,$kmertable,$ksize,$NAME,$slidingunit,$genometable);
		}
	}
}

sub GNM_IsSameScaffolds{
	my ($scafftable,$a,$b)  = @_ ;
	$a =~ s/_(A|B|C)//;
	$b =~ s/_(A|B|C)//;
	$a =~ s/\.MER\..*//;
	$b =~ s/\.MER\..*//;
	die "$a not in scafftable" if(! exists $scafftable->{$a});
	die "$b not in scafftable" if(! exists $scafftable->{$b});

	my $sa = $scafftable->{$a};
	my $sb = $scafftable->{$b};
	return 1 if($sa eq $sb);
	return 0 ;
}
