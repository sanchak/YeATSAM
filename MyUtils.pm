package MyUtils;
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

my $verbose = 0 ;

  local $SIG{__WARN__} = sub {};

@ISA = qw(Exporter);
@EXPORT = qw( util_execute_chop util_exit_if_doesnt_exist util_get_tag util_print_n_log util_get_user 
            util_get_cmdline util_num_lines util_percentage util_percentages util_printAndDo util_get_pwd util_log_base util_check_machine
			util_debug_printNList
        util_get_tech util_printAndWrite2Script 
        util_get_time_from_string util_isDesignSuccessful util_parseTimelog util_diff util_round
        parseFPocketFile 
        util_print_mgc_home util_get_mgc_home util_write_list_in_file util_fullline 
		util_GetOrf util_WriteORFSingleFasta
		util_ConcatStringNTimes
		util_maketablefromlist
util_read_list_firstitem
util_mapID2fullStringInFastaFile
        util_table_merge  util_table_diff     
		util_FindAlphabetStrings
		util_ParseContactFile
		util_mapFullLinetofirst
		util_ConvertNucleotide2Amino
		util_subtract2tables
        util_MapOneTRS2Scaffold
		util_ParseBlastPW
        util_read_list_numbers
        util_GetSequentialSetResidues
        util_ceil util_floor
		util_maptofirst
        util_getComplimentaryString
        util_get_max 
        util_GetTerminalStrings
        util_GetDistancesBetween2SetsOfResidues
		util_HELIXParseAHBS 
        util_GetDistancesBetween2SetsOfAtoms
        util_PrintNewAtom
        ParseSingleMatch
		util_SortResiduesBasedonNumber
util_CreateTableAndListFromSymbolSeperated

util_get_maxinlist
util_splitstring 
util_ProcessOneORF
util_ProcessOneORFGivenExactCDNA
util_GetCodonBiasFromNucleotide
util_SplitSliding util_SplitHopping

util_R_extractIndex
        util_extractSliceFromFastaString
        util_extractSliceFromFasta
		util_GetStatOfList
        util_GetCentreOfMassFromSet
        util_open_or_append
        util_Blastout_isrev
        util_parsePDBSEQRESNEW
        util_GetCentreOfMass
        util_ConcatTwoStringsWithCommonBeginAndEnd
        util_ReadPremonOut
        util_GetDisulphide
        util_ParsePremonIn
		util_FixFastaNames
        util_IsPro util_IsMSE 
        util_ReadCodonBiasFile
        util_EmitCodonBiasfromInfo
        util_helixWheel
        util_ReadPremonOutAndGiveFasta
        util_IsResisdueThisType
        util_verilog2eqn util_runSis util_two2thepowerof util_emitGrid
        util_wait_on_lockfiles util_writePsp util_get_grid_cmd util_print_synp_file
        util_writelist2file
        util_readAPBSPotentialFromStart
        util_sortsingleString
        util_readfasta
        util_ReadAnnotateFile util_ReadAnnotateFileFRAGAL
        utils_parseBlastOut util_ParseWebBlast util_PARSEBLAST
        util_parseHolo2Apo
        vprint vprintheader vSetVerbose vIncrVerbose vDecrVerbose
        util_copy_table
        util_GetClosestAtoms_intwoPDBs
        util_ParseDSSP
        util_Ann2Simple
        util_AlignAndMatchRemaining
        util_GetMeanSD
        parseCastp
        util_AlignAndMatchRemainingAtoms
util_AreResiduesContinuous
util_FindRmsd
util_FindRmsdAllAtoms
util_ProcessSingleLine
util_parseHETPDB
util_ProcessSingleLineAtomlist

util_RunHelixWheel
util_ParseAAGroups
util_GetPotForAtom
util_GetPotDiffForAtoms
util_CheckSequence
util_GetPotDiffForResidues
util_ReadPdbs
util_WriteClustalAln
        util_WriteFastaFromAtoms
        util_WriteFastaFromResidueNumbers
        
        util_print_vars util_uniq2
        util_SAVEDIR util_INITNAMES
        util_read util_append util_write util_pick_random_from_list
        util_make_list util_make_table
        util_GetFasta
        util_wget util_get_pdb 
        util_makeCSH
        util_parse_pdbseqres
        util_readPDB
        util_read_list_words
        util_read_list_sentences
        util_enter_maxcutoff
        util_is_integer util_is_float util_EnterName util_EnterNumber util_EnterTwoNumbers util_EnterTwoNames 
        util_AddBeforeEach
        util_ignoreIfDistanceIsLessthan util_ignoreIfAnyAreEqual util_ignoreIfAnyAreFurtherThanCutoff
        util_split_list_numbers
        util_read_Mapping_PDB_2_SWISSPROT util_filter_basedon_EC util_getPDBID_basedon_SP util_getECfromPDB
        util_ReadLine
        util_ParseBlocks util_ParseBlockForString
        util_ProcessRowAndColumnsForMean
        util_SetEnvVars
        util_ScoreDistance util_ScorePD
        util_format_float
        util_mysubstrDontStripSpace util_mysubstr
        util_printResult
        util_printTablePre util_printTablePost util_printTableLine util_PrintMeanAndSD
        util_Annotate util_ReadAnnfile
        util_table_print
        util_pick_n_random_from_list
        util_ExtractSliceFromFasta
        util_GetFastaFiles 
        util_maketablefromfile  util_maketablefromfile_firstentry

util_NeedleSeq
util_NeedleFiles
util_NeedlePDBNamesFromSeq
util_NeedlePDBNamesFromFASTADIR
        util_readAPBSPotential
        util_usage util_CmdLine
        util_parsePDBSEQRES
        util_SortTwoStrings
        util_GetEnv
        util_IsZero
        util_Banner util_PrintInfo

        util_getTmpFile
        util_PrintOutConf

        util_sortuniqArray

ParseAPBSResult
        util_printHtmlHeader util_printHtmlEnd util_HtmlizeLine
        util_HtmlTableHead util_HtmlTableEnd util_HtmlTableCell
        util_MakeLink
        util_EC_CreateLevels util_EC_AddPDB util_EC_CorrelatePDBS
        util_PairwiseDiff
util_readPeptideInfo
util_getECID

util_ExtractOneIdxFromFile
        );

use strict ;
use FileHandle ;
use Getopt::Long;


my $havetokeepthispostive = 13 ;

sub util_SetEnvVars{
   my @vars = qw ( RESULTDIR PDBDIR FASTADIR APBSDIR FPOCKET SRC MATCH3D ANNDIR UNIPROT PREMONITION HELIXDIR DSSP CONFIGGRP BLASTOUT BLASTDB);
   my @ret ;
   print STDERR "=============\n" if($verbose);
   foreach my $var (@vars){
       my $v  = $ENV{$var} or die "$var not set";
       #print  STDERR " $var = $v\n " ;
       push @ret, $v ;
   }
   print STDERR  "\n===============\n" if($verbose);
   return @ret ;
}
sub util_sortuniqArray{
    my (@l) = @_ ;
    my $t = util_make_table(\@l);
    my @r = sort keys %{$t};
    return @r ;

}
sub util_CreateTableAndListFromSymbolSeperated{
    my ($str,$symbol) = @_ ;
    $str =~ s/$symbol/ /g;
    my @l = split " ", $str ;
    my $t = util_make_table(\@l);
    return ($t,@l) ;
}


sub util_GetStatOfList{
	my (@l) = @_ ;
	my @sort = sort {$a <=> $b} @l ;
	my $N = @l -1 ;
	my $first = $sort[0]; ;
	my $last = $sort[$N]; ;
	my $range = $last - $first;
	my ($mean,$sd) = util_GetMeanSD(\@l);
	return ($N+1,$first,$last,$range,$mean,$sd);
}





sub util_GetEnv{
    my ($l) = @_ ;
    my $ret = $ENV{$l} or die "Need to set environment variable $l"; 
    return $ret ;
}

sub util_make_list{
    my (@l) = @_ ;
    return \@l ;
}

sub util_debug_printNList{
	my ($str,@list) = @_ ;
	my $N = @list ;
	print "$str $N\n";
}

# make table from list
sub util_maketablefromlist{
    my ($l) = @_ ;
    my $t = {};
	my @l = @{$l} ;
	foreach my $i (@l){
		$t->{$i} = 0  if(! exists $t->{$i});
		$t->{$i}  =  $t->{$i} + 1 ;
	}
    return $t ;
}
        
sub util_make_table{
    my ($l) = @_ ;
    my $t = {};
	my @l = @{$l} ;
	foreach my $i (@l){
		$t->{$i} = 0  if(! exists $t->{$i});
		$t->{$i}  =  $t->{$i} + 1 ;
	}
    return $t ;
}
sub util_table_print{
    my ($table,$infile) = @_ ;
    my $N = keys %{$table};
    if(defined $infile){
       my $ofh = util_write($infile);
       foreach my $k (sort { $a <=> $b} keys %{$table}){
           print $ofh "{$k}\n";
       }
       close($ofh);
    }
    return ;

    print "There are $N entries \n";
    print "===========\n";
    foreach my $k (sort { $a <=> $b} keys %{$table}){
        print "Tableprint $k $table->{$k} \n";
    }
}

sub util_table_merge{
    my ($t1,$t2) = @_ ;
    my $t = $t1 ;
    foreach my $k (keys %{$t2}){
        my $v = $t2->{$k} ;
        $t->{$k} = $v ;
    }
    my $N = (keys %{$t});
    return ($t,$N) ;
}


sub util_mapFullLinetofirst{
	my ($fname) = @_ ;
	my $fh = util_read($fname);
	my $table = {};
	while(<$fh>){
	   next if(/^\s*$/);
	   next if(/^\s*#/);
	   chomp ;

		my (@l) = split ;
		my $first = shift @l ;
		$table->{$first} = $_;
	}
	close($fh);
	return $table ;
}
sub util_maketablefromfile{
    my ($file) = @_ ;
    my $fp = util_read($file);
    my $table ;
    while(<$fp>){
         next if(/^\s*$/);
         next if(/^\s*#/);
         my ($a,$b) = split ; 
         $b = 1 if(!defined $b);
         $table->{$a}= $b ;
    }
    close($fp);
    my $N = (keys %{$table});
    return ($table,$N) ;

}

## SAME as util_mapFullLinetofirst
sub util_maketablefromfile_firstentry{
    my ($file) = @_ ;
    my $fp = util_read($file);
    my $table ;
    while(<$fp>){
         next if(/^\s*$/);
         next if(/^\s*#/);
		 chomp ;
         my ($a,$b) = split ; 
         die "$file has no second" if(!defined $b);
         $table->{$a}= $_ ;
    }
    close($fp);
    my $N = (keys %{$table});
    return ($table,$N) ;

}


sub util_table_diff{
    my ($t1,$t2) = @_ ;
    my $common = 0 ;
    my $inAbutnotinB = 0 ;
    my $inBbutnotinA = 0 ;
    my $x ;
	my @common ;
	my @inAbutnotinB ;
	my @inBbutnotinA ;
    foreach my $k (keys %{$t1}){
        $x = exists $t2->{$k}? push @common, $k : push @inAbutnotinB, $k;
    }
    foreach my $k (keys %{$t2}){
        $x = exists $t1->{$k}? $x++ : push @inBbutnotinA, $k ;
    }

    return (\@common,\@inAbutnotinB,\@inBbutnotinA);
}


sub util_two2thepowerof{
     my ($number) = @_;
     my $ret = 1; 
     foreach my $i (1..$number){
        $ret = 2*$ret;
     }
    return $ret ;
}

sub util_round {
    my($number) = shift;
    return int($number + .5);
}

sub util_log_base {
    my ($base, $value) = @_;
    if($value < $base){
        return 1;
    }
     my $val =  log($value)/log($base);
     my $finalval =ceil($val);
     return $finalval;
}

sub util_ceil{
    my ($value) = @_;
     my $finalval =ceil($value);
     return $finalval;
}
sub util_floor{
    my ($value) = @_;
     my $finalval =floor($value);
     return $finalval;
}

sub util_get_pwd{
   my $PWD = getcwd ;
   return $PWD ;
}


sub util_execute_chop{
  my ($exec,$fname) = @_;
  my  $ret = ` $exec $fname `;
  chomp $ret ;
  $ret ;
}

sub util_printAndDo{
   my ($what,$dry) = @_ ;
   my $comment = defined $dry ? " Will run ( this is dry run ) " : "Running" ;
   print STDERR "$comment $what ...\n";
   system($what) if(!defined $dry);
}
sub util_printAndWrite2Script{
   my ($what,$fh) = @_ ;
   ## 
   #$what =~ s/\"/\\\"/g ;
   croak " undefined file handle " if(!defined $fh);
   print STDERR "$what ...\n";
   print $fh "$what ; ";
}


sub util_exit_if_doesnt_exist{
  my ($fname) = @_;
  croak "File $fname does not exist. Quitting " if(!-e $fname);
}

sub util_get_tag{
  return "" if(!-e "CVS/Tag");
  my $fh = new FileHandle("CVS/Tag",O_RDONLY) or croak ;
  my $tag = <$fh>;
  $tag =~ s/^T//;
  $tag ;
}

sub util_print_n_log {
    my ($ofh,$msg) = @_ ;
    print $msg ;
    print $ofh $msg ;
}

sub util_get_user{
   my $user = `whoami` ;
   chomp $user ;
   $user ;
}
sub util_get_cmdline{
    my ($exec , $list ) = @_ ;
    map { $exec = $exec . " $_ " ; } @$list ;
    $exec ;

}

sub util_percentage{
   my ($a,$b,$justval) = @_ ;  
   croak "a is undefined" if(!defined $a);
   croak "b is undefined" if(!defined $b);
   return $a if(defined $justval);
   $a = 1 if($a eq 0 || !defined $a);
   $b = 1 if($b eq 0 || !defined $b);
   my $percent = ($a - $b)/$a ;  # Changing the diff so that we are in sync with harness eqn
   #my $percent = 1 - ($a/$b);
   #my $percent = ($a/$b)-1 ; 
   #int($percent*100) . "%" ;
   int($percent*100);
}

sub util_is_integer {
   defined $_[0] && $_[0] =~ /^[+-]?\d+$/;
}

sub util_is_float {
      defined $_[0] && $_[0] =~ /^[+-]?\d+(\.\d+)?$/;
}


sub util_read_list_numbers{
    my @list = ();
    my ($file) = @_ ;
    my $fp = new FileHandle($file,O_RDONLY) or croak " Error for file: $file $!" ;
    while(<$fp>){
         next if(/^\s*$/);
         s/\s*//g;
         chomp ;
         croak " just expect a number " if(!/^[+-]?\d+(\.\d+)?$/);
         my ($num) = $_ ;
         push @list, $num ;
    }
    return @list ;
}


sub util_read_list_words{
    my @list = ();
    my ($file) = @_ ;
    my $fp = new FileHandle($file,O_RDONLY) or croak " $!" ;
    while(<$fp>){
         next if(/^\s*$/);
         next if(/^\s*#/);
         chomp ;
         my @l = split ; 
         push @list, @l ;
    }
    return @list ;

}
sub util_read_list_sentences{
    my @list = ();
    my ($file) = @_ ;
    my $PWD = getcwd ;
    my $fp = new FileHandle($file,O_RDONLY) or croak " util_read_list_sentences :  $! $file $PWD" ;
    while(<$fp>){
         next if(/^\s*$/);
         chomp ;
         push @list, $_ ;
    }
    return @list ;

}

sub util_read_list_firstitem{
    my @list = ();
    my ($file) = @_ ;
    my $PWD = getcwd ;
    my $fp = new FileHandle($file,O_RDONLY) or croak " util_read_list_sentences :  $! $file $PWD" ;
    while(<$fp>){
         next if(/^\s*$/);
         chomp ;
         my @l = split ;
         push @list, $l[0] ;
    }
    return @list ;

}

sub util_get_abs_path{
    my ($fname) = @_ ; 
    my $pwd = util_get_pwd();
    return $fname if($fname =~ /^\s*\// || $fname =~ /^\s*\~/);
    my $abspath = $pwd . "/" . $fname ;
    return $abspath ;
}

sub util_get_time_from_string{
    my ($str) = @_ ;
    my ($a,$b,$c,$d) = ($str =~ /(\d+)\.(\d+)u\s*(\d+)\.(\d+)s/);
    return $a + $c ;

}

sub util_isDesignSuccessful{
    my ($logname) = @_ ;
    my $fp = new FileHandle($logname,O_RDONLY) or croak " $! $logname" ;
    my $foundSynthesized = "FAILED";
    while(<$fp>){
        if(/Finished synthesizing design/){
            $foundSynthesized = "PASSED";
        }
    }
    return $foundSynthesized;
}

sub util_parseTimelog {
   my ($logname) = @_ ;
   my $fp = new FileHandle($logname,O_RDONLY) or croak " $logname $!" ;
   my $lastline;
   while(<$fp>){
     $lastline = $_;
   }
   return util_get_time_from_string($lastline);
}

sub util_open_or_append{
     my ($outfile)= @_;
     croak "Blank file name " if($outfile =~ /^\s*$/);
     my $fh ;
     if(! -e $outfile){
         $fh = util_write($outfile);
     }
     else{
         $fh = util_append($outfile);
     }
     return $fh ;
}

sub util_append{
     my ($outfile)= @_;
     die "Blank file name " if($outfile =~ /^\s*$/);
     #unlink $outfile;
     my $fh = new FileHandle($outfile,O_WRONLY|O_APPEND) or croak " could not write file $outfile as $!" ;
     return $fh ;
}

sub util_write{
     my ($outfile)= @_;
     croak "not defined" if(!defined $outfile);
     unlink $outfile;
     my $fh = new FileHandle($outfile,O_CREAT|O_WRONLY) or croak " could not write file $outfile as $!" ;
     #print "Writing to $outfile\n";
     return $fh ;
}


sub util_writelist2file{
     my ($outfile,@list)= @_;
     my $fh = util_write($outfile);
     foreach my $x (@list){
         print $fh "$x\n";
     }
	 close($fh);
}

sub util_read{
     my ($outfile)= @_;
     my $fh = new FileHandle($outfile,O_RDONLY) or croak " could not read file $outfile as $!";
     return $fh ;
}

sub util_wget{
     my ($file)= @_;
     util_printAndDo("wget --no-proxy $file");
}

sub util_makeCSH{
     my ($ofh)= @_;
     print $ofh "#!/bin/csh -f\n";
}

                                                                                                                                                             
sub util_pick_random_from_list{
    my ($list) = @_ ;
    my @temporaries = @{$list};
    my $num = @temporaries ;
    my $r = floor($num*rand());      
    my $operator = $temporaries[$r] or croak ;
    return $operator ;
}

sub util_pick_n_random_from_list{
    my ($list,$n) = @_ ;
    my @l ; 
    my $done ; 
    my $N = @{$list} ; 
    foreach my $i (1..$n){
        my $over = 0 ; 
        while (!$over){
            my $x = util_pick_random_from_list($list);
            #print "$x $i $N\n";
            if(! exists $done->{$x}){
                push @l, $x ; 
                $over =  1 ; 
                $done->{$x} = 1 ; 
            }
        }
    }
    return \@l ; 
}




sub util_enter_maxcutoff{
    print STDERR " Enter max cutoff number \n";
    my $qnum = <> ;
    chomp $qnum ; 

    my $ret = 1000000 ;
    $ret = $qnum if($qnum !~ /^\s*$/);
    return $ret ;
}


sub util_EnterName{
    my ($default) = @_ ; 
    if(defined $default){
        print "Press enter to choose default name $default\n";
    }
    my $qname = <> ;
    chomp $qname ; 
    if(defined $default){
        if($qname =~ /^\s*$/){
            return $default ; 
        }
    }
    my @l = split " ",$qname; 
    if(@l != 1){
        print "Warning: Need a single name. Try again\n";
        return util_EnterName(); 
    }
    return $l[0]; 
}

sub util_EnterNumber{
    my $qnum = <> ;
    chomp $qnum ; 
    my @l = split " ",$qnum; 
    if(@l != 1 || !util_is_integer($l[0])){
        print "Warning: Need a single number. Try again\n";
        return util_EnterNumber(); 
    }

    return $l[0]; 
}
sub util_EnterTwoNumbers{
    my $qnum = <> ;
    chomp $qnum ; 
    my @l = split " ",$qnum; 
    if(@l != 2 || !util_is_integer($l[0]) || !util_is_integer($l[1])){
        print "Warning: Need 2 numbers. Try again\n";
        return util_EnterTwoNumbers(); 
    }

    return @l ;
}
sub util_EnterTwoNames{
    my $qname = <> ;
    chomp $qname ; 
    my @l = split " ",$qname; 
    if(@l != 2){
        print "Warning: Need 2 Names. Try again\n";
        return util_EnterTwoNames(); 
    }

    return @l ;
}


sub util_AddBeforeEach{
    my ($expr,@l) = @_ ; 
    my $str  = join " -$expr ",@l ;
    $str  = " -$expr $str";
    return $str ; 
}

sub util_ignoreIfDistanceIsLessthan{
    my ($num,@l) = @_ ; 
    my $iter = combinations(\@l, 2);
    while (my $c = $iter->next) {
        my @combo = @{$c} ; 
        if(abs($combo[0] - $combo[1]) < $num){
            return 1 ; 
        }
    }
    return 0 ;
}
sub util_ignoreIfAnyAreEqual{
    my (@l) = @_ ; 
    my $iter = combinations(\@l, 2);
    while (my $c = $iter->next) {
        my @combo = @{$c} ; 
        if($combo[0] == $combo[1]){
            return 1 ; 
        }
    }
    return 0 ;
}

sub util_split_list_numbers{
    my ($l,$point) = @_ ; 
    my @l = @{$l} ; 

    my @lte ;
    my @gt ;
    foreach my $e (@l){
        if($e <= $point){
            push @lte, $point ;  
        }
        else{
            push @gt, $point ;  
        }
    }
    my $n1 = @lte ;
    my $n2 = @gt ;
    return (\@lte, \@gt,$n1,$n2);
}


sub util_read_Mapping_PDB_2_SWISSPROT{
    my ($infile) = @_;
    my $info = {};
    my $ifh = util_read($infile);
    print "Reading mapping $infile\n";

    #ignore first two lines
    <$ifh>;
    <$ifh>;
    my $uniqueEC ; 
    my $uniqueSP ; 
    while(<$ifh>){
         next if(/^\s*$/);
         s/|//g;
         s/,//g;
         s/,//g;
         my ($nm,$j1,$chainid,$j2,$n1,$j3,$n2,$j4,$swissprot,$jjj,$ec) = split ;
         #print " $swissprot $ec\n";
    
             if($chainid eq "A"){
                if(!defined $info->{$nm}){
                   $info->{$nm} = {} ;
               }
                   #print "$nm $chainid $swissprot  \n";
                 $info->{$nm}->{SWISSPROT} = $swissprot ;
                 $info->{$nm}->{EC} = [] if(!defined $info->{$nm}->{EC});
                 $info->{SWISSPROT2EC}->{$swissprot} = [] if(!defined $info->{SWISSPROT2EC}->{$swissprot}) ;
                 if($ec ne "0.0.0.0"){
                    push @{$info->{$nm}->{EC}}, $ec ; 
                    push @{$info->{SWISSPROT2EC}->{$swissprot}},$ec ;
                    $uniqueSP->{$swissprot} = 1 ; 
                    $uniqueEC->{$ec} = 1 ; 
                 }
            }
    }
    close($ifh);
    return ($info,$uniqueEC,$uniqueSP) ; 
}

sub util_getECfromPDB{
    my($info,$nm) = @_ ; 
    $nm = lc($nm);
    my $ec =  $info->{$nm}->{EC} ; 
    return undef if(!defined $ec);
    return $ec ; 
}


sub util_filter_basedon_EC{
    my ($info,$uniqueEC,$uniqueSP,$list,$ofh,$ignore0000) = @_ ; 
    my @list = @{$list};
    my $N = @list  ;
    my $ecdone = {};
    my $spdone = {};
    my $pdbTable = {};
    my $cnt = 0 ;
    my $ignored = util_write("ignored");
    foreach my $UCPDBID  (@list){

       my $pdbid = lc($UCPDBID);
       if(!defined $info->{$pdbid}){
             print $ignored "PDB $pdbid not found in mapping \n";
          next ;
       }
       my @ec = @{$info->{$pdbid}->{EC}};
       my $sp = $info->{$pdbid}->{SWISSPROT};
       my $N = @ec ; 
       if($N > 1){
             print $ignored "PDB $pdbid ignored as it has more than 1 :$N ECS\n";
          next ;
       }
       if(@ec == 0){
            my @ECCC = @{$info->{SWISSPROT2EC}->{$sp}}; 
            if(@ECCC == 0){
                print $ignored "Ignoring pdb $pdbid with swiss $sp as there are no ECS\n";
                next ; 
            }
               print "Did not find EC number for $pdbid\n";
            @ec = @ECCC ;
       }
          print $ofh "$UCPDBID\n" if(defined $ofh);
       $spdone->{$sp} = $pdbid ; 

       my $added = 0 ; 
       print "did not find EC for $pdbid \n" if(@ec == 0);
       foreach my $ec (@ec){
                 #if(exists $ecdone->{$ec}){
                    #print $ignored "$ec exists already for $pdbid $ecdone->{$ec} \n";
                  #next ;
              #}
              $added = 1 ;
              $ecdone->{$ec} = [] if(!defined $ecdone->{$ec});
              push @{$ecdone->{$ec}}, $pdbid ; 
        }

        $pdbTable->{$pdbid} = $sp ; 
        $cnt++ if($added);

   }
   print "Wrote ignored PDBS in ignored\n";

   return ($ecdone,$spdone,$pdbTable,$cnt); 
}


sub util_getPDBID_basedon_SP{
    my ($info,$uniqueEC,$list,$id) = @_ ; 
    my ($ecdone,$spdone,$pdbTable,$cnt)= util_filter_basedon_EC($info,$uniqueEC,$list);
    if(exists $pdbTable->{$id}){
        return $id ; 
    }
    else{
       my $sp = $info->{$id}->{SWISSPROT};
       if(exists $spdone->{$id}){
           return $spdone->{$id}  ;  
       }
       else{
           return undef ; 
       }
    }
}


###########################################################################
# Parses a file and returns blocks between start and end
###########################################################################

sub util_ParseBlocks{
    my ($infile,$start,$end) = @_;
    my @blocks ;
    my $ifh = util_read($infile);

    while(<$ifh>){
         if(/$start/){
             my @block ;
            push @block,$_;
            while(!/$end/){
                $_ = <$ifh>;
                if(!$_){
                    print "NULL - hnece lasting\n";
                     last ;
                }
                if(/$start/){
                    @block = ();
                    push @block,$_;
                }
                else{
                    push @block,$_;
                }

            }
            push @blocks, \@block ;
         }
    }
    close($ifh);
    return \@blocks ;
}

sub util_ParseBlockForString{
    my ($block,$str) = @_;
    my @lines ;
    foreach my $i (@{$block}){
        $_ = $i ;
         if(/$str/){
            push @lines,$_;
         }
    }
    return \@lines ;
}



#############################################################
## Calculate the mean and SD of a matrix, columnwise 
## Also takes a range, and tells how many lie outside that range.
## So this range ideally makes sense for one column only
#############################################################
sub util_ProcessRowAndColumnsForMean{

   my ($rows,$top,$low)= @_ ;
   my @rows = @{$rows};
   my $nrows = @rows - 1;
   my @cols ; 


    #print "Number of rows -1 = $nrows\n";
    my $once = 1 ;
    foreach my $i (0..$nrows){
        my $row = $rows[$i];
        my @row = @{$row};
        my $ncol = @row - 1 ;
        if($once){
           print "Number of column -1 = $ncol\n";
           $once = 0 ;
        }
        foreach my $j (0..$ncol){
            if(!defined $cols[$j]){
                my @l = ();
                $cols[$j] = \@l;
            }
            push @{$cols[$j]}, $row[$j];
        }
    }
    
    my @means ; 
    my @sds ; 
    my $NN = @cols ;
    foreach my $container (@cols){
        my $mean = Math::NumberCruncher::Mean($container) or warn "Mean not found" ;
        my $sd = Math::NumberCruncher::StandardDeviation($container) or warn "sd not found" ;
        next if(!defined ($mean && $sd));
        push @means, util_format_float($mean,1) ;
        push @sds, util_format_float($sd,1) ;
        my $incnt = 0 ;
        my $outcnt = 0 ;
        foreach my $i (@{$container}){
            $i = -($i);
            if($i > $low && $i < $top){
                $incnt++;
            }
        else{
                $outcnt++;
            }
        }
        print "mean = $mean sd = $sd \n";
        print "incnt = $incnt outcnt = $outcnt \n";
    }
    return (\@means,\@sds);

}


sub util_format_float{
    my ($d,$v) = @_; 
    croak if(!defined $d);
    if (defined $v && $v eq 3) {return sprintf("%8.3f", $d);}
    if (defined $v && $v eq 1) {return sprintf("%8.1f", $d);}

    return sprintf("%8.3f", $d); 
}

sub util_mysubstr {
    my ($len,$str,$start,$end) = @_ ; 
    if($len < $start ){
        die "dying $str" if(!($str =~ /TER/));
        return "" ;
    }
    my $diff = $end == 1 ? $end :  $end - $start + 1 ; 
    $start-- ; 
    #print "($str,$start,$diff)\n";
    my $s =  substr($str,$start,$diff);
    $s =~ s/\s*//g ; 
    return $s ; 
}

sub util_mysubstrDontStripSpace {
    my ($len,$str,$start,$end) = @_ ; 
    if($len < $start ){
        die "dying $str" if(!($str =~ /TER/));
        return "" ;
    }
    my $diff = $end == 1 ? $end :  $end - $start + 1 ; 
    $start-- ; 
    return substr($str,$start,$diff);
}


sub ParseAPBSResult{
   my ($size,$infile,$listfile) = @_ ;
   my $ifh = util_read($infile);
   my $list = {};
   if(defined $listfile){
       my @list= util_read_list_sentences($listfile);
       map { $list->{$_} = 1 ; } @list ;
   }
   

   my $blocks = util_ParseBlocks($infile,"Starting read","Ending read");
   my @blocks = @{$blocks};

   my $len = @blocks ; 
   print STDERR "\%There were $len blocks\n";




    my @diffs = ();
    foreach my $block (@blocks){

    chomp $block ;
        my $lines = util_ParseBlockForString($block,"Resultfile is");
        if(@{$lines} == 0){
            die ;
        }
        my $line = shift @{$lines}; 
        my ($pdb)  = ($line =~ /\/(....)\.pdb.out/);
        ($pdb)  = ($line =~ /(....)\.pdb.out/) if(!defined $pdb);
        if(defined $listfile && defined $pdb ){
             if(!exists $list->{$pdb}){
                print "nexting\n";
                #return undef ;
                 next;
             }
        }
    
    
        $lines = util_ParseBlockForString($block,"potential");
        if(@{$lines} != $size){
            print "nexting as size $size doestn match\n";
                return undef ;
            next ;
        }
    
       my @vals  ; 
       my $somethingwrong = 0 ; 
       while($line = shift @{$lines}){
            my ($val) = ($line =~ /potential\s*=\s*(.*)/);
            if(!defined $val || $val =~ /^\s*$/){
                $somethingwrong = 1 ;
                last ; 
            }
            #print $ofh "$val ";
            push @vals, $val;
        }
    
        next if($somethingwrong);
    
    
        my $iter = combinations(\@vals, 2);
        my @diff  ;
        while (my $c = $iter->next) {
                my @combo = @{$c} ; 
                my ($a,$b) = @combo ; 
                my $d= $a - $b ;
                #print " $d= $a - $b ; \n";
                $d=util_format_float($d,1);
                push @diff , $d; 
        }
        push @diffs , \@diff ;
    }
    
    #my ($means,$sds) = util_ProcessRowAndColumnsForMean(\@diffs,250,150);
    #util_PrintMeanAndSD($ofh,$means,$sds);


    return \@diffs ;
}

sub util_PrintMeanAndSD{
    my ($ofh,$means,$sds) = @_ ;
    #print $ofh "\\rowcolor{orange} \n";
    print $ofh "Mean & " ;
    util_printTableLine($ofh,$means);
    #print $ofh "\\rowcolor{orange} \n";
    print $ofh "SD & ";
    util_printTableLine($ofh,$sds);
}


sub util_printResult{
    my ($ofh,$bestScore,$bestList,$cnt) = @_ ; 
    carp "Best score not defined" if(! defined $bestScore);
    print $ofh  "#RESULT $cnt  SCORE - $bestScore\n";
    print  $ofh "# ";
    foreach my $atom (@{$bestList}){
          my ($res,$num,$type) = split "/", $atom ;
          print $ofh  "-  $atom ";
    }
    print  $ofh "\n";
}


sub util_printTablePre{
    my ($ofh,$caption) = @_ ; 
    $caption = "XXX" if(!defined $caption);
    print $ofh "\\begin{center} \n";
    print $ofh "\\begin{table*} \n";
   print $ofh "\\caption { $caption }  \n";
    #print $ofh "\\rowcolors{1}{tableShade}{white} \n";
    
    print $ofh "\\begin{tabular}{ |c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c            } \n";
    #print $ofh "\\rowcolor{red!50}  \n";
}

sub util_printTablePost{
    my ($ofh,$caption) = @_ ; 
    $caption = "XXX" if(!defined $caption);

   #print $ofh "\\rowcolor{green!50}  \n";
   #print $ofh "\\rowcolor{orange!50}  \n";
   print $ofh "\\end{tabular}  \\label{label} \n";
   #print $ofh "\\caption {\\bf $caption }  \n";
   print $ofh "\\end{table*} \n";
   print $ofh "\\end{center} \n";
}

sub util_printTableLine{
    my ($ofh,$list,$NNN) = @_ ; 
    my @l ; 
    map { 
       if(util_is_float($_)){
             $_ = util_format_float($_,$NNN);
       }
       push @l, $_ ; 
    } @{$list} ;

    my $str = join " & ", @l;
    print $ofh " $str \\\\ \n";
    #print $ofh "\\hline \n";
}

sub util_readAPBSPotentialFromStart{
    my ($protein,$apbsdir) = @_ ; 
    my $pqrfile = "$apbsdir/$protein/$protein.pqr";
    my $pqr = new PDB();
    $pqr->ReadPDB($pqrfile,"hackforpqr");
    my @pots = ();
    my $potential = "$apbsdir/$protein/pot1.dx.atompot";
    util_readAPBSPotential(\@pots,$potential);
    return ($pqr,\@pots);
}

sub util_readAPBSPotential{
    my ($pots,$potential) = @_ ; 
    my $ifh = util_read($potential);
    print "util_readAPBSPotential : Reading file $potential \n";
    while(<$ifh>){
         next if(/^\s*$/);
         chomp ;
         my @l = split ",",$_;
         #print "$l[3] \n";
         push @{$pots}, $l[3] ;
    }
    close($ifh);
    
}

sub util_CmdLine{
    my ($nm,$var) = @_ ;
    #$var = ${$nm} ;
    die "Error: In command line parsing. Needed command line ==> $nm" if(!defined $var);
}

sub util_parsePDBSEQRESNEW{
    my ($infile,$all,$writedata,$rev,$FASTADIR) = @_ ; 
	die if(! $infile && $all && $writedata && $rev && $FASTADIR);
    my $ifh = util_read($infile);
    #my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT) = util_SetEnvVars();

    my $seenNm ; 
    my $info  = {}; 
    my $infoSeq2PDB = {} ; 
    my $mapChainedName2Name = {} ; 
    my $uniqueseqs = {} ; 
    
    my $nm ;
    my $CNT = 0 ;
    while(<$ifh>){

         next if(/^\s*$/);
         next if(/^\t*$/);
         if(/^>/){
             s/^>//;
             my @l = split ;
             $nm = $l[0];
            $nm = uc($nm);
            die if(!defined $nm); 
            warn "Repeat $nm " if(exists $info->{$nm});
            $CNT++;
            $info->{$nm} = "" ;
            next ;
        }


        ## add to seq
        die "What is this $_"  if(/>/);
        if(!defined $_){
            die "$nm $CNT " ;
        }
        chomp;
        die $_ if(!defined $nm);
		s/\s*//g;
	#print "KKKKKKKKKK $infile\n";
        $info->{$nm} = $info->{$nm} . $_ ;

     }
     print "There were $CNT seqeunces to start with\n";

     foreach my $k (keys %{$info}){
         my $seq = $info->{$k};
        my $length = length($seq);
        #print "$length $k \n";
        $infoSeq2PDB->{$seq} = [] if(!defined $infoSeq2PDB->{$seq}) ;
        push @{$infoSeq2PDB->{$seq}}, $k ;
     }

     my $N = (keys %{$infoSeq2PDB});
	 my $DIFFNSEQ = $CNT - $N ;
     print "There are $N unique (diff = $DIFFNSEQ) \n";

     if(defined $writedata && $writedata){
         my $ofh = util_write("list.unique");

         my $SORTING = {};
         foreach my $k (keys %{$infoSeq2PDB}){
             my $len = length($k);
            $SORTING->{$k} = $len ;
         }


         my @sl = sort { $SORTING->{$b} <=> $SORTING->{$a}} (keys %{$SORTING});
         my $CNTwrite = 0 ;
         foreach my $k (@sl){
            my @pdbswithseq  = sort @{$infoSeq2PDB->{$k}} ;
             my $len = length($k);
            #if($CNTwrite < 3 ){
            if(1){
               #print "$len \n";
               foreach my $pdb (@pdbswithseq){
                   ($pdb) = split " ", $pdb ;
                   #print "KKKKKKKK $pdb\n";
                   my $fastanm = "$FASTADIR/$pdb.ALL.1.fasta";
                   $CNTwrite++;
                   my $FH = util_write($fastanm);
                  print $FH "\>$pdb\n";
                  print $FH "$k\n";
                  close($FH);
				  if($rev){
                      my $fastanm = "$FASTADIR/$pdb.ALL.1.fasta.comp.fasta";
					  my $rev = util_getComplimentaryString($k) ;
                      my $FH = util_write($fastanm);
                      print $FH "\>${pdb}_rev\n";
                      print $FH "$rev\n";
                      close($FH);
				  }
       
               }
               my $v = $pdbswithseq[0];
                print $ofh "$v\n";
            }
         }
         print "Wrote $CNTwrite finally\n";
             close($ofh);


         $ofh = util_write("list.allpdbschained");
         foreach my $k (keys %{$mapChainedName2Name}){
             print $ofh "$k\n";
         }
         close($ofh);

         $ofh = util_write("list.allpdbs");
         foreach my $k (keys %{$info}){
             $k = uc($k);
             print $ofh "$k\n";
         }
         close($ofh);
     }


     return ($info,$infoSeq2PDB,$mapChainedName2Name) ;
}


sub util_parsePDBSEQRESOLLLLLL{ die ; }
  



sub parseSingleLinePDBSEQRES{
    my ($line) = @_ ; 
    my ($nm,$chain,$type,$len,$fullnm) = ($line =~ /^.(....).(.)\s*mol:(\w+)\s*length:(\d+)\s*(.*)/);
    return ($nm,$chain,$type,$len,$fullnm) ;

}
sub util_Annotate{
    my ($file) = @_ ; 
    my ($outfile) = "$file". ".annotate";
    util_printAndDo("annotate.pl -in ~/pdb_seqres.txt -lis $file -out $outfile -cutoff 100 -anndis 10 ");

}

sub util_IsZero{
    my ($num) = @_ ; 
    if(abs($num) < $EPSILON){
        return 1 ; 
    }
    else{
        return 0 ; 
    }

}

sub util_Banner{
    my ($str) = @_ ; 
    print STDERR "===============================================\n";
    print STDERR "$str\n";
    print STDERR "===============================================\n";

}
sub util_PrintInfo{
    my ($str) = @_ ; 
    print STDERR "Info: $str\n";

}


sub util_getTmpFile{
    my $time = int(time() *rand());
    my $tmpfile =  "sandeeptmp.$time";
    return $tmpfile ;
}

sub util_percentages{
    my (@values) = @_ ;
    my $sum = 0 ; 
    #print " KKKKKKKKKKKKK " , @values , "\n";
    foreach my $v (@values){
       $sum = $sum + $v ; 
    }
    my @l ; 
    foreach my $v (@values){
       push @l, util_format_float(($v*100)/$sum,3) ; 
    }
    return @l ;
}



sub util_readPeptideInfo{
   my ($info,$nm,$infile) = @_ ;
   $info->{$nm} = {};
   die "ERROR: util_readPeptideInfo $infile does not exist" if(! -e $infile);
   my $ifh = util_read($infile);
   while(<$ifh>){
        next if(/^\s*$/);

        if(/Average Residue Weight/){
           my ($charge) = (/Charge\s*=\s*(.*)/) or die;
           $info->{$nm}->{CHARGE} = $charge ;
       }
       if(/Residues = /){
               my ($nres) = (/Residues = (\d+)/) or die;
            $info->{$nm}->{NRES} = $nres ;
       }
        if(/^(Molecular|Basic|Acidic|Polar)/i){
           my (@l) = split ;
           my $N = @l -1 ;
           $info->{$nm}->{$l[0]} = $l[$N]; ;
        }
   }
   if(defined $info->{$nm}->{Acidic} && defined $info->{$nm}->{Basic}){
       $info->{$nm}->{AcidBasic} = $info->{$nm}->{Acidic} + $info->{$nm}->{Basic} ;
       $info->{$nm}->{PAB} = $info->{$nm}->{Acidic} + $info->{$nm}->{Basic} + $info->{$nm}->{Polar} ;
   }
   close($ifh);
}


sub util_ReadLine{
           my ($LINE) = @_ ;
           my $len = length($LINE);
           my ($atomstr , $serialnum , $atomnm , $alt_loc , $resname , $chainId , $resnum , $codeforinsertion , $x , $y , $z );

 #1 -  6        Record name     "ATOM  "                                            
 #7 - 11        Integer         Atom serial number.                   
 #13 - 16        Atom            Atom name.                            
 #17             Character       Alternate location indicator.         
 #18 - 20        Residue name    Residue name.                         
 #22             Character       Chain identifier.                     
 #23 - 26        Integer         Residue sequence number.              
 #27             AChar           Code for insertion of residues.       
 #31 - 38        Real(8.3)       Orthogonal coordinates for X in Angstroms.                       
 #39 - 46        Real(8.3)       Orthogonal coordinates for Y in Angstroms.                            
 #47 - 54        Real(8.3)       Orthogonal coordinates for Z in Angstroms.                            
 #55 - 60        Real(6.2)       Occupancy.                            
 #61 - 66        Real(6.2)       Temperature factor (Default = 0.0).                   
 #73 - 76        LString(4)      Segment identifier, left-justified.   
 #77 - 78        LString(2)      Element symbol, right-justified.      
 #79 - 80        LString(2)      Charge on the atom.       

           $atomstr = util_mysubstr($len,$LINE,1 ,  6);
           $serialnum = util_mysubstr($len,$LINE,7 , 11);
           $atomnm = util_mysubstrDontStripSpace($len,$LINE,13 , 16);
           $alt_loc = util_mysubstr($len,$LINE,17,1);
           $resname = util_mysubstr($len,$LINE,18 , 20);
           $chainId = util_mysubstr($len,$LINE,22,1);
           $resnum = util_mysubstr($len,$LINE,23 , 26);
           $codeforinsertion = util_mysubstr($len,$LINE,27,1);
           $x = util_mysubstr($len,$LINE,31 , 38);
           $y = util_mysubstr($len,$LINE,39 , 46);
           $z = util_mysubstr($len,$LINE,47 , 54);
		   #print " $x $y $z lllllllllllllllllllllll\n";
		   my $singlespace = " ";
		  if($atomstr =~ /ATOM/){
		  	 $atomstr = "ATOM  ";
		  }
          #printf( "%6s%5s %4s%1s%3s %1s%4s%1s   %8s%8s%8s\n", 
		        #$atomstr,$serialnum,$atomnm,$alt_loc,$resname,$chainId,$resnum,$codeforinsertion,$x,$y,$z);
          return ($atomstr , $serialnum , $atomnm , $alt_loc , $resname , $chainId , $resnum , $codeforinsertion , $x , $y , $z ); 
}



sub util_getECID{
    my ($ec,$level) = @_ ; 
    $ec =~ s/\./YYY/g ; 
    my @l = split "YYY",$ec ; 
    my $str = "";
    my $first = 1 ; 
    foreach my $i (1..$level){
        if(!$first){
            $str = $str . ".";
        }
        $str = $str . $l[$i -1];
        $first = 0 ; 
    }
    return $str ; 
}

sub util_uniq2 {
    my ($list) = @_ ;
    my %seen = ();
    my @r = ();
    foreach my $a (@{$list}) {
        unless ($seen{$a}) {
            push @r, $a;
            $seen{$a} = 1;
        }
    }
    return @r;
}
sub vSetVerbose{
    my ($val) = @_ ;
    my $ret = $ENV{VERBOSE};
    if(defined $ret){
        $ENV{VERBOSE} = $val ;
    }
}
sub vIncrVerbose{
    my $ret = $ENV{VERBOSE};
    if(defined $ret){
        $ENV{VERBOSE} = $ENV{VERBOSE} + 1 ;
    }
}
sub vDecrVerbose{
    my $ret = $ENV{VERBOSE};
    if(defined $ret){
        $ENV{VERBOSE} = $ENV{VERBOSE} - 1 ;
    }
}

sub vprint{
    my ($str) = @_ ;
    my $ret = $ENV{VERBOSE};
    if(defined $ret && $ret > 0 ){
        my $tab = "";
        foreach my $i (1..$ret){
            $tab = $tab. "\t";
        }
        print STDERR "$tab$str\n";
    }
}
sub vprintheader{
    my ($str) = @_ ;
    my $ret = $ENV{VERBOSE};
    if(defined $ret && $ret > 0){
        my $tab = "";
        foreach my $i (1..$ret){
            $tab = $tab. "\t";
        }
        print STDERR "\n\n======================================================================\n";
        print STDERR "$tab$str\n";
        print STDERR "======================================================================\n";
    }
}
sub util_copy_table{
    my ($table) = @_ ; 
    my $ret = {};
    foreach my $k (keys %{$table}){
        $ret->{$k} = $table->{$k}; 
    }
    return $ret; 
}
sub util_GetFastaFiles{
    my ($FASTADIR,@pdbs) = @_ ; 
    my $mapinfo ;
    my $fhnotfound = util_write("list.notfound");
    my @files = ();
    foreach my $i (@pdbs){
        $i = uc($i);
        $i =~ /^\s*$/;
        #my @f = <$FASTADIR/$i.*fasta>;
        my $j = "$FASTADIR/$i.ALL.1.fasta";
        if(!defined $j){
            warn "Did not find fasta file $j for pdb $i in dir $FASTADIR \n";
            print $fhnotfound "$i\n";
            #push @ignored, $i;
            next ;
        }
        print " pushing $j \n" ;
        push @files, $j;
        $mapinfo->{$j} = $i ;
    }
    return ($mapinfo,@files) ; 
}

sub util_WriteFastaFromAtoms{
    my ($pdb,$allatoms,$fastafh,$origpdb) = @_ ; 
    my @allresidues ;
    my $allres  = {};
    foreach my $atom (@{$allatoms}){
        my $nm = $atom->GetName();
        next if($nm =~ /HOH/);
        my $resnum = $atom->GetResNum();
        my ($res) = $pdb->GetResidueIdx($resnum);
        push @allresidues ,$res ; 
        $allres->{$res->GetResNum()} = $res ; 
    }
    my $allresiduesN = @allresidues ;
    my $s2 = "";
    my $s1 = "";
    foreach my $i (sort  { $a <=> $b }  keys %{$allres}){
        my $r = $allres->{$i} ; 
        my $s = $r->PrintSingleLetter($pdb);
        $s1 = $s1 .  "$s$i," ;
        $s2 = $s2 . "$s" ;
    }
    print $fastafh "\>$origpdb.$s1;\n";
    print $fastafh "$s2\n";
}


sub util_WriteFastaFromResidueNumbers{
    my ($pdb,$resnumbers,$fastafh,$origpdb) = @_ ; 
    my @allresidues ;
    my $allres  = {};
    foreach my $resnum (@{$resnumbers}){
        my ($res) = $pdb->GetResidueIdx($resnum);
        push @allresidues ,$res ; 
        $allres->{$res->GetResNum()} = $res ; 
    }
    my $allresiduesN = @allresidues ;
    my $s2 = "";
    my $s1 = "";
    foreach my $i (sort  { $a <=> $b }  keys %{$allres}){
        my $r = $allres->{$i} ; 
        my $s = $r->PrintSingleLetter($pdb);
        $s1 = $s1 .  "$s$i," ;
        $s2 = $s2 . "$s" ;
    }
    print $fastafh "\>$origpdb.$s1;\n";
    print $fastafh "$s2\n";
}

sub util_GetPotForAtom{
    my ($a,$pqr,$pots) = @_ ;
    croak "not defined" if(!defined $a);
    my $number = $a->GetResNum();
    #print STDERR "kkkkkkkkk $number \n";
    my $atomnm = $a->GetType();
    my ($aPqr) = $pqr->GetAtomFromResidueAndType($number,$atomnm) or croak "No PQR?" ;
    croak "could not find $number $atomnm" if(!defined $aPqr);

    my ($x,$y,$z) = $a->Coords();
    my ($x1,$y1,$z1) = $aPqr->Coords();
    if(1 && !util_IsZero($x-$x1+$y -$y1+$z-$z1)){
            # this will mismatch as rotation is done
             #warn "Warning: $x,$y,$z $x1,$y1,$z1 do not match" ; 
    }

    #my ($i1) = $a->GetIdx();
    my ($i2) = $aPqr->GetIdx();
    #imp -1 
    my $NPots = @{$pots};
    #$a->Print();
    my $pot = $pots->[$i2-1] or die "Expected to find potential for residue number $number i2=$i2 and $NPots";
    return $pot ; 
}

sub util_ReadPdbs{
    my ($PDBDIR,$APBSDIR,$readpotential,@P) = @_ ; 
    my @ret ; 
    die "Expected at least one protein" if(!@P);
    foreach my $p1 (@P){
         my $file1 = "$PDBDIR/$p1.pdb";
         my $PPP = new PDB();
         $PPP->ReadPDB($file1);
         print STDERR "Reading $p1\n";

         my $pqrfile = "$APBSDIR/$p1/$p1.pqr";
         my $pqr = new PDB();
         my @pots = ();
         if($readpotential){
              #$pqr->SetReadCharge(1);
              $pqr->ReadPDB($pqrfile,"hack");
              #$pqr->ReadPDB($pqrfile);
              my $potential = "$APBSDIR/$p1/pot1.dx.atompot";
              util_readAPBSPotential(\@pots,$potential);
         }

         

         my $info->{PDBNAME} = $p1 ;
         $info->{PDBOBJ} = $PPP ;
         $info->{PQR} = $pqr ;
         $info->{POTS} = \@pots ;
         push @ret , $info ; 
    }
    return @ret ;
}

sub util_WriteClustalAln{
    my ($protein1,$matches,$matchedProteins,$alnfh,$map,$annMap) = @_ ; 
    my @matchedProteins = @{$matchedProteins} ; 
    my $N = @matchedProteins ; 
    my $halfN = $N/2 ;
    my $quarterN = $halfN/2 ;
    my @matches = @{$matches} ; 
    my $first = 1 ; 
    my $EXTEND = {};
    while(@matches){
        my $protein2 = shift @matchedProteins ; 
        my $MATCH = shift @matches ; 
        if($first){
              print $alnfh "CLUSTAL 2.1 multiple sequence alignment\n\n\n";
              my $pdb1 = $map->{$protein1} ;
              my $mappedname = exists $annMap->{$protein1} ? $protein1. "." . $annMap->{$protein1} : $protein1 ;
              my $XXX =  sprintf ( "%-12s",  $mappedname); 
              print $alnfh "$XXX";
              my $CNT = 0 ; 
              foreach my $resnum (sort {$a <=> $b}  keys %{$MATCH}){
                   $CNT++ ; 
                 $EXTEND->{$CNT} = {};
                 $EXTEND->{$CNT}->{CNT} = 1 ; 
                 $EXTEND->{$CNT}->{MATCHES} = {};
                 my ($res) = $pdb1->GetResidueIdx($resnum);
                 my $s = $res->PrintSingleLetter($pdb1);
                 $EXTEND->{$CNT}->{MATCHES}->{$s} =  1 ; 
                 print $alnfh "$s";
              }
              print $alnfh "\n";
              $first = 0 ;
        }


        my $pdb2 = $map->{$protein2} ;
        my $mappedname = exists $annMap->{$protein2} ? $protein2. ".". $annMap->{$protein2} : $protein2 ;
        my $XXX =  sprintf ( "%-12s",  $mappedname); 
        print $alnfh "$XXX";
        my $CNT = 0 ; 
        foreach my $k (sort {$a <=> $b}  keys %{$MATCH}){
             $CNT++ ; 
             my $resnum = $MATCH->{$k} ;
             my $s = "-";
             if($resnum ne "-"){
                  $EXTEND->{$CNT}->{CNT} =  $EXTEND->{$CNT}->{CNT} + 1 ; 
                  my ($res) = $pdb2->GetResidueIdx($resnum);
                  $s = $res->PrintSingleLetter($pdb2);
                  $EXTEND->{$CNT}->{MATCHES}->{$s} =  1 ; 
             }
             print $alnfh "$s";
        }
        print $alnfh "\n";
    }

    foreach my $CNT (sort {$a <=> $b}  keys %{$EXTEND}){
        my $cnt = $EXTEND->{$CNT}->{CNT} ;
        if($cnt >= (3*$quarterN)){
           my $res = $EXTEND->{$CNT}->{MATCHES} ;
           my @keys = (keys %{$res});
           my $str = join "|" , @keys ;
           $str = "(" . $str . ")";
           print "$CNT $cnt $str ========= \n";
        }
    }

    
}

sub util_GetClosestAtoms_intwoPDBs{
    my ($p1,$p2,$PDBDIR,$maxdist) = @_ ; 
    my $file1 = "$PDBDIR/$p1.pdb";
    my $file2 = "$PDBDIR/$p2.pdb";
    my $pdb1 = new PDB();
    $pdb1->ReadPDB($file1);
    my $pdb2 = new PDB();
    $pdb2->ReadPDB($file2);
    
    my @reslist = $pdb1->GetResidues();
    
    my @atoms ;
    while(@reslist){
        my $r = shift @reslist ;
        my @aaa = $r->GetAtoms();
        push @atoms, @aaa;
    }
    
    my @resultsall ; 
    my $DONE ;
    foreach my $atom (@atoms){
         my $list = util_make_list($atom);
         my ($junk,$neighatoms)  = $pdb2->GetNeighbourHoodAtom($list,$maxdist);
         foreach my $r (@{$neighatoms}){
               my $atomstr = $r->GetAtomStr();
               next if($atomstr eq "HETATM");
    
                my $d = $atom->Distance($r) ;
                my $nm = $atom->GetName() . " -> " . $r->GetName();
                my $info = {};
                $info->{NAME} = $nm;
                $info->{SCORE} = $d;
                push @resultsall, $info ;
         }
    
    }
    
    my @resultssorted = sort { $a->{SCORE} <=> $b->{SCORE} } @resultsall ;
    my $CNT = 0 ; 
    foreach my $r (@resultssorted){
        $CNT++;
        my $nm = $r->{NAME};
        my $score = $r->{SCORE};
        print STDERR "NAME = $nm score = $score \n";
        last if($CNT eq 10);
    }
    return @resultssorted ;
}


sub util_Ann2Simple{
    my ($infile,$outfile) = @_; 
    my $ofh ;
    $ofh = util_write($outfile) if(defined $outfile);
    my $ifh = util_read($infile);
    my @retlist  ; 
    while(<$ifh>){
         next if(/^\s*$/);
         chomp ;
         if(/^POINTS/){
             s/POINTS//g;
            s/\// /g;
            my @l = split " ", $_ ; 
            while(@l){
                my $a = shift @l ;
                my $b = shift @l ;
                #$a =~ s/\s*//g;
                #$b =~ s/\s*//g;
                print  $ofh "$b$a " if(defined $ofh);
                push @retlist, $a ;
            }
            print $ofh "\n" if(defined $ofh);
         }
    }
    return @retlist ;
}

sub util_GetPotDiffForAtoms{
     my ($pdb1,$pqr1,$pots1,$a,$b) = @_ ; 
     my $pota = util_GetPotForAtom($a,$pqr1,$pots1) ;
     my $potb = util_GetPotForAtom($b,$pqr1,$pots1) ;
     my $diff = $pota - $potb ;
     return $diff;
}

sub util_GetPotDiffForResidues{
    my ($pdb1,$pqr1,$pots1,$res1,$res2,$what) = @_ ; 
    my $a = $pdb1->GetAtomFromResidueAndType($res1->GetResNum(),$what);
    my $b = $pdb1->GetAtomFromResidueAndType($res2->GetResNum(),$what);
    my $diff = util_GetPotDiffForAtoms($pdb1,$pqr1,$pots1,$a,$b);
}

sub util_ProcessSingleLine{
     my ($pdb1,$pqr1,$pots1,$line,$isnew) = @_ ; 
     my $atomlist = $pdb1->ParseResultLine($line,$isnew);
     return util_ProcessSingleLineAtomlist($pdb1,$pqr1,$pots1,$atomlist);
}

sub util_ProcessSingleLineAtomlist{
     my ($pdb1,$pqr1,$pots1,$atomlist) = @_ ; 
     my @names ; 
     foreach my $a (@{$atomlist}){
        my $pot = util_GetPotForAtom($a,$pqr1,$pots1) ;
        push @names, $a->GetName();
     }
     my $name = join ",",@names ; ;
     my @dist = @{$pdb1->DistanceInGivenSetOfAtoms($atomlist)};
     my @pots = @{$pdb1->PotInGivenSetOfAtoms($atomlist,$pqr1,$pots1)};
     return (\@dist,\@pots,$name);
}

sub util_FindRmsd{
    my ($pdb1,$pdb2) = @_ ;
    
    my @res = $pdb1->GetResidues();
    my $N = @res;
    my $cnt = 0 ;
    my $sum = 0 ;
    my $cntmatch = 0 ; 
    foreach my $res (@res){
        next if($res->GetAtomStr() ne "ATOM");
        my $resnumA = $res->GetResNum() ;
        my $resnumB = $res->GetResNum() ;
        my $CAatom1 = $pdb1->GetAtomFromResidueAndType($resnumA,"CA");
        my $CAatom2 = $pdb2->GetAtomFromResidueAndType($resnumB,"CA");
        my $d = util_format_float($pdb1->DistanceAtoms($CAatom2,$CAatom1),1);
        $cnt++;
        $sum = $sum + $d * $d ; 
    }
    
    my $rmsd = util_format_float(sqrt($sum/$cnt),3) ; 
    print  " $rmsd $cnt\n";
    return $rmsd ;
}

sub util_GetMeanSD{
    my ($container) =@_ ;
    my $mean = util_format_float(Math::NumberCruncher::Mean($container),0) or warn "Mean not found" ;
    my $sd = util_format_float(Math::NumberCruncher::StandardDeviation($container),0) or warn "sd not found" ;
    return ($mean,$sd);
}

sub util_FindRmsdAllAtoms{
    my ($pdb1,$pdb2) = @_ ;
    
    my @atoms = $pdb1->GetAtoms();
    my $N = @atoms;
    my $cnt = 0 ;
    my $sum = 0 ;
    my $cntmatch = 0 ; 
    foreach my $atom1 (@atoms){
        my $type = $atom1->GetType();
        my $resnum = $atom1->GetResNum() ;
        my $atom2 = $pdb2->GetAtomFromResidueAndType($resnum,$type);
        next if (!defined $atom2);
        my $d = util_format_float($pdb1->distanceAtoms($atom1,$atom2),1);
        $cnt++;
        $sum = $sum + $d * $d ; 
    }
    
    my $rmsd = util_format_float(sqrt($sum/$cnt),3) ; 
    print  " $rmsd $cnt\n";
    return $rmsd ;
}
sub util_sortsingleString{
   my $s = shift ;
   my @sl = split "", $s ;
   my @XX = sort @sl ; 
   
   my $rev = join "", @XX ;
   return $rev ; 
}

sub util_AreResiduesContinuous{
    my ($pdb1,$res1,$res2,$eitherispro) = @_ ;
    my $CA1 = $pdb1->GetAtomFromResidueAndType($res1->GetResNum(),"CA");
    my $CA2 = $pdb1->GetAtomFromResidueAndType($res2->GetResNum(),"CA");
    my $d1 = util_format_float($pdb1->DistanceAtoms($CA1,$CA2),1);
    my $diff1 = abs ($d1 - 3.8);
    if($eitherispro){
        my $XXX = abs ($d1 - 3);
        $diff1 = $XXX if($XXX < $diff1);
    }
    return (0,$d1,$diff1) if($diff1 > 0.3);
    return (1,$d1,$diff1) ;
}

sub util_IsPro{
    my ($res) = @_ ;
    my $name = $res->GetName();
    return 1 if($name eq "PRO");
    return 0 ; 
}

sub util_IsResisdueThisType{
    my ($res,$type) = @_ ;
    my $name = $res->GetName();
    return 1 if($name eq "$type");
    return 0 ; 
}
sub util_IsMSE{
    my ($res) = @_ ;
    my $name = $res->GetName();
    return 1 if($name eq "MSE");
    return 0 ; 
}

## Spaces are removed...
sub util_readfasta{
    my ($infile) =@_ ; 
    my $str = "";
    my $ifh = util_read($infile);
    my $firstline ; 
    while(<$ifh>){
         if(/^\s*>/){
             $firstline = $_ ; 
             next ; 
         }
         next if(/^\s*$/);
         chomp ;
         $str = $str . $_ ; 
    }
    $str =~ s/\s*//g;
	close($ifh);
    return ($str,$firstline) ; 
}

sub util_ExtractSliceFromFasta{
    my ($ofh,$infile,$start,$end) = @_ ; 
    my ($str,$firstline) = util_readfasta($infile);

    chomp $firstline ;
    my $time = int(time() *rand());
    print "$time";
    $time = "$start.$end";
    
    print $ofh ">$time.$firstline  \n";
    
    my @l = split "", $str ; 
    my $N = @l ;
    print "$N sll \n";
    foreach my $i ($start..$end){
        $i = $i - 1 ; 
        print $ofh "$l[$i]";
    }
    print $ofh "\n";
}


## THIS NEEDS TO BE FIXED FOR FRAGMENTCOMPARE.PL ###
sub util_ReadAnnotateFile{
    my ($infile,$DIFF) = @_ ; 
    $DIFF = 50 if(!defined $DIFF);
    my $ifh = util_read($infile);
    my $info = {};
    while(<$ifh>){
         chomp; 
         next if(/^\s*$/);
         if(/^\s*(Repeat|Region)/){
             my (@l) = split ; 
             my $start = $l[1];
             my $end = $l[3];
             my $diff = abs($start - $end);
             next if($diff > $DIFF);
    
             my $anno = $_ ;
             print "$start $end $anno \n";
             $info->{$start} = $anno ;
             $info->{$end} = $anno ;
         }
         if(/^\s*Modified/){
             my (@l) = split ; 
             my $start = $l[2];
             my $anno = $_ ;
             $info->{$start} = $anno ;
         }
             
    }
    return $info ;
}

sub util_ReadAnnotateFileFRAGAL{
    my ($infile,$DIFF) = @_ ; 
    $DIFF = 50 if(!defined $DIFF);
    my $ifh = util_read($infile);
    my $inforepeatregion = {};
    my $infomodified = {};
    my $CNT = 0; 
    while(<$ifh>){
         chomp; 
         next if(/^\s*$/);
         $CNT++;
         if(/^\s*(Repeat|Region)/){
             my (@l) = split ; 
             my $start = $l[1];
             my $end = $l[3];
             my $diff = abs($start - $end);
             if($diff > $DIFF){
                 print "DIFF $diff is freater than $DIFF\n";
                next ;
             }
    
             my $anno = $_ ;
             $inforepeatregion->{$CNT}->{START} = $start;
             $inforepeatregion->{$CNT}->{END} = $end;
         }
         if(/^\s*Modified/){
             my (@l) = split ; 
             my $start = $l[2];
             my $anno = $_ ;
             $infomodified->{$CNT}->{START} = $start ;
         }
             
    }
    return ($inforepeatregion,$infomodified) ;
}


sub util_SortTwoStrings{
    my ($a,$b) = @_ ; 
    if($a lt $b){
        return ($a, $b);
    }
    else{
        return ($b, $a);
    }
    
}

sub util_ParseAAGroups{
    my ($in) = @_ ; 
    my $ifh = util_read($in);
    my $grp = 0 ; 
    my $info ={};
    my $map3to1 ={};
    my $map1to3 ={};
    my @grps ;
    while(<$ifh>){
        next if(/^\s*$/);
        #print ; 
        my (@l) = split ; 
        $grp++;
        my $NM ; 
        while(@l){
            my $single = shift @l ; 
           if(!defined $NM){
               $NM = $single ;
               push @grps, $single ;
           }
           #$info->{$single} = $grp ; 
           $info->{$single} = $NM ; 
           #print "$single \n";
            shift @l ; 
            shift @l ; 
            my $three = shift @l ; 

           $three =~ s/\(//;
           $three =~ s/\)//;
           $three = uc($three);
           $map1to3->{$single} = $three ;
           $map3to1->{$three} = $single ;
        }
    }
    return ($info,\@grps,$map3to1,$map1to3);
}

sub util_PairwiseDiff{
    my ($l1,$l2) = @_ ; 
    my @l1 = @{$l1};
    my @l2 = @{$l2};
    my $N1 = @l1 ;
    my $N2 = @l2 ;
    die "$N1 $N2 not equal" if($N1 ne $N2);
    my $diff = 0 ;
    while(@l1){
       my $a = shift @l1 ;
       my $b = shift @l2 ;
       $diff = $diff + abs ($a-$b);
    }
    return ($diff);
}

sub util_parseHETPDB{
    my ($infile) = @_ ;
    my $NOHET = {};
    my $YESHET = {};
    my $HET2PDB = {};
    my $HET2PDBSIZE = {};
    my $ifh = util_read($infile);
    while(<$ifh>){
         if(/^\s*NOHET/){
              my ($nm,@v) = split ; 
              $NOHET = util_make_table(\@v);
         }
         elsif(/^\s*YESHET/){
              my ($nm,@v) = split ; 
              foreach my $v (@v){
                     my ($protein,$num) = split ":", $v ;
                   $YESHET->{$protein} = $num ;
              }
        }
        else{
              my ($nm,$size,@v) = split ; 
              $HET2PDB->{$nm} = (\@v);
              $HET2PDBSIZE->{$nm} = $size ;  
        }
    }
    return ($NOHET,$YESHET,$HET2PDB,$HET2PDBSIZE);
}

sub util_parseHolo2Apo{
    my ($infile) = @_ ;
    my $info = {};
    my $ifh = util_read($infile);
    while(<$ifh>){
        my ($a,$b,$num) = split ;
        $info->{$a} = $b ; 
    }
    return $info ;
}

sub util_ParsePremonIn{
    my ($infile,$errfh) = @_ ;
    my $ifh = util_read($infile);
    my $info = {};
    my ($statusStr,$statusD,$statusPD,$statusREALSTRING,$resultstr);
    while(<$ifh>){
        print if($verbose);
        if(/^\s*STR /){
            $statusStr = 1 ;
            my ($junk,$str) = split ;
            $info->{STR} = [] if(!defined  $info->{STR});
            push @{$info->{STR}}, $str ;
        }
        if(/^\s*D /){
            $statusD = 1 ;
            my ($junk,@v) = split ;
            $info->{D} = \@v 
        }
        if(/^\s*PD /){
            $statusPD = 1 ;
            my ($junk,@v) = split ;
            $info->{PD} = \@v ;
            print "@v ......\n";
        }
        if(/^\s*REALSTRING /){
            $statusREALSTRING = 1 ;
            my ($junk,@v) = split ;
            $info->{REALSTRING} = \@v 
        }
        if(/^\s*RESULTSTYLE /){
            $resultstr = 1 ;
            my ($junk,@v) = split ;
            $info->{RESULTSTYLE} = \@v 
        }
    }
    my $status = 1 ; 
    if(!defined ($statusStr && $statusPD && $statusD && $statusREALSTRING && $resultstr)){
        print $errfh "Could not read $infile\n" ; 
        die ;
        $status = 0 ;
    }
    close($ifh);
    return ($status,$info);
}

sub util_ReadPremonOut{
   my ($premon,$size) = @_ ;
   my $info = {};
   my $ifh = util_read($premon);
   while(<$ifh>){
          next if(/^\s*$/);
          next if(/^\s*#/);
          my @l  = split ; 
          my $nm = shift @l ;
          my $len = length($nm);
          die "len = $len $nm " if($len ne $size);
          $info->{$nm} = \@l ;
   }
   close($ifh);
   return $info ;
}

sub util_ReadPremonOutAndGiveFasta{
   my ($premon,$size) = @_ ;
   my $info = {};
   my $ifh = util_read($premon);
   while(<$ifh>){
          next if(/^\s*$/);
          next if(/^\s*#/);
          my @l  = split ; 
          my $nm = shift @l ;
          my (@aacodes) = split "", $nm ;
          my $len = length($nm);
          die "len = $len $nm " if($len ne $size);
          foreach my $match (@l){
                $match =~ s/\./ /g;
                my (@numbers) = split " ", $match ;
                foreach my $idx (0..$size-1){
                    my $aa = $aacodes[$idx];
                    my $num = $numbers[$idx];
                    $info->{$num} = $aa ;
                }
                
          }  
   }

   my $fasta = "";
   foreach my $i (sort {$a <=> $b} keys %{$info}){
        my $s = $info->{$i};
        $fasta = $fasta . $s ;
   }
   close($ifh);

   return $fasta ;
}

## onus is onto caller to delete seq file
sub util_GetFasta{
    my ($pdb1,$name) = @_ ; 
    my $pdb = new PDB();
    $pdb->ReadPDB($pdb1);

    my $outfile = util_getTmpFile();
    my $ofh = util_write($outfile);
    my $seq = $pdb->WriteFasta("$name",$ofh);
    print STDERR "Wrote fasta in $outfile\n" if($verbose);

    close($ofh);
    return ($seq,$outfile) ;
}

sub util_ScorePD{
    my ($l1,$l2) = @_ ;
    my @l1 = @{$l1};
    my @l2 = @{$l2};
    my $N = @l1 ;
    die if($N ne @l2);
    my @listofdiffs ;
    foreach my $i (0..$N-1){
        my $p1 = $l1[$i];
        my $p2 = $l2[$i];
        my ($realdiff,$thr);
        my $absdiff = util_format_float(abs($p1 - $p2),3);
    
        my $ABSREF = abs($p1);
        my $ABSQUERY = abs($p2);
        $thr = 150 ;
            
        if($ABSREF > 400){
                 $thr = 170  ;
        }
        elsif($ABSREF > 300){
                 $thr = 150  ;
        }
        elsif ($ABSREF > 200){
                 $thr = 125  ;
        }
        #$thr = 75   if($ABSREF > 100 && $ABSREF < 131 && $ABSQUERY < $ABSREF);
        $thr = 75   if($ABSREF > 100 || $ABSQUERY > 100);
    
    
        if(abs($p1) < 99.9 && abs($p2) < 99.9){
            $thr = "bothbelow100";
            $realdiff = 0 ; 
         }
         else{
              $realdiff = $absdiff < $thr ? 0 : $absdiff;
         }
         push @listofdiffs, $realdiff ;

    }

    my @resultssorted = sort { $a <=> $b } @listofdiffs ; 
    my $indexformax = @resultssorted ;
    my $maxdiff = int($resultssorted[$indexformax - 1]);

    return ($maxdiff);
}

sub util_ScoreDistance{
    my ($l1,$l2,$MAXALLOWEDDISTDEV) = @_ ;
    my @l1 = @{$l1};
    my @l2 = @{$l2};
    my $N = @l1 ;
    die if($N ne @l2);
    my $sumsq =  0 ;
    my $maxdiff  = 0 ;
    foreach my $i (0..$N-1){
        my $a = $l1[$i];
        my $b = $l2[$i];

        my $diff = $a - $b;
        my $diffsq = $diff * $diff ;

        my $absdiff = util_format_float(abs($a - $b),1);
        $maxdiff = $absdiff if($absdiff > $maxdiff);

        $sumsq = $sumsq + $diffsq ;
    }
    my $rmsd = util_format_float(sqrt($sumsq/$N),1);
    if($maxdiff  > $MAXALLOWEDDISTDEV){
        $rmsd = $rmsd + $maxdiff ;
    }
    return ($rmsd,$maxdiff);
}

sub parseFPocketFile{
    my ($protein,$pdb1,$addtopdb,$infile) = @_ ;
    my $info = {};
    my $infoperatom = {};
    my $ifh = util_read($infile);
    my $volume  = 0 ;
    while(<$ifh>){
         next if(/^\s*$/);
         chomp ;

         if(/Real volume/){
             ($volume) = (/Real volume \(approximation\)\s*:(.*)/);
         }
    
        if(/^ATOM/){
           my ($atomstr , $serialnum , $atomnm , $alt_loc , $resname , $chainId , $resnum , $codeforinsertion , $x , $y , $z ) = util_ReadLine($_);

           $atomnm =~ s/\s*//g;
           my $I = $resname . "/" . $resnum ;
           my $ATOM = $I . "/" . $atomnm ;
           $info->{$resnum} = $I ; 
           $pdb1->AddPocket($resnum) if($addtopdb);
           $ATOM =~ s/\s*//g;
           $infoperatom->{$ATOM} = 1  ; 
           #print "$ATOM \n";
        }
    }
    return ($info, $infoperatom,$volume);
}

sub parseCastp{
    my ($infile,$subtract) = @_ ;
    my $info = {};
    my $infoperatom = {};
    my $ifh = util_read($infile);
    my $volume  = 0 ;
    while(<$ifh>){
         next if(/^\s*$/);
         chomp ;

        if(/^ATOM/){
           my ($atomstr , $serialnum , $atomnm , $alt_loc , $resname , $chainId , $resnum , $codeforinsertion , $x , $y , $z ) = util_ReadLine($_);
           my @l = split ;
           my $N = @l ;
           my $POCKNUMER = $l[$N-2];

           $atomnm =~ s/\s*//g;
           my $I = $resname . "/" . $resnum ;
           my $ATOM = $I . "/" . $atomnm ;
           if(!defined $info->{$POCKNUMER}){
                $info->{$POCKNUMER} = {};
           }
           $info->{$POCKNUMER}->{$resnum} = $I ; 

           #$ATOM =~ s/\s*//g;
           #$infoperatom->{$ATOM} = 1  ; 
           #print "$ATOM \n";
        }
    }

    my $biggestnumber ;
    foreach my $POCKNUMER ( sort { $a <=> $b} keys %{$info}){
        my @res = (keys %{$info->{$POCKNUMER}});
        my $N = @res ;
        #print "$POCKNUMER has $N res - castp \n";
        $biggestnumber = $POCKNUMER ;

    }
    my $infobiggest = $info->{$biggestnumber - $subtract};


    return ($info,$infobiggest);
}

sub util_AlignAndMatchRemaining{

   my ($p1,$p2,$infile,$maxdist) = @_ ;
   my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();

   
   my $file1 = "$PDBDIR/$p1.pdb";
   my $file2 = "$PDBDIR/$p2.pdb";
   my $pdb1 = new PDB();
   $pdb1->ReadPDB($file1);
   my $pdb2 = new PDB();
   $pdb2->ReadPDB($file2);

   my ($atoms1,$atoms2) = pymolin_getResultsLine($infile,$pdb1,$pdb2);
   
   return util_AlignAndMatchRemainingAtoms($pdb1,$pdb2,$atoms1,$atoms2,$maxdist);
       
   
}

sub util_AlignAndMatchRemainingAtoms{

   my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
   
   my ($pdb1,$pdb2,$atoms1,$atoms2,$maxdist) = @_ ;

   #print "Aligning geom_Align3PointsToXYPlane \n";
   my ($done1,$remainingAtoms1) = MyGeom::geom_Align3PointsToXYPlane($pdb1,$atoms1,$verbose);
   my ($done2,$remainingAtoms2) = MyGeom::geom_Align3PointsToXYPlane($pdb2,$atoms2,$verbose);
   
   
   my $i = shift(@{$remainingAtoms1}) or die "Need an extra residue";
   my $JJ = shift(@{$remainingAtoms2}) or die "Need an extra residue";
   my $NM2 = $JJ->GetName();
   
       my @tmp1 = (@{$done1}, $i);
       $i->Print() if($verbose);
       my @atomlist ;
       push @atomlist, $i ;
       my ($results,$combined) = $pdb2->GetNeighbourHoodAtom(\@atomlist,$maxdist);
       print STDERR  "Atoms close to this one\n" if($verbose);
       my $sort ;
       foreach my $j (@{$combined}){
           my $NM = $j->GetName();
           if($NM eq $NM2){
           my $d = util_format_float($i->Distance($j),1) ;
              print "$d $NM ;;;;;;;;;;\n" if($verbose);
              return $d; ;
           }
       }
   
       return 1000 ;
       
   
}

sub util_ParseDSSP{
   my ($protein,$pdb1,$dsspinfile,$what,$finaloutfile,$writeIndivual,$force) = @_ ;
   print "Parsing DSSP for $protein for $what\n";

   my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP) = util_SetEnvVars();



   ## make the dssp file is not there 
   if($force ||  ! -e $dsspinfile){
       print "\tWriting DSSP $dsspinfile\n";
       system("mkdssp -i $PDBDIR/$protein.pdb -o $DSSP/$protein.dssp");
   }
   
   die "Error: Need dsspinfile" if(! -e $dsspinfile);
   my $ifh = util_read($dsspinfile);
   while(<$ifh>){
        last if(/#  RESIDUE AA STRUCTURE BP1 BP2/);
   }
   
   my @arr = ();
   foreach my $i (0...10000){
       push @arr, 0 ;
   }
   
   my $final = 0 ;
   while(<$ifh>){
        next if(/^\s*$/);
        chomp ;
        my (@l) = split ; 
        my $n = $l[1];
        $final= $n;
        my $h = $l[4];
        if($what eq "HELIX"){
           #if($h eq "H" || $h eq "G" || $h eq "I"){
           if($h eq "H" ){
                    $arr[$n] = 1 ;
           }
        }
        else{
            #if($h eq "B" || $h eq "E"){
            if( $h eq "E"){
                $arr[$n] = 1 ;
            }
         }
   }
   $final++;
   #print "Setting 2 value to $final\n";
   $arr[$final]= 2;
   
   my $idx = 0 ;
   my $prevend = 0 ;
   my $prevstart = 0 ;
   my $start = 0 ;
   my $helixnumber = 1 ;

   ## Figure out the HTH - if diff is less than hardcoded now
   #my $ofhhth = util_open_or_append("HTH");
   while($idx < 10000){
      ($idx,$start) = __GetNextStretch(\@arr,$idx,$helixnumber,$protein,$pdb1,$finaloutfile,$what,$writeIndivual);
      last if($idx eq -1);
      my $diff = $start - $prevend ;
      if($diff < 5 && $helixnumber ne 1){
          my $prevhnumber = $helixnumber -1 ;
          #my $QQQ = "$protein.${what}$prevhnumber";
          #my $PPP = "$protein.${what}$helixnumber";
          my $QQQ = "${what}$prevhnumber";
          my $PPP = "${what}$helixnumber";
          #print $ofhhth "$QQQ $PPP $protein $diff $prevstart $idx   \n";
      }
   
      $prevend = $idx; 
      $prevstart = $start; 
      $idx++;
      $helixnumber++;
   
   }
   

   $helixnumber--;

   if($helixnumber){
       print "\tSplit in $helixnumber of $what\n";
   }
   else{
        my $donefh = util_open_or_append($finaloutfile);
		print $donefh "$protein NULL $what\n";
	}
   
   close($ifh);
   return $finaloutfile ;
}
   
sub __GetNextStretch{
       my ($l,$idx,$hnumber,$protein,$pdb1,$finaloutfile,$what,$writeIndivual) = @_;
    my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR) = util_SetEnvVars();
       my @l = @{$l};
       my @ret ; 
       my $start;
       my $end;
       foreach my $i ($idx...10000){
           my $x = $l[$i];
           return -1 if($x eq 2);
           if($x eq 1){
               $start = $i;
               last ;
           }
       }
       return -1 if(!defined $start) ;
   
       foreach my $i ($start...10000){
           my $x = $l[$i];
           if($x eq 0 || $x eq 2){
               $end = $i -1 ;
               last ;
           }
       }
           die "start = $start " if(!defined $end);
           my ($val,$tableofres,@listofres) = util_GetSequentialSetResidues($pdb1,$start,$end) ;
		   if(!$val){
             my $donefh = util_open_or_append($finaloutfile);
             print $donefh "$protein ${what}$hnumber $start $end broken\n ";
             return ($end,$start) ;
		   }

        if($writeIndivual){
              my $PPP = "$HELIXDIR/$protein.${what}$hnumber";
              my $ofhPDBHELIX = util_open_or_append($PPP. "${what}LIST");
              print $ofhPDBHELIX "${what}$hnumber\n";
              my $OOO = "$PPP.pdb";
              my $ofhhelix= util_write($OOO);
              print "Writing hnumber $hnumber to file $OOO : residues $start to $end \n" ;
              my $pdb = "$PDBDIR/$protein.pdb";
              $pdb1->ReadPDBAndWriteSubset($pdb,$tableofres,$ofhhelix);
              close($ofhhelix);
        }


        my $donefh = util_open_or_append($finaloutfile);


        if($what eq "HELIX"){
           $pdb1->AddHelix($hnumber,$start,$end);
           print $donefh "$protein ${what}$hnumber $start $end ";
           util_helixWheel ($protein,$pdb1,\@listofres,$donefh,$what,$hnumber,$writeIndivual);
        }
        else{
           $pdb1->AddBeta($hnumber,$start,$end);
           print $donefh "$protein ${what}$hnumber $start $end \n";
        }

       return ($end,$start) ;
}


sub util_helixWheel{
   my ($protein,$pdb1,$listofres,$donefh,$what,$hnumber,$writeIndivual,$justseqeunce) = @_ ; 
   die if(!defined $listofres && !defined $justseqeunce);

   my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP) = util_SetEnvVars();

    
    my ($tableATOM,$HYDROVAL,$colortable,$value,$chargedtable) = Config_Helix();

    my $hydroval = {};
    foreach my $k (keys %{$HYDROVAL}){
        my $single = $pdb1->GetSingleLetter($k);
        $hydroval->{$single} = $HYDROVAL->{$k};
    }
    
    
    
   my ($outfile,$ofh,$ofhpymol,$ofhcommands) ;
   if($writeIndivual){
       $outfile = "$protein.wheel.tikz.tex";
       $ofh = util_write("$outfile");
       $ofhpymol= util_write("$protein.pymol.helix.in");
       $ofhcommands = util_open_or_append("$protein.commands");
       print $ofh "% for protein $protein \n";
   }
        
   
  my @singles ;
  my @numbers ;
  my $N ;
  if(! defined $justseqeunce){
    my @res = @{$listofres};
    $N = @res;
    my $cnt = 0 ;
      foreach my $r (@res){
            #die "There are HET Atoms, dd you not divide into helices?" if($r->GetAtomStr() eq "HETATM");
            #die if($r->GetName() eq "HOH");
			next if($r->GetName() eq "HOH");
			next if($r->GetName() eq "HETATM");
			my $NM = $r->GetName()  ;
            my $single =  $r->PrintSingleLetter($pdb1);
            push @singles, $single ;
            my $num = $r->GetResNum();
			my $x  = $r->PrintSingleLetter($pdb1);
			#print "$NM $num $x ikkkkkkkkkkkkkk\n";
            push @numbers, $num ;
       }
   }
   else{
       my $ifh = util_read($justseqeunce);
       my $sequence = "";
       while(<$ifh>){
                 next if(/^\s*$/);
                 next if(/^\s*>/);
                 chomp ;
                 $sequence = $sequence . $_ ;
       }
       $sequence =~ s/\s*//g ;
       @singles = split "", $sequence ;
       $N = @singles;
	   my $OFFSET = 426;
	   print "HACKOFFSET = $OFFSET\n";
       foreach my $i (1..$N){
	             $i = $i + $OFFSET -1 ;
                push @numbers, $i ; 
       }
       close($ifh);
   }
    
    
        #my @lol ;
        #$lol[0] = []; $lol[1] = []; $lol[2] = []; $lol[3] = [];
        my $done = {};
        my $initval = -270 ; 
        my $rad = 5 ; 
        my $loopcnt = 0 ;
    
        my @centre = qw (0 0 0);
        my $finalveccharged; 
        my $finalvechydro; 
    
        my $TOTALCHARGED = 0 ;
        my $TOTALPOS = 0 ;
    
    
        my $bins = {};
        my $TOTALNUM = 0 ; 
        my $circularseq = {};
        my $circularval = {};

		#my $NNNNNNN = @singles ;
		#print "$NNNNNNN lllllllllllllll @singles \n";

        while(@singles){
            my $single =  shift @singles ;
			next if($single =~ /^\s*$/);
            $TOTALNUM++;
            my $num = shift @numbers ;
    
    
            my $idxcnt = $loopcnt % 4 ;
            if(!defined $bins->{$idxcnt}){
                $bins->{$idxcnt} = [];
            }
    
    
            my $x = $pdb1->{SINGLE2THREE}->{$single} ;
			if(!defined $x){
				die "XXX $protein @singles $single $x something wrong in util_helixWheel \n";
				die ;
			}
            my $what = $tableATOM->{$x};
			if(!defined $what){
				die "$protein $single $x something wrong in util_helixWheel \n";
				die ;
			}
            #print "$x/$num/$what $idxcnt $x\n";
            my $SSSS = $single . $num;
    
            push @{$bins->{$idxcnt}}, "$x/$num/$what" ;
    
            my $val = $initval - $loopcnt * 100 ; 
            
            my $cval = (deg2rad($val ));
            while( exists $circularval->{$cval}){
                $cval = $cval + 0.01;
            }
            $circularseq->{$SSSS} = $cval ; 
            $circularval->{$cval} = 1 ;
    
            $x = util_format_float($rad * cos(deg2rad($val)),1) + 0 ;
            my $y = util_format_float($rad * sin(deg2rad($val)),1) + 0 ;
            my $str = "$x.$y";
    
            my @thispoint ; 
            push @thispoint, $x ;
            push @thispoint, $y ;
            push @thispoint, 0 ;
    
                my $VEC = MyGeom::MakeVectorFrom2Points(\@centre,\@thispoint) ;
                my $tmpvec = $VEC->norm * $hydroval->{$single} ;
                if(!defined $finalvechydro){
                    $finalvechydro = $tmpvec ;
                }
                else{
                    $finalvechydro = $finalvechydro + $tmpvec ;
                }
    
            if(exists $chargedtable->{$single}){
                $TOTALCHARGED++;
                my $factor = $chargedtable->{$single} ;
                if($factor eq 1){
                    $TOTALPOS++;
                }
                #my $tmpvec = $VEC->norm * $hydroval->{$single} *$factor;
                my $tmpvec = $VEC->norm * $factor;
                if(!defined $finalveccharged){
                    $finalveccharged = $tmpvec ;
                }
                else{
                    $finalveccharged = $finalveccharged + $tmpvec ;
                }
            }
    
            while(exists $done->{$str}){
                #print "$str existed\n";
                $y = $y + 1 ;    
                $str = "$x.$y";
            }
            $done->{$str} = 1 ;
            my $color = $colortable->{$single};
            my $VAL = $value->{$single};
            

            if($writeIndivual){
               print $ofh "\\node[fill=$color!$VAL,draw=blue,very thick] (FinalNode) at ($x, $y) {$single$num} ; \n";
               print $ofhpymol "select A, /PDBB//A/$num/CA\n";
               print $ofhpymol "# $single\n";
               print $ofhpymol "color $color, A\n";
            }
    
    
            $loopcnt++;
            #push @{$lol[$idxcnt]}, $r ;
        }

    my $Lcharged = 0  ;
    my $cx = 0 ;
    my $cy = 0 ;
    if(defined $finalveccharged){
        $Lcharged = util_format_float($finalveccharged->length,1) ;
        ($cx) = $finalveccharged->x * 1;
        ($cy) = $finalveccharged->y * 1;
    }
    my $Lhydro = util_format_float($finalvechydro->length,1) ;
    
    my ($hx) = $finalvechydro->x * 1;
    my ($hy) = $finalvechydro->y * 1;
    
    my $PERCENTPOS = -1 ;
    if($TOTALCHARGED){
        $PERCENTPOS = util_format_float($TOTALPOS/$TOTALCHARGED,1);
    }
    
    
    
   my $writecircular = 0 ;
		if($writecircular){
            my $ofhcircular = util_open_or_append("CIRCULARFASTA");
            my $fastacircular = "$protein.${what}$hnumber";
            print $ofhcircular ">$fastacircular\n";
            foreach my $k (sort { $circularseq->{$a} <=> $circularseq->{$b}} keys %{$circularseq}){
                my $v = $circularseq->{$k};
                my ($firstletter) = ($k =~ /(.)/);
                print $ofhcircular "$firstletter";
        
            }
            print $ofhcircular "\n";
		}
    
    
      if($writeIndivual){
           #### FOR CLASP #########
           my ($REALPDB) = ($protein =~ /(....)/) ;
           $REALPDB = $REALPDB . "A";
           foreach my $k (sort  keys %{$bins}){
               my $l = $bins->{$k} ;
               my $CONCATSTR = "";
               my $cntconcat = 0 ;
               foreach my $i (@{$l}){
                   $cntconcat++;
                   $CONCATSTR = $CONCATSTR . " $i ";
               }
               if($cntconcat > 1){
                   print $ofhcommands "echo $CONCATSTR > ! IN \n";
                   my $CCC = "printPairwiseOneAfterAnother.pl";
                   print $ofhcommands "$CCC -out out -c \$CONFIGGRP -ra 222 -pr $REALPDB -in IN -helix $protein \n";
               }
           }
    

            print $ofh "\\node[fill=black!10,draw=yellow,very thick] (HPH) at ($hx, $hy) {\$\\mu_{h}\$} ; \n" ;
            #print $ofh "\\node[fill=white!100,draw=blue,very thick] (Ch) at ($cx, $cy) {Ch} ; \n" ;
            print $ofh "\\node[fill=blue!100,draw=blue] (centre) at (0, 0) {} ;\n" ;
        
            print $ofh "\\begin{scope}[every path/.style=line]\n" ;
            print $ofh "\\path          (centre) -- node [near end] {$Lhydro} (HPH);\n" ;
            #print $ofh "\\path          (centre) -- node [near end] {$Lcharged} (Ch);\n" ;
            print $ofh "\\end{scope}\n";
    
           close($ofh) ;
           system("unlink in ; ln -s $outfile in");
    } ## writeIndivual

   # this appends to the line already printed
   print $donefh " $N $Lhydro $PERCENTPOS $TOTALCHARGED \n";

	$outfile = "NA" if(!defined $outfile); ## can happen for writeIndivual=0
    print "TIKZ o/p written to $outfile Lcharged = $Lcharged and Lhydro = $Lhydro PERCENTPOS = $PERCENTPOS TOTALNUM = $TOTALNUM \n" if($verbose);

}


sub util_GetDisulphide{
    my ($protein,$pdb1,$plsprint) = @_ ;
    my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP) = util_SetEnvVars();
    my @res = $pdb1->GetResidues();
    my $N = @res;

    my $CNT = 0 ;
    while(@res){
        
        my $res1 = shift @res ;
        foreach my $res2(@res){
            next if($res1->GetAtomStr() eq "HETATM");
            next if($res2->GetAtomStr() eq "HETATM");
            next if($res1->GetName() eq "HOH");
            next if($res2->GetName() eq "HOH");
    
            my $bothcys = util_IsResisdueThisType($res1,"CYS") && util_IsResisdueThisType($res2,"CYS") ;
            next if(!$bothcys);
    
    
            my $num1 = $res1->GetResNum();
            my $num2 = $res2->GetResNum();
            my $SG1 = $pdb1->GetAtomFromResidueAndType($res1->GetResNum(),"SG");
            my $SG2 = $pdb1->GetAtomFromResidueAndType($res2->GetResNum(),"SG");
    
            my $distSG = util_format_float($pdb1->DistanceAtoms($SG1,$SG2),1);
            if($distSG < 2.2){
                   print "Definite! dist = $distSG between SG for  resnum=$num1 resnum=$num2 \n" if($plsprint);
                   $CNT++;
                   $pdb1->AddDisulphide($num1,$num2);
            }
            elsif($distSG < 4){
                   print "possible? dist = $distSG between SG for  resnum=$num1 resnum=$num2 \n" if($plsprint);
                   $CNT++;
                   $pdb1->AddDisulphide($num1,$num2);
            }

         }
    }


    return $CNT;

}

sub util_NeedleFiles{
    my ($outfile,$f1,$f2,$arg) = @_ ;
    my $end = " -endopen 1 -endweight 1 -endextend 0.2";
    my $execinit = ParseArgFile($arg);
    #print "INIT = $execinit \n";
    #$exec = "needle -gapopen 25  -gapex 0.5 -ou $outfile $file1 $file2  2> /dev/null ";
    my $exec = $execinit . "-ou $outfile $f1 $f2  2> /dev/null "; 
    system("touch $outfile");
    system($exec);

   my ($iden,$simi); 
   my $ifh = util_read($outfile);
   while(<$ifh>){
        next if(/^\s*$/);
        if(/Identity/){
            ($iden) = (/.*\((.*)\)/);
            $iden =~ s/\%//;
        }
        if(/Similarity/){
            ($simi) = (/.*\((.*)\)/);
            $simi =~ s/\%//;
        }
        
   }
    return ($iden,$simi);
}

sub util_NeedleSeq{
    my ($outfile,$seq1,$seq2,$arg,$pdb1name,$pdb2name) = @_ ;
        
    my $f1 = (util_getTmpFile());
    my $f2 = (util_getTmpFile());
    my $ofh1 = util_write($f1);
    my $ofh2 = util_write($f2);

    print $ofh1 ">$pdb1name";
    print $ofh1 "$seq1\n";

    print $ofh2 ">$pdb2name\n";
    print $ofh2 "$seq2\n";

    close($ofh1);
    close($ofh2);

    return util_NeedleFiles($outfile,$f1,$f2,$arg);
}

sub util_NeedlePDBNamesFromSeq{
    my ($outfile,$pdb1name,$pdb2name,$arg) = @_ ;
    my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP) = util_SetEnvVars();
    my $pdb1 = "$PDBDIR/$pdb1name.pdb";
    my $pdb2 = "$PDBDIR/$pdb2name.pdb";

    my ($seq1,$f1) = util_GetFasta($pdb1,$pdb1name);
    my ($seq2,$f2) = util_GetFasta($pdb2,$pdb2name);

    
    my $N1 = length($seq1);
    my $N2 = length($seq2);
    if(!$N1 || !$N2){
        print "empty sequence in one $N1 or $N2\n";
        return (-1,-1)  ;
    }



   my ($iden,$simi)=     util_NeedleFiles($outfile,$f1,$f2,$arg,$pdb1name,$pdb2name);
   unlink($f1);
   unlink($f2);
   return ($iden,$simi);
}



sub util_NeedlePDBNamesFromFASTADIR{
    my ($outfile,$pdb1name,$pdb2name,$arg) = @_ ;
    my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP) = util_SetEnvVars();
    my $f1 = "$FASTADIR/$pdb1name.ALL.1.fasta";
    my $f2 = "$FASTADIR/$pdb2name.ALL.1.fasta";

    if(! -e $f1 || ! -e $f2){
        die "$f1 or $f2 does not exist \n";
    }

    return util_NeedleFiles($outfile,$f1,$f2,$arg);
}

sub ParseArgFile{
    my ($arg) = @_ ;
    my $exec  = "needle " ;
    my $ifh = util_read($arg);
    while(<$ifh>){
         next if(/^\s*$/);
         my ($nm,$junk) = split ; 
         $exec = $exec . "-". $nm . " $junk " ;
    }
    return $exec ; 
}

sub utils_parseBlastOut{
    my ($infile,$cutoff) = @_ ;
    my $ifh = util_read($infile);
    my @l ; 
    my $push = 0 ; 
    while(<$ifh>){
         next if(/^\s*$/);
         chomp ;
         if(/Sequences producing significant alignments/){
            $push = 1 ; 
            next ;
         }
         if($push && !/(chromosome|hypothe|unnamed|uncharacterized)/i){
             push @l, $_ 
         }
    }
    my $v = $l[0] ;
    die "$infile problem" if(!defined $v);
    
    if(defined $cutoff){
         my @lll = split " ", $v ;
         my $N = @lll ;
         my $evalue = $lll[$N-1];
         if($evalue > $cutoff){
             print "Warning; $infile has evalue $evalue more than cutoff $cutoff\n";
         }
    }
    return $v ;
}

## removed - better way
sub util_Blastout_isrev{ die ; }


sub util_ConcatTwoStringsWithCommonBeginAndEnd{
    my ($A,$B,$N) = @_;

    foreach my $i (5..$N){
        my ($junk1,$Ax) = util_GetTerminalStrings($A,$i);
        my ($Bx,$junk2) = util_GetTerminalStrings($B,$i);
        if($Ax eq $Bx){
            $A =~ s/$Ax//;
            my $ret = $A . $B ;
            return $ret ;
        }
    }
    return -1 ;

}

sub util_GetTerminalStrings{
    my ($seq,$size) = @_ ;

	my $STR2MATCH = "";
	foreach my $i (1..$size){
		$STR2MATCH = $STR2MATCH . ".";
	}
    my ($begin) = ($seq =~/^($STR2MATCH)/);
    my ($end) = ($seq =~/($STR2MATCH)$/);
    return ($begin,$end) ;
}

sub util_ParseWebBlast{
   my ($infile) = @_ ;
   my $ifh = util_read($infile);
   my $info = {};
   my ($length,$percent);
   my $push = 0 ; 
   my $pushident = 0 ; 
   my @STRS;
   my $fulllength  ;
   my $identities  ;
   my $str ; 
   
   my $SSSS ;
   my @LLLs ;
   my @TRSnm ;
   my @SCORES ;
   my @EVALUES ;
   my @IDENTITY ;
   while(<$ifh>){
        next if(/^\s*$/);
        chomp ;
        if(/No hits found/){
            return ($fulllength,undef) ;
        }
        if(/Length/ && ! defined $fulllength){
           ($fulllength) = (/Length=(\d+)/);
        }
        elsif(/Identities/){
            ($identities,$length,$percent) = (/Identities\s*=\s*(\d+)\/(\d+)\s*\((\d+)\%\)/);
           #$str = "$identities / $length $percent";
            push @IDENTITY, $percent;
            #push @STRS, "$SSSS  $SCORES[@SCORES-1]";

           my $evalue = $EVALUES[@EVALUES-1] + 0  ;
            push @STRS, "\t$SSSS\t$evalue";
            push @TRSnm, $SSSS ;
           push @LLLs, $length ;
        }
        elsif(/Score =/){
            s/,//g;
            my @lll = split;
            push @SCORES, $lll[2];
            push @EVALUES, $lll[7];
        }
        if(/^>/){
            s/>//;
           $SSSS = $_ ;
        }
   }
   if(!defined $fulllength){
            return (undef) ;
    }
   return ($fulllength,\@TRSnm,\@STRS,\@LLLs,\@SCORES,\@IDENTITY,\@EVALUES);
}

sub util_extractSliceFromFasta{
   my ($str,$start,$end) = @_;
   return util_extractSliceFromFastaString($str,$start,$end);
}


sub util_extractSliceFromFastaString{
   my ($str,$start,$end) = @_;
   die if(!defined ($start && $end));
   print "util_extractSliceFromFasta $start $end\n" if($verbose);

   my @l = split "", $str ; 
   my $N = @l -1 ;
   my $diff = abs ($start - $end) + 1 ;
   print "$N+1 length, asked for $diff  \n" if($verbose);
   die "start > end $start > $end " if($start > $end);
   

   my $retstr = "";
   foreach my $idx ($start...$end){
       #print "$idx $start $end\n";
       $idx = $idx -1 ;
       die "$idx asked from list of length $N" if(! defined $l[$idx]);
       $retstr = $retstr . $l[$idx] ;
   }
   return $retstr ;

}




sub util_R_extractIndex{
    my ($INF,$var,$idx) = @_ ;
    my $ifh = util_read($INF);
    my @retl;
    while(<$ifh>){
         next if(/^\s*$/);
         my (@l) = split " " ,$_ ;
          push @retl, $l[$idx];
    }
    close($ifh);
    return @retl ;
}

sub util_get_maxinlist{
    my (@l) = @_ ;
	return util_get_max(\@l);
}

sub util_get_max{
    my ($l) = @_ ;
    my @l = @{$l};
    my $max = -100000;
    foreach my $i (@l){
        $max = $i if($i > $max) ;
    }
    return $max ;
}



sub util_getComplimentaryString{
   my ($STR) = @_ ;
   die if($STR =~ /u/i);
   $STR =~ tr/ACGTacgt/TGCAtgca/;
   my $rev = reverse $STR ;
   return $rev ;
}

sub util_ParseBlastPW{
    my ($infile) = @_ ;
    my $ifh = util_read($infile);
    while(<$ifh>){
        if(/Expect/){
            my ($val) = (/Expect = (.*),/);
            if(!defined $val){
            ($val) = (/Expect = (.*)/);
            }
            return $val ;
        }
    }
    close($ifh);
    return 10000000000 ; 
}

sub util_ReadCodonBiasFile{
   my ($codonfile) = @_ ;
   my $ifh = util_read($codonfile);
   
   my $sort = {};
   while(<$ifh>){
        next if(/^\s*$/);
        my ($nm,@l) = split ; 
       $sort->{$nm} = $_ ;
   }
   my $info = {};
   my $mapValues = {};
   foreach my $k (sort keys %{$sort}){
       $_ = $sort->{$k};
       my $at = 0 ;
       my $gc = 0 ;
       my $py = 0 ;
       my $pu = 0 ;
       next if(/^\s*$/);
       my ($nm,@l) = split ; 
       my $sum = 0 ;
   
       $info->{$nm} = {};
       $mapValues->{$nm} = {};
   
       my $start = 0 ;
       foreach my $i (@l){
            $i =~ s/=/ /;
           my ($id,$junk) = split " ",$i ;
           $mapValues->{$nm}->{$id} = $junk ;
   
   
           my $val = 10* $junk ;
           my $end = $val+$start ;
           #print "$id $junk $val $end \n";
           foreach my $idx ($start..$end){
               $info->{$nm}->{$idx} = $id ;
           }
           $start = $end;
       }
   }
   return ($info,$mapValues) ;
}

sub util_ExtractOneIdxFromFile{
    my ($infile,$idx,$varname,$fix) = @_ ;
   my @list ;
   my $ifh = util_read($infile);
   while(<$ifh>){
     next if(/^\s*$/);
     my (@l) = split ;
     my $N = @l -1 ;
     next if($idx > $N);

     my $NNN = $l[$idx];

     $NNN =~ s/\(//;
     $NNN =~  s/\)//;

     if(defined $fix){
         $NNN = $NNN/$fix;
     }
     push @list, $NNN ;
   }
   close($ifh);


   my ($mean,$sd) = util_GetMeanSD(\@list);
   my $join = join " ,", @list ;
   my $assign = " $varname <- c ($join) \n";
   return (\@list,$mean,$sd,$assign) ;
}




### gets the codon bias from ORF and nt fasta file
sub util_ProcessOneORF{
  my ($trs,$orfile,$fastafile,$orf,$info,$cutoff) =@_ ;
  my ($actuallprocessed,$ignoredduetounknown,$ignoredduetosize);
  $actuallprocessed = $ignoredduetounknown = $ignoredduetosize = 0 ;
  my $ifh = util_read($orfile);
  my @NT_RET ;
  my @AA_RET ;
  my $exactcodingseq = "";
  my $found = 0 ;
  while(<$ifh>){
      if(/$orf/){
         s/-//g;
         s/\[//g;
         s/\]//g;
         my @l = split ;
         my $start = $l[1];
         my $end = $l[2];

         my $isreverse = 0 ;
         if(/REVERSE/){
             $isreverse = 1 ;
         }

         my $SSS = "";
         while(<$ifh>){
             chomp;
             if(/^\s*>/){
                last ;
            }
            $SSS = $SSS . $_ ;
            
         }
         my $LLL = length($SSS);
         if(defined $cutoff && $LLL < $cutoff){
             #print "Ignoring $trs\n";
            $ignoredduetosize++;
            next ;
         }
         
         my ($str,$firstline) = util_readfasta($fastafile);
		 my $ignoreunknowninfasta = Config_NTFastaUnknown(); # "R|Y|K|M|S|W|B|D|H|V|N"
         if($str =~ /($ignoreunknowninfasta)/){
             $ignoredduetounknown++;
            next ;
         }
         $actuallprocessed++;


         if($isreverse){
             my $tmp = $start ;
            $start = $end ;
            $end = $tmp ;
         }


         my $retstr = util_extractSliceFromFasta($str,$start,$end);
         my $len = length($retstr);



         if($isreverse){
            $retstr = util_getComplimentaryString($retstr) ;
            print "Reverse - so changing start-end, and also complimenting\n" if($verbose);
         }

         my $rlen = length($retstr);
         die "Something wrong" if($LLL*3 ne $len);
         print "Extracting $start $end , got $LLL for aa and $len for nt. isreverse = $isreverse and $rlen = rlen\n" if($verbose);
         if($verbose){
             $found = 1 ;
             my $FFFFFF = util_open_or_append("$trs.$orf.FFFFFF.fasta");
             print $FFFFFF ">$trs\n" ;
             print $FFFFFF "$retstr\n"  ;
         }
         $exactcodingseq = $retstr ;


         @NT_RET = ($retstr =~ /(...)/g);
         @AA_RET = ($SSS =~ /(.)/g);


         util_GetCodonBiasFromNucleotide($info,$retstr,$SSS);
         last ;

      }
  }
  die "did not find $orf" if(!$found && $verbose);
  close($ifh);
  return (\@NT_RET,\@AA_RET,$exactcodingseq,$actuallprocessed,$ignoredduetounknown,$ignoredduetosize);
}

sub util_ProcessOneORFGivenExactCDNA{
         my ($trs,$fastafile,$info,$cutoff) =@_ ;
         my ($actuallprocessed,$ignoredduetounknown,$ignoredduetosize);
         $actuallprocessed = $ignoredduetounknown = $ignoredduetosize = 0 ;
         my @NT_RET ;
         my @AA_RET ;
         my $exactcodingseq = "";
         my ($str,$firstline) = util_readfasta($fastafile);

		 my $ignoreunknowninfasta = Config_NTFastaUnknown(); # "R|Y|K|M|S|W|B|D|H|V|N"
         if($str =~ /($ignoreunknowninfasta)/ || !($str =~ /^ATG/)){
             $ignoredduetounknown++;
            return (0);
         }
         $actuallprocessed++;



         my $retstr = $str ;
         my $len = length($retstr);


         my $SSS = util_ConvertNucleotide2Amino($retstr,1);

         my $LLL = length($SSS);
         die "Something wrong $trs $LLL $len \n$retstr \n $SSS\n " if($LLL*3 ne $len);
		 #print "$SSS \n";
         $exactcodingseq = $retstr ;


         @NT_RET = ($retstr =~ /(...)/g);
         @AA_RET = ($SSS =~ /(.)/g);


         util_GetCodonBiasFromNucleotide($info,$retstr,$SSS);
         return (\@NT_RET,\@AA_RET,$exactcodingseq,$actuallprocessed,$ignoredduetounknown,$ignoredduetosize);
}


sub util_ConvertNucleotide2Amino{
         my ($retstr,$replaceendcodonwitX) = @_ ;
         my @NT = ($retstr =~ /(...)/g);
         my $codontable = Config_getCodonTable();
		 my $SSS = "";
         foreach my $i (@NT){
                 my $a = $codontable->{$i} ;
                if(!defined $a){
					if($replaceendcodonwitX){
                        $a = "X";
					}
					else{
                         die "$i not defined, in $retstr";
					}
                }
                $SSS = $SSS . $a ;
         }
		 return $SSS ;
}

sub util_GetCodonBiasFromNucleotide{
         my ($info,$retstr,$SSS) = @_ ;
         my @NT = ($retstr =~ /(...)/g);
         my $codontable = Config_getCodonTable();
         if(!defined $SSS){
             $SSS = "";
             foreach my $i (@NT){
                 my $a = $codontable->{$i} ;
                if(!defined $a){
                    die "$i not defined";
                }
                $SSS = $SSS . $a ;
             }
         }

         my @AA = ($SSS =~ /(.)/g);
         my $N = @AA ;
         my $NTN = @NT ;
         die "$N =n an $NTN = NTN" if($N ne $NTN);
         while(@NT){
             my $a = shift @NT ;
             my $b = shift @AA ;
            if(!defined $codontable->{$a}){
				#next ;
                die "$a";
            }

            die "$a $b while this is in $codontable->{$a} " if($codontable->{$a} ne $b);
            #print "$a $b\n";
            if(!exists $info->{$b}){
                $info->{$b} = {};
            }
            if(!exists $info->{$b}->{$a}){
                $info->{$b}->{$a} = 0 ;
            }
            $info->{$b}->{$a} = $info->{$b}->{$a} + 1 ;
         }
         return ($SSS);
}

sub util_EmitCodonBiasfromInfo{
    my ($ofh,$info) = @_ ;
    foreach my $aa (sort keys %{$info}){
        my $tab = $info->{$aa} ;
        print $ofh "$aa ";
        foreach my $k (keys %{$tab}){
            my $v = $tab->{$k} ;
            print $ofh " $k=$v ";
        }
        print $ofh "\n";
    }
}



sub util_GetCentreOfMass{
   my ($pdb1,$isCAOnly) = @_;
   my @res = $pdb1->GetResidues();
   return util_GetCentreOfMassFromSet($pdb1,$isCAOnly,@res);
}

sub util_GetCentreOfMassFromSet{
   my ($pdb1,$isCAOnly,@res) = @_;
   
   my @allatoms ;
   foreach my $res (@res){
        if($isCAOnly){
           my $resnum = $res->GetResNum();
           my $CAatom1 = $pdb1->GetAtomFromResidueAndType($resnum,"CA");
           push @allatoms, $CAatom1 ;
        }
        else{
           my $resnum = $res->GetResNum();
           my @atoms1 = $res->GetAtoms();
           push @allatoms, @atoms1;
        }
   }
   my @X ; my @Y ; my @Z ;
   while(@allatoms){
       my $atom = shift @allatoms ;
       my ($x,$y,$z) = $atom->Coords();
       push @X,$x; push @Y,$y; push @Z,$z;
   }
   my ($meanX) = util_GetMeanSD(\@X);
   my ($meanY) = util_GetMeanSD(\@Y);
   my ($meanZ) = util_GetMeanSD(\@Z);
   return ($meanX,$meanY,$meanZ);
}

sub util_GetDistancesBetween2SetsOfResidues{
   my ($p1,$p2,$pdb1,$pdb2,$l1,$l2,$cutoff) = @_ ;
   my @allRes1 = @{$l1};
   my @allRes2 = @{$l2};
   my @allatoms1 ;
   my @allatoms2 ;
   foreach my $r (@allRes1){
        my @aaa = $r->GetAtoms();
        push @allatoms1, @aaa;
   }
   foreach my $r (@allRes2){
        my @aaa = $r->GetAtoms();
        push @allatoms2, @aaa;
   }
   return util_GetDistancesBetween2SetsOfAtoms($p1,$p2,$pdb1,$pdb2,\@allatoms1,\@allatoms2,$cutoff);
}

sub util_GetDistancesBetween2SetsOfAtoms{
   my ($p1,$p2,$pdb1,$pdb2,$l1,$l2,$cutoff) = @_ ;
   my @allatoms1 = @{$l1};
   my @allatoms2 = @{$l2};
   my $A = {};
   my $B = {};
   my $AB = {};
   my $sorted = {};
   while(@allatoms2){
               my $a2 = shift @allatoms2 ;
               foreach my $a1 (@allatoms1){
                   die if(!defined $pdb1);
                   my ($x1,$y1,$z1) = $a1->Coords();
                   my ($x2,$y2,$z2) = $a2->Coords();
                      my $d = util_format_float($a1->Distance($a2),3);
                   #my $d =  geom_Distance($x1,$y1,$z1,$x2,$y2,$z2) ;
                      next if($d > $cutoff);
                      my $nm1 = $a1->GetNameSlashSep();
                      my $nm2 = $a2->GetNameSlashSep();
   
                      my @l1 = split "/", $nm1 ;
                      my @l2 = split "/", $nm2 ;
   
                      # Ignore atoms that have H...
                      next if($l1[2] =~ /^H/i || $l2[2] =~ /^H/i);
   
                      my $N1 = $a1->GetResNum();
                      my $N2 = $a2->GetResNum();
   
                      my $X = $a1->GetResName();
                      my $Y = $a2->GetResName();
   
                      $A->{$N1} = $X ;
                      $B->{$N2} = $Y ;
                      $AB->{$X.$N1.$Y.$N2} = 1 ;
                      #next if($N1 eq $N2);
   
                      ####my $nm = "$nm1.$nm2";
                      my $nm = "$nm1 $nm2";
                      my $final = "$p1 $p2 $nm $d";
                      $sorted->{$final} = $d ;
               }
   
   }
   my $minvalue = 1000 ;
   foreach my $k (sort {$sorted->{$a} <=> $sorted->{$b}} keys %{$sorted}){
       $minvalue = $sorted->{$k};
       last ;
   }
   return ($sorted,$minvalue,$A,$B,$AB);
}



sub util_GetSequentialSetResidues{
    my ($pdb1,$start,$end) =@_ ;
           my $tableofres = {};
           my @listofres ;
           foreach my $i ($start..$end){
               $tableofres->{$i} = 1 ;
               my ($res) = $pdb1->GetResidueIdx($i);
			   ## broken 
			   if(!defined $res){
			   	  return 0;
			   }
			   my $x  = $res->PrintSingleLetter($pdb1);
               push @listofres, $res ;
           }
        return (1,$tableofres,@listofres);
}

sub util_PrintNewAtom{
   my ($X,$Y,$Z) = @_;
   my $atomstr = "HETATM";
   my $serialnum = 90000 + int(rand()*999) ;
   my $BUFF1 = " ";
   my $atomnm = "JUNK";
   my $alt_loc = " ";
   my $resname = "AAA";
   my $BUFF2 = " ";
   my $chainId = "A";
   my $resnum = 9000 + int(rand()*999) ;
   my $codeforinsertion = " ";
   my $BUFF3 = "   ";
   my $x = bufferStringWithSpaces($X,8);
   my $y = bufferStringWithSpaces($Y,8);
   my $z = bufferStringWithSpaces($Z,8);

   return  "${atomstr}${serialnum}${BUFF1}${atomnm}${alt_loc}${resname}${BUFF2}${chainId}${resnum}${codeforinsertion}${BUFF3}${x}${y}${z}\n";
}


sub bufferStringWithSpaces{
    my ($N,$bufflen) = @_ ;
    $N =~ s/ //g;
    my $LEN = length($N);
    die "$N ...$bufflen $LEN" if ($bufflen < $LEN);
    return $N if($bufflen eq $LEN);

    my $cnt = 0 ;
    while(length($N) < $bufflen){
         $N = " " . $N ;
         $cnt++;
    }
    #print "add $cnt\n";
    return $N ;
}
sub util_subtract2tables{
	my ($table1,$table2) = @_ ;
	my $table = {};
	foreach my $k (keys %{$table1}){
		if(! exists $table2->{$k}){
			$table->{$k} = 1 ;
		}
	}
	return $table ;
}
### given a set of ids - maps all to the first id...:w
## except the first, which is in the values
sub util_maptofirst{
	my ($fname) = @_ ;
	my $fh = util_read($fname);
	my $table = {};
	while(<$fh>){
	   next if(/^\s*$/);
	   next if(/^\s*#/);

		my (@l) = split ;
		my $first = shift @l ;
		foreach my $i (@l){
			$table->{$i} = $first ;
		}
	}
	return $table ;
}


sub util_HELIXParseAHBS{
  my ($infile) = @_ ;
  my $ifh = util_read($infile);
  my $infoAHBS = {};
  while(<$ifh>){
     next if(/^\s*$/);
	 next if(/^\s*#/);
	 next if(/NULL/);

	 my ($pdb1,$hnumber,$start,$end) = split ; 
	 $infoAHBS->{$pdb1} = {} if(!defined $infoAHBS->{$pdb1}) ;
	 foreach my $i ($start..$end){
	     $infoAHBS->{$pdb1}->{$i} = $hnumber ;
      }
  }	
  	return($infoAHBS);

}


## TODO: dont need the nletter variable - can parse out using string matching...
## like in createclaspinput..
sub util_SortResiduesBasedonNumber{
   my ($nletter,@l) = @_ ;
   my $t = {};
   foreach my $l (@l){
   	   my $orig = $l ;
	   if($nletter eq 3){
	       $l =~ s/...//;
	   }
	   else{
	       $l =~ s/.//;
	   }
	   $t->{$l} = $orig ;
   }
   my @ret ;
   foreach my $k (sort {$a <=> $b} keys %{$t}){
   	   my $v = $t->{$k} ;
	   push @ret, $v ;
   }
   return @ret ;
}

sub util_ParseContactFile{
	my ($X) = @_ ;
	my ($x,$p1,$p2,$A,$B,$H1,$H2,$res1,$res2) = split " ",$X ;
	return ($x,$p1,$p2,$A,$B,$H1,$H2,$res1,$res2) ;
}


## find different number of alphabets in a given string
sub util_FindAlphabetStrings{
	my ($j) = @_ ;
      my @l = split "",$j ;
      my $letters = {};
      foreach my $l (@l){
            $letters->{$l} = 1 ;
      }
      my $N = (keys %{$letters});
     return $N ;
}

sub util_FixFastaNames{
    my ($nm) = @_ ;
	$nm =~ s/\(/\./g;
	$nm =~ s/\)/\./g;
	$nm =~ s/\|/\./g;
	$nm =~ s/'/\./g;
	return ($nm);
}

##start KMERS 
sub util_SplitSliding{
    my ($string,$size,$name,$del,$kmertable,$expr2ignore) = @_ ;
	croak if(! defined $kmertable);
	my $len = length($string);
	return $kmertable if($len < $size);
	
	my $start = 0 ;
	my $end = $size ;
	my $mapped = 0 ;
	my $cnt = 0 ;
	while(1){
		$cnt++;
	    my $s1 = substr ($string , $start, $size)   ;
	    #next  if($s1 =~ /$expr2ignore/);


		## Why did I use cnt???
		#$kmertable->{$s1} = $cnt;

		$kmertable->{$s1} = [] if(! exists $kmertable->{$s1});
		push @{$kmertable->{$s1}}, $name ;

		last if($end eq $len);
		$start = $start + $del;
		$end = $start + $size;
	}
	return $kmertable;
}

## delfactor is the amount of hopping in terms of the ksize
sub util_SplitHopping{
    my ($str,$ksize,$delfactor,$name,$kmertable,$expr2ignore) = @_ ;
	my $len = length($str);

	return if($len < $ksize);

	my $start = 0 ;
	my $end = $ksize ;

	my $laststart = $len - $ksize ;

	## when hopping you mostly would need the name ...
	while(1){
	    my  $s1 = $end > $len ? substr ($str , $laststart) : substr ($str , $start, $ksize)   ;
	    #next  if($s1 =~ /$expr2ignore/);
		$kmertable->{$s1} = [] if(! exists $kmertable->{$s1});
		push @{$kmertable->{$s1}}, $name ;

		last if($end > $len );

		$start = $start + $ksize/$delfactor;
		$end = $start + $ksize/$delfactor;
	}
	return $kmertable ;
}

## Do we need this version?
sub util_splitstring{
      my ($line,$incr,$size) = @_ ; 
      my $len = length ($line);
      
      my $do = 1 ; 
      my $start = 0 ; 
      my $cnt = 0 ; 
      my @strings ;
      while($do){
           my $x = abs($start - $len) ;
           ## process last one
           if($x < $size){
               $size = $x  ;
               $do = 0 ; 
           }
           else{
                $cnt++ ; 
                my $s =  substr($line,$start,$size);
                push @strings, $s ;
           }
            $start = $start + $incr ;
      }
      return @strings ;
}

## KMERS end


##start PDB functions 

sub util_readPDB{
    my ($pdb1) = @_ ;
    my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP,$BLASTOUT,$BLASTDB) = util_SetEnvVars();
    my $origpdb = $pdb1 ;
    my $pdbfile = "$PDBDIR/$pdb1.pdb";
    my $pdb = new PDB($pdbfile);
    return $pdb ;
}

sub util_mapID2fullStringInFastaFile{
    my ($infile) = @_ ; 
    my $map = {} ; 
    my $info = {} ; 
    
	my $ifh = util_read($infile);
    while(<$ifh>){
		 chomp ;
         if(/^>/){
		 	my $orig = $_ ;
            s/^>//;
            my @l = split ;
            my $nm = $l[0];
            $nm = uc($nm);
            die if(!defined $nm); 
            warn "Repeat $nm " if(exists $info->{$nm});
            $info->{$nm} = "" ;
			$map->{$nm} = $orig ;
        }
	}
	close($ifh);
	return $map ;

}

sub util_GetOrf{
	my ($sequence,$maxnumber) = @_ ;
	croak if(!defined $maxnumber);
	my $retable = {};
	GetORFForSequence($retable,1,$sequence);
	my $complimentarysequence = util_getComplimentaryString($sequence) ;
	GetORFForSequence($retable,0,$complimentarysequence);
	my $cnt = 0 ;
	my $topORFs = {};
	foreach  my $key(sort {$retable->{$b} <=> $retable->{$a}} keys %{$retable}){
	    last if($cnt eq $maxnumber);
		$cnt++;
		my $length = $retable->{$key};
		my ($ntsequence,$fwdorrev,$start,$end) = split " ",$key ;
		$topORFs->{$cnt} = $key ;
		if($verbose){
		    #my $AA = util_ConvertNucleotide2Amino($ntsequence, 1);
		    #print "$cnt $ntsequence $length $AA start=$start end=$end\n";
		}
	}
	return $topORFs ;
}

sub GetORFForSequence{
	my ($retable,$fwdorrev,$sequence) = @_ ;
	print "Getting ORF for $sequence\n" if($verbose);
	$_ = $sequence ;
	my $ORIG = $sequence ;

	### Get the stop codons for each frame in one pass
	### for one direction - have to do twice for 
	### complimentary
	my @f1 ; my @f2 ; my @f3;
    while ( /T(?:AA|AG|GA)/g){
        my $pos = pos();

	    my $frame = $pos%3 ;
	    
	    if($frame eq 0){
	   	    push @f1, $pos ;
	    }
	    elsif($frame eq 1){
	   	    push @f2, $pos ;
	    }
	    else{
	   	    push @f3, $pos ;
	    }
    } 

	ProcessFrame($retable,$fwdorrev,$ORIG,0,@f1);
	ProcessFrame($retable,$fwdorrev,$ORIG,1,@f2);
	ProcessFrame($retable,$fwdorrev,$ORIG,2,@f3);
}

sub ProcessFrame{
	my ($retable,$fwdorrev,$ORIG,$frame,@f) = @_ ;
	my $len = length($ORIG);
	if(@f eq 0){
		my $TMP = $ORIG;
		if($frame eq 1){
			$TMP =~ s/.//;
		}
		elsif($frame eq 2){
			$TMP =~ s/..//;
		}
		my @l = ($TMP =~ /(...)/g);
		#my @ll = ($ORIG =~ /(...)/g);
		my $N = @l ;
		#my $NN = @ll ;
		#my $LLL =length($TMP);
		#print "KKKKKKKKKKK $N $NN $len $LLL\n";

		my $COMPLEXFACTOR = $frame eq 2? 3:2;
		my $pos = ($N+3) * 3 - $frame ;
		push @f, $pos ;

		my $description = $fwdorrev ? "FORWARD" : "REVERSE" ;
		print "frame $frame in $description has no stop codons. POS = $pos\n" if($verbose);
		die ;
	}
	my $prev = $frame ;
	foreach my $pos (@f){
	   	my $ntsequence = "XXXX";
	   	die if($prev eq -1);
		  
		## from one past the last stop codon to below 3 nt the next stop codon
		my $end = ($pos -3) ;
		my $start = $prev + 1 ;

		my $del = $end - $start ;
		if($del eq -1){
			 ## this means two end codons next to each other
             $prev = $pos ;
			 next ;
		}
	    $ntsequence = substr ($ORIG, $start -1 , $del+1);
	    if($fwdorrev eq 0){
		    $start = $len - $start + 1 ;
		    $end = $len - $end + 1 ;
	    }
	    print "fwdorrev=$fwdorrev start=$start end=$end frame=$frame prev=$prev , str=$ntsequence   \n" if($verbose);
		my $key = "$ntsequence $fwdorrev $start $end" ;
		$retable->{$key} = length($ntsequence);
        $prev = $pos ;
	}
}

sub util_ConcatStringNTimes{
	 my ($str,$N) = @_ ;
	 my $ret = "";
	 foreach my $i (1..$N){
	 	$ret = $ret . $str ;
	 }
	 return $ret ;
}
sub util_WriteORFSingleFasta{
	my ($nm,$sequence,$ofh,$ORIG,$maxnumber) = @_ ;
	### since we change $_ here 
	#print "Processing sequence $nm $sequence\n";
    my $topORFS = util_GetOrf($sequence,$maxnumber);
	my $N = keys %{$topORFS};
	foreach my $i (1..$N){
		my $v = $topORFS->{$i};
        my ($ntsequence,$fwdorrev,$start,$end) = split " ",$v ;
        my (@l) = split " ",$v ;
		if(@l ne 4){
			die "Error: somethign wrong not defined $nm $v \n";
		}
		my $AA = util_ConvertNucleotide2Amino($ntsequence, 1);
		my $description = $fwdorrev ? "FORWARD" : "REVERSE" ;
		print $ofh ">$nm.ORF_$i [$start-$end] $description \n";
		print $ofh "$AA \n";
	}
	return $ORIG ;
}

sub util_CheckSequence{
 	my ($name,$sequence,$isNT) = @_ ;
	$name =~ s/\(//g;
	$name =~ s/\)//g;
	#$nm =~ s/\.//g;
	$name =~ s/\|//g;
	$name =~ s/'//g;
	my $uc = uc($name);
	die "Need to make $name uppercase or make fix characters" if($uc ne $name);

    my $NElem = util_FindAlphabetStrings($sequence);
	if($isNT){
		    if($NElem > 6){
			    die "Error: isNT=$isNT, but NElem=$NElem for sequence=$sequence";
			}
		}
		else{
		    if($NElem < 5){
			    die "Error: isNT=$isNT, but NElem=$NElem for sequence=$sequence";
			}
	}
}
