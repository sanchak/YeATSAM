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
use MyGNM;


use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($idx,$infile,$listfile,$postfix,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$protein);
my ($ignorefile,@listfiles);
my $howmany = 100000 ;
my $verbose = 0 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "postfix=s"=>\$postfix ,
            "infile=s"=>\$infile ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "ignorefile=s"=>\$ignorefile ,
            "outfile=s"=>\$outfile ,
            "listfiles=s"=>\@listfiles,
            "howmany=i"=>\$howmany ,
            "idx=i"=>\$idx ,
            "verbose=i"=>\$verbose ,
            "cutoff=f"=>\$cutoff ,
           );
my $PWD = cwd;
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
#usage( "Need to give a input file name => option -infile ") if(!defined $infile);
#my $ifh = util_read($infile);
#usage( "Need to give a input file name => option -fastafile ") if(!defined $fastafile);
#my ($tablefasta,$infoSeq2PDB,$mapChainedName2Name) = util_parsePDBSEQRESNEW($fastafile,0,0);

 my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP,$BLASTOUT,$BLASTDB) = util_SetEnvVars();

my $MINCOUNT = 50 ;
my $MINLENGHT = 100 ;
my $MAXAADIFF = 45 ;
my $CNT = 0 ;

my $countStr = Config_getCountStr();
# for reading fasta
#my ($fe,$x) =  util_readfasta($infile);
#chomp $x ;


my ($tableTRS) = util_maketablefromfile("list.transcriptome.clean");
$CNT = (keys %{$tableTRS});
print "There are $CNT transcripts \n";


my ($tableaAADIFF) = util_maketablefromfile("DATA/AAD.txt");
my ($tableORFLength) = util_maketablefromfile("DATA/length.trsORF");
my ($tableRRNA) = util_maketablefromfile("DATA/RRNA.anno.1000.anno.real");

my ($tableTRSlength) = util_maketablefromfile("DATA/length.trs");
my ($tableScaffLength) = util_maketablefromfile("DATA/length.scaff");
my ($tableRepeatAA) = util_maketablefromfile("DATA/repeatAA.txt");
my ($infoAnno) = GNM_parseANNFile("genome.ann.realgood",$tableRRNA);
$CNT = (keys %{$infoAnno});
print "There are $CNT annotate TRS \n";
my ($tableNotanno) = util_subtract2tables($tableTRS,$infoAnno);
$CNT = (keys %{$tableNotanno});
print "There are $CNT non annotated TRS \n";

## 
my ($tableCountSummed,$allinfocounts) = GNM_parseCountSummed("bwa_counts_run1.txt.0.CLEAN.usingMeanSD.summed",$MINCOUNT);
my @sortedCounts = sort { $tableCountSummed->{$b} <=> $tableCountSummed->{$a}} (sort keys %{$tableCountSummed});

#save the complete count as a table with the TRS
my ($tableCount) = util_maketablefromfile_firstentry("bwa_counts_run1.txt.0.CLEAN.usingMeanSD");
$CNT = (keys %{$tableCountSummed});
print "There are $CNT trs with cnt > $MINCOUNT \n";


my $tableNovelgenes = {};
foreach my $trs (sort keys %{$tableNotanno}){
	#$trs =~ s/_A//;
	#$trs =~ s/_B//;
	#$trs =~ s/_C//;
	#next if (-e "REPEATINFO/$trs.palin" || ! -z "REPEATINFO/$trs.outfile");
	if(exists $tableCountSummed->{$trs} && (! exists $tableRepeatAA->{$trs})){
	
		my $aadiff = $tableaAADIFF->{$trs} ;
		my $lenORF = $tableORFLength->{$trs} ;
		my $lenTRS = $tableTRSlength->{$trs} ;
		my $hack ;
		if(!defined $lenORF){
			$hack = $trs . "_A";
		    $lenORF = $tableORFLength->{$hack} ;
		}
		die "$trs $hack" if(!defined $lenORF);
		if($aadiff <= $MAXAADIFF && $lenORF> $MINLENGHT){
		   $tableNovelgenes->{$trs}= 1 ;	
		}
	}
}
$CNT = (keys %{$tableNovelgenes});
print "There are $CNT novel genes \n";

my ($tableNCRNA) = {};
foreach my $trs (sort keys %{$tableNotanno}){
	#$trs =~ s/_A//;
	#$trs =~ s/_B//;
	#$trs =~ s/_C//;
	#next if (-e "REPEATINFO/$trs.palin" || ! -z "REPEATINFO/$trs.outfile");
	
		my $aadiff = $tableaAADIFF->{$trs} ;
		my $lenORF = $tableORFLength->{$trs} ;
		my $lenTRS = $tableTRSlength->{$trs} ;
		my $hack ;
		if(!defined $lenORF){
			$hack = $trs . "_A";
		    $lenORF = $tableORFLength->{$hack} ;
		}
		die "$trs $hack" if(!defined $lenORF);
	    if(exists $tableCountSummed->{$trs} && $aadiff > 45 && $lenORF < 80){
		   $tableNCRNA->{$trs}= 1 ;	
		}
}

my $WHAT = $tableNCRNA ;
PRINTRESULTS($tableNovelgenes,"novel","P");
PRINTRESULTS($tableNCRNA,"ncrna","N");
#my $WHAT = $tableNovelgenes ;


system("extractindexfromfile.pl -in novel.sortedMax -idx 1");
system("extractindexfromfile.pl -in ncrna.sortedMax -idx 1");
system("createTexTable.pl -in ncrna.sortedMax");
system("createTexTable.pl -in novel.sortedMax");

#### print counts ...

sub PRINTRESULTS{
    my ($WHAT,$name,$prefix) = @_ ; 
    my $ofhsorted = util_write("$name.sortedMax");
    print $ofhsorted "$countStr \n";
    my $ofhstat = util_write("$name.stat");
    print $ofhstat "MAX\n";
	my $CNTN = 0 ;

	my $DONETRSPREFIX = {};
    foreach my $trs (@sortedCounts){
    	if(exists $WHAT->{$trs}){
            my $v = $tableCount->{$trs};
    
			my $TTT = $trs ;
			$TTT =~ s/G.*//;
			next if(exists $DONETRSPREFIX->{$TTT});
			$DONETRSPREFIX->{$TTT} = 1;
			$CNTN++;
			my $ID = $prefix. $CNTN;
    	   GNM_printCountWithCutoff($ofhsorted,$v,$MINCOUNT,$ID);
    
            my $aadiff = $tableaAADIFF->{$trs} ;
            my $lenORF = $tableORFLength->{$trs} ;
    		my $lenTRS = $tableTRSlength->{$trs} ;
            if(!defined $lenORF || !defined $aadiff){
    		     my $hack = $trs . "_A";
    		     $lenORF = $tableORFLength->{$hack} ;
    		     $aadiff = $tableORFLength->{$hack} ;
            }
            if(!defined $lenORF || !defined $aadiff){
    			warn "$trs lllll\n";
    		}
            print $ofhstat "$ID $trs $aadiff $lenTRS $lenORF \n";
    	}
    }
    
    print $ofhstat "SINGULAR\n";
    
    my $ofh = util_write($name);
    foreach my $trs (sort keys %{$WHAT}){
       my $v = $tableCount->{$trs};
       print $ofh "$v\n";
       my $aadiff = $tableaAADIFF->{$trs} ;
       my $lenORF = $tableORFLength->{$trs} ;
    	my $lenTRS = $tableTRSlength->{$trs} ;
       if(!defined $lenORF ){
    		my $hack = $trs . "_A";
    		$lenORF = $tableORFLength->{$hack} ;
       }
       if(!defined $aadiff ){
    		my $hack = $trs . "_A";
    		$aadiff = $tableaAADIFF->{$hack} ;
       }
       if(!defined $aadiff ){
          die "$trs llll\b";
       }
       print $ofhstat "$trs $aadiff $lenTRS $lenORF\n";
    }
    
    
    
    system("countsTissueNumbers.pl -in $name -outf $name.sortedUNIQ -start 0");
}


###############


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
		 if($nm =~ /_A/){
		 	$nm =~ s/_A//;
		    $table->{$nm} = {};
		    $table->{$nm}->{EVAL} = $eval ;
		    $table->{$nm}->{BLASTSCORE} = $blastscore ;
		    $table->{$nm}->{PERCENT} = $percent ;
		 }
	}
	return $table ;
    
}


sub GNM_printCountWithCutoff{
	my ($ofhsorted,$v,$cutoff,$ID) = @_ ;
		my ($nm, @l) = split " ", $v ;
        print $ofhsorted "$ID ";
		foreach my $i (@l){
			my $VV = $i > $cutoff ? $i : 0 ;
            print $ofhsorted "$VV ";
		}
		print $ofhsorted "\n";
}
sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
