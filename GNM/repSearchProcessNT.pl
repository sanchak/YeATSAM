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


use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($idx,$infile,$postfix,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$protein);
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
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);
#usage( "Need to give a input file name => option -infile ") if(!defined $infile);
#my $ifh = util_read($infile);
#usage( "Need to give a input file name => option -fastafile ") if(!defined $fastafile);
#my ($tablefasta,$infoSeq2PDB,$mapChainedName2Name) = util_parsePDBSEQRESNEW($fastafile,0,0);

usage( "Need to give a listfiles pdb id -option -listfiles  ") if(!@listfiles);
my $CNT = 0 ; 
#my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
 my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP,$BLASTOUT,$BLASTDB) = util_SetEnvVars();

# for reading fasta
#my ($fe,$x) =  util_readfasta($infile);
#chomp $x ;
my $PWD = cwd;


my $ignoretable ;
if(defined $ignorefile){
    ($ignoretable) = util_maketablefromfile($ignorefile);
}

my ($aadifftable) = util_maketablefromfile("AAD.txt");
my ($orflengthtable) = util_maketablefromfile("length.trsORF");

my $ofhpalinmap = util_write("palin.map");
my $ofhrepmap = util_write("rep.map");
my $ofhgcisland = util_write("gcisland");
my $ofhactuallyprotein = util_write("actuallyprotein");
my $ofhwarning = util_write("warning");

my $name2len = {};
my ($infoPalin0,$infoRep0) = ParseOneList($listfiles[0]);
my ($infoPalin1,$infoRep1) = ParseOneList($listfiles[1]);



print "Mapping palin\n";
my $mappedpalin = Map2Tables($infoPalin0,$infoPalin1,1);
print "Mapping Rep\n";
my $mappedrep = Map2Tables($infoRep0,$infoRep1,0);

system(" sortOnLast.pl -in palin.map -rev ");
system(" sortOnLast.pl -in rep.map -rev ");

my $DONEPRINTACTUALLY = {};
sub Map2Tables{
   my ($info0,$info1,$ispalin) = @_ ;
   my $mapped = {};
   foreach my $k (keys %{$info1}){
   	    my $tab1 = $info1->{$k} ;
		foreach my $kk (keys %{$tab1}){
	      if(exists $info0->{$k}){
		  	foreach my $v0 (keys %{$info0->{$k}}){
		    #my $v0 = $info0->{$k}->{$kk};
		       my $v1 = $kk ;
		       my $l0 = $name2len->{$v0};
		       my $l1 = $name2len->{$v1};
		       my $len = length($k);

			   

			   my $VAL0 = $info0->{$k}->{$v0};
			   my $VAL1 = $info1->{$k}->{$v1};
		       $mapped->{$v0} = $v1 ;

		       my $FH = $ispalin ?  $ofhpalinmap : $ofhrepmap ;
		       my $str2print = "mapped $v0 ($l0) to $v1 ($l1) through $k $VAL0 $VAL1 $len  \n";



			   die "" if(!exists $aadifftable->{$v1}) ;
			   die "$v1 not in orf " if(!exists $orflengthtable->{$v1});
			   my $aadiffval = $aadifftable->{$v1};
			   my $orflengthval = $orflengthtable->{$v1};
			   


			   my $v0ORG = $v0 ;
			   my $v1ORG = $v1 ;
			   $v0 =~ s/_G.*//;
			   $v0 =~ s/\s*//;

			   $v1 =~ s/\s*//;
			   $v1 =~ s/_G.*//;
			   #print "$v0 $v1 \n";
			   if($v1 eq $v0 || ($aadiffval < 45 && $orflengthval >100) ){
			   	    my $str = $v1 eq $v0 ? "eq" : " $aadiffval $orflengthval";
					if(! exists $DONEPRINTACTUALLY->{$v1}){
			   	        print $ofhactuallyprotein "$v1ORG aadiff=$aadiffval orf=$orflengthval $str2print";
						$DONEPRINTACTUALLY->{$v1}  = 1 ;
					}
			   }
			   else{
		            print $FH $str2print ;
			   }
			  }
	      }
	  }
   }
   return $mapped ;
}

sub ParseOneList{
	my ($i) = @_ ;
	my $l = ReadListFile($i);
	my $CNT = 0 ;
	my $CNTPALIN = 0 ;
	my $IGNPALIN = 0 ;
	my $CNTREP = 0;
	my $CNTREP60PLUS = 0;
    my $infoPalin = {};
    my $infoRep = {};
    foreach my $i (keys %{$l}){
		$i =~ s/_A//;
		$i =~ s/_B//;
		$i =~ s/_C//;
	    my $out = "REPEATINFO/$i.outfile";
	    my $palinfile = "REPEATINFO/$i.palin";
		if(! -e $out){
			print "$out not found \n";
			next ;
		}
		$CNT++;


		
	     my $ifh = util_read($out);

		 my $donestring = {};
		 while(<$ifh>){
                next if(/^\s*$/);
	            next if(/^\s*#/);
				my ($protein,$s,$len,$fulllen,$isPalin,$NF,$NR,$NUMFINAL,$LLL_l,$MMM_l ) = split ;
				next if(exists $donestring->{$s});
				my $rev = util_getComplimentaryString($s) ;
				$donestring->{$s} = 1 ;
				$donestring->{$rev} = 1 ;


				$name2len->{$protein} = $fulllen ;
				if($len =~ /^C/){
					die "$out .. $_ \n";
				}
				next if($s =~ /(ACACACAC|CACACACA|GAGAGAGA|ATATATAT|TATATATA|AAAAAA|TTTTTT|AGAGAGAG|GTGTGTGT|CTCTCTCT|TCTCTCTC)/);
				if($len >= 60){
                    $CNTREP60PLUS++;
				}
				else{
					$infoRep->{$s} = {} if(!defined $infoRep->{$s});
					$infoRep->{$s}->{$i} = $NUMFINAL;
				     $CNTREP++;
				}
				last ;
		 }
		 close($ifh);


		if(-e $palinfile){
	       my $ifh = util_read($palinfile);
		   my $internalpalincnt = 0 ;
		   
           while(<$ifh>){
                next if(/^\s*$/);
	            next if(/^\s*#/);
    
					s/fromend=//;
					s/fromstart=//;
	            my ($nm,$junk) = split ; 
				my ($protein,$s,$len,$NF,$fulllen,$LLL_l,$YYY) = split ;
				$name2len->{$protein} = $fulllen ;
				if($junk =~ /(ATATATAT|TATATATA|AAAAAA)/){
					my $min = $LLL_l > $YYY ? $YYY : $LLL_l ;
					if($min > 30){
					   print "$protein $s $min - ignoring, but not at edge \n" if($verbose);
					   $IGNPALIN++;
					}
				}
				elsif($junk =~ /(GCGCGC)/){
					   print $ofhgcisland  "$nm $junk \n";
				}
				else{
		           $internalpalincnt++;
				   $infoPalin->{$junk} = {} if(!defined $infoPalin->{$junk});
				   $infoPalin->{$junk}->{$nm} = 1 ;
				}
			 #print "$nm $junk \n";
           }
		   if($internalpalincnt){
		   	 $CNTPALIN++;
		   }
		   close($ifh);
		}
    }
	print "Found $CNT outfile with CNTREP $CNTREP, and CNTREP60PLUS $CNTREP60PLUS  and $CNTPALIN with palin, and IGNPALIN = $IGNPALIN\n";
	return ($infoPalin,$infoRep);
}

sub ReadListFile{
   my ($f) = @_ ; 
   my @list= util_read_list_sentences($f);
   my $list = {};
   map { s/\s*//g ; $list->{$_} = 1 ; } @list ;
   return $list ;
}

chmod 0777, $outfile ;

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
