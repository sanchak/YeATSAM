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
my ($fastafile,$idx,$infile,$postfix,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$protein);
my ($ignorefile,@expressions);
my $howmany = 100000 ;
my $repeatdir = "REPEATINFO" ;
my $verbose = 0 ;
GetOptions(
            "repeatdir=s"=>\$repeatdir ,
            "fastafile=s"=>\$fastafile ,
            "protein=s"=>\$protein ,
            "postfix=s"=>\$postfix ,
            "infile=s"=>\$infile ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "idx=i"=>\$idx ,
            "verbose=i"=>\$verbose ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);
usage( "Need to give a input file name => option -listfile ") if(!defined $listfile);
usage( "Need to give a input file name => option -fastafile ") if(!defined $fastafile);


system("mkdir -p  $repeatdir");
#my $STARTLENGHT = 21 ; 
my $STARTLENGHT = 21 ; 
my $FINALLENGTH = 6 ; 

print STDERR "Reading fastafile ...will do from $STARTLENGHT to $FINALLENGTH \n";
my ($tablefasta,$infoSeq2PDB,$mapChainedName2Name) = util_parsePDBSEQRESNEW($fastafile,0,0);

#my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
 my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP,$BLASTOUT,$BLASTDB) = util_SetEnvVars();

my @list ; 
if( -e $listfile){
     @list = util_read_list_words($listfile);
}
else{
   die "$listfile does not exist in the fasta" if(! -exists $tablefasta->{$listfile});
   push @list, $listfile ;
}


my $map2str = {};
my $map2strrev = {};
foreach my  $i (@list){
   die "$i does not exist" if(! -exists $tablefasta->{$i});
   ## reintialize all these lines
   my $ofhprotein = util_write("$repeatdir/$i.outfile");
   system("rm -f $repeatdir/$i.palin");
   my $fe = $tablefasta->{$i} ;
   $map2str->{$i} = $fe ;
   my $rev = util_getComplimentaryString($fe) ;
   $map2strrev->{$i} = $rev ;
}
print "Done storing...\n";

my $DONE = {};

my $NUMOFMATCHREQ = 1 ;
my $palintable = {};

my @NUMBERS = qw (60 30);
while($STARTLENGHT > $FINALLENGTH){
	$STARTLENGHT--;
	push @NUMBERS, $STARTLENGHT;
}

while(@NUMBERS){
	$STARTLENGHT = shift @NUMBERS ;

	## order is imp :)
	if($STARTLENGHT < 6){
		$NUMOFMATCHREQ = 40 ;
	}
	elsif($STARTLENGHT < 8){
		$NUMOFMATCHREQ = 8 ;
	}
	elsif($STARTLENGHT < 10){
		$NUMOFMATCHREQ = 4 ;
	}
	elsif($STARTLENGHT < 13){
		$NUMOFMATCHREQ = 2 ;
	}
  
  foreach my $protein (@list){
      print "processing $protein for $STARTLENGHT\n" if($verbose);
	  $palintable->{$protein} = "" if(!defined $palintable->{$protein});

	  next if(exists $DONE->{$protein});

	  my $fe = $map2str->{$protein} ;
	  my $rev = $map2strrev->{$protein} ;
	  my $fulllen = length($fe);

      my @strings = util_splitstring($fe,1,$STARTLENGHT);

	  my $donestring = {};
	  
	  my $palinstrings = $palintable->{$protein};
	  my $table2print = {};
      foreach my $s (@strings){
	      my $len = length($s);
	      next if( exists $donestring->{$s});
	      $donestring->{$s} = 1 ;
		  my $revstr = util_getComplimentaryString($s);
		  my $isPalin = 0 ;
          my $NF = ProcessSingleString($fe,$s,$DONE,$STARTLENGHT,$NUMOFMATCHREQ,$protein,"F");
          my $NR = ProcessSingleString($rev,$s,$DONE,$STARTLENGHT,$NUMOFMATCHREQ,$protein,"R");

		  my ($LLL) = ($fe =~ /(.*)$s/);
		  my $LLL_l = length($LLL);
		  my ($MMM) ;
		  if($rev =~ /$s/){
		      ($MMM) = ($rev =~ /(.*)$s/);
		  }

	      if($s eq $revstr){
	  	     $isPalin = 1 ;

			 ## ignore subsets - that have only one occurance
			 if($palinstrings =~ /$s/ && ($NF eq $NR) && ($NF eq 1)){
			 	   next ;
		     }

			 $palinstrings = $palinstrings . "XXX" . $s ;
			 ## write palins only if >= 14
			 if($len > 13){
		         my $ofhpalin = util_open_or_append("$repeatdir/$protein.palin");
			     my $YYY = $fulllen - $LLL_l ;
			     print $ofhpalin "$protein $s $len $NF $fulllen fromstart=$LLL_l fromend=$YYY \n";
			 }
	      }
		  my $NUMFINAL = $isPalin ? $NF : $NF+$NR;
	      if($NUMFINAL > $NUMOFMATCHREQ){
			  my $MMM_l = defined $MMM ?  $fulllen - length($MMM) : -11111111;
			  $table2print->{$s} = "$protein $s $len $fulllen $isPalin $NF $NR $NUMFINAL $LLL_l $MMM_l \n";
	     }
     } ## foreach string

	 ## printing 
	 foreach my $s (keys %{$table2print}){

	 	  my $v = $table2print->{$s};
		  my ($protein,$s,$len,$fulllen,$isPalin,$NF,$NR,$NUMFINAL,$LLL_l,$MMM_l ) = split " " ,$v ;
		  #print "$s  $palinstrings\n";
		  #$s =~ s/^.//;
		  #$s =~ s/.$//;
		  if($palinstrings =~ /$s/ && ($NR eq $NR) && ($NF eq 1)){
				 	next ;
		  }
	      my $ofhprotein = util_open_or_append("$repeatdir/$protein.outfile");
	      print $ofhprotein "$v";
	 }
	 $palintable->{$protein} = $palinstrings;
  } ## forrach protein

}

sub ProcessSingleString{
	my ($fe,$s,$DONE,$LEN,$NUMOFMATCHREQ,$protein,$what) =@_ ;
	my @l = ($fe =~ /($s)/g);
	my $N = @l ;
	return $N ;
}

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
