#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($all,$infile,$outfile,$or,$silent,$groupinfo);
my ($in1,$in2,$size,$DIR,$listfile,$ignorefile,$mergedir);
my $howmany = 600000 ; 
my $cutofflength = 0 ; 
my @types = (); 
my $protein ;
my @motifs = (); 
my $createnewname = 0 ;
my $writedata = 0 ;
my $reverse = 0 ;
my $printbande = 0 ;
GetOptions(
            "all"=>\$all ,
            "groupinfo"=>\$groupinfo ,
            "silent"=>\$silent ,
            "in1=s"=>\$in1 ,
            "in2=s"=>\$in2 ,
            "mergedir=s"=>\$mergedir ,
            "dir=s"=>\$DIR ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            "createnewname=i"=>\$createnewname ,
            "printbande=i"=>\$printbande ,
            "size=i"=>\$size ,
            "or=i"=>\$or ,
            "cutofflength=i"=>\$cutofflength ,
            "type=s"=>\@types,
            "protein=s"=>\$protein,
            "motif=s"=>\@motifs,
            "outfile=s"=>\$outfile 
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a input file name => option -in1 ") if(!defined $in1);
usage( "Need to give a input file name => option -in2 ") if(!defined $in2);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write("$outfile");
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT) = util_SetEnvVars();
print "INFO: Parameters are writedata=$writedata, createnewname=$createnewname and reverse=$reverse,FASTADIR =$FASTADIR\n";
print STDERR "Info: parsing file $infile - might take some time\n";
my ($info1,$infoSeq2PDB1,$mapChainedName2Name1) = util_parsePDBSEQRESNEW($in1,0,$writedata,$reverse);
my ($info2,$infoSeq2PDB2,$mapChainedName2Name2) = util_parsePDBSEQRESNEW($in2,0,$writedata,$reverse);


my $LLL = 2 ;
my $done = {};
my $doneB = {};
my $doneE = {};
my $ofhmatches = util_write("MMM");
my $ofhALL = util_write("ALLWITHerrstoo");
my $ERR = {};
foreach my $LLL (0..10){
   my $size = 30 - $LLL ;
   my $table1 = util_GetUniqueEntries($infoSeq2PDB1,$LLL);
   my $table2 = util_GetUniqueEntries($infoSeq2PDB2,-1*$LLL);
   
   my $N1 = (keys %{$table1});
   my $N2 = (keys %{$table2});
   print "$N1 $N2 llllllllll\n";
   
   my $cnt = 0 ;
   foreach my $k (keys %{$table1}){
	   if(exists $table2->{$k}){
	   	   my @l1 = @{$table1->{$k}};
	   	   my @l2 = @{$table2->{$k}};
		   my $nm1 = $l1[0] ;
		   my $nm2 = $l2[0] ;
		   my $concat = $nm1 . " " .  $nm2 ;
		   next if(exists $done->{$concat});
		   next if(exists $ERR->{$nm1} ||exists $ERR->{$nm2});
		   my $str = " ";
		   if(exists $doneB->{$nm1} ||exists $doneE->{$nm2}){
		   	  $ERR->{$nm1}  =1 ;
		   	  $ERR->{$nm2}  =1 ;
			  $str = "ERR";
		   }
	       print $ofhALL "$size $nm1 $nm2 $str \n";
		   $doneB->{$nm1} = 1 ;
		   $doneE->{$nm2} = 1 ;
		   $done->{$concat} = $size ;
		   $cnt++;
	   }
   }
   print "There were $cnt unique matches ...\n";
}

my $ignorecnt = 0 ;
foreach my $k (keys %{$done}){
	my ($nm1,$nm2) = split " ", $k ;
	if(exists $ERR->{$nm1} ||exists $ERR->{$nm2}){
		$ignorecnt++;
	}
	my $size = $done->{$k};
	print $ofhmatches "$size $nm1 $nm2 \n";
}
print STDERR "Ignored $ignorecnt\n";


## Get entries that have only match
## remove $N chars from the begin or end - depending on the sign
sub util_GetUniqueEntries{
	my ($infoSeq2PDB,$N) = @_ ;
	my $ret  = {};
    foreach my $seq (keys %{$infoSeq2PDB}){
	    my @l = sort @{$infoSeq2PDB->{$seq}};
		my $val = abs($N);
		if($N > 0){
			foreach my $id (1..$val){
			    $seq =~ s/.$//;
			}
		}
		else{
			foreach my $id (1..$val){
			    $seq =~ s/^.//;
			}
		}
		my $L = length($seq);
	    if(@l eq 1){
			$ret->{$seq} = \@l ;
	    }
	}
	return $ret ;
}


system("wc -l $outfile.*");


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
