#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($all,$infile,$outfile,$or,$silent,$groupinfo);
my ($tag,$size,$DIR,$listfile,$ignorefile,$mergedir);
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
            "infile=s"=>\$infile ,
            "mergedir=s"=>\$mergedir ,
            "dir=s"=>\$DIR ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            "createnewname=s"=>\$createnewname ,
            "writedata=i"=>\$writedata ,
            "printbande=i"=>\$printbande ,
            "reverse=i"=>\$reverse ,
            "size=i"=>\$size ,
            "or=i"=>\$or ,
            "cutofflength=i"=>\$cutofflength ,
            "type=s"=>\@types,
            "protein=s"=>\$protein,
            "motif=s"=>\@motifs,
            "outfile=s"=>\$outfile 
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
#usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
$outfile = "$infile.uniq" if(!defined $outfile);
my $ofh = util_write("$outfile");
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT) = util_SetEnvVars();
print STDERR "Info: parsing file $infile - might take some time\n";
my ($info,$infoSeq2PDB,$mapChainedName2Name) = util_parsePDBSEQRESNEW($infile,0,$writedata,$reverse,$FASTADIR);
print "INFO: Parameters are writedata=$writedata, createnewname=$createnewname and reverse=$reverse,FASTADIR =$FASTADIR\n";
system("mkdir -p $FASTADIR");

my $ofhmappingname = util_write("$outfile.mappingname");
my $ofhmappingsame = util_write("$outfile.mappingsame");
my $ofhlength = util_write("$outfile.length");
my $ofhlist = util_write("$outfile.list");
my $ofhlistALL = util_write("$outfile.list.all");
my $ofhmappingALL = util_write("$outfile.mappingALL");

foreach my $i (sort keys %{$info}){
	print $ofhlistALL "$i\n";
}


my $ofhb100 ; my $ofhe100 ; my $ofhb30 ; my $ofhe30 ;

if(0){
   $ofhb100 = util_write("$outfile.begin.100");
   $ofhe100 = util_write("$outfile.end.100");
   $ofhb30 = util_write("$outfile.begin.30");
   $ofhe30 = util_write("$outfile.end.30");
}

## map seq and names to numbers
my $TRSID = 0 ;
my $MAPCNT2NAME = {};
my $MAPCNT2SEQ = {};
my $MAPSEQ2CNT = {};

my $shortlen = 0 ;

my $mapfullline = util_mapID2fullStringInFastaFile($infile);
foreach my $seq (sort keys %{$infoSeq2PDB}){
	$TRSID++;
	my $len = length($seq);

	if($cutofflength && $len < $cutofflength){
		$shortlen++;
		next ;
	}
	$TRSID++;

	my @l = sort @{$infoSeq2PDB->{$seq}};

	my $name = shift @l ;

	my ($fixName) = util_FixFastaNames($name);
	if($fixName ne $name){
		die "Error: Name $name contains illegal character, run \n fixFastaNm.pl -in $infile -out xxx \n";
	}


	## create new name? else mappingname remains the same
	my $ID = $createnewname ? $createnewname. "_" . $TRSID: $name ;


	print $ofhlength "$ID $len \n";
	print $ofhlist "$ID\n";

	die if(! exists $mapfullline->{$ID});
	my $mapfulllinechar = $mapfullline->{$ID};
	print $ofh "$mapfulllinechar\n"; ## no need for >
	print $ofh "$seq\n";
	print $ofhmappingname "$ID $name \n";


	## mutiple
	if(@l){
	   print $ofhmappingALL "$name @l $seq \n";


	   print $ofhmappingsame "$ID ";
	   foreach my $l (@l){
	       print $ofhlength "$l $len \n";
	       print $ofhmappingsame "$l\n";
	   }
	}

	if($printbande){
	    PrintBandE($ID,$seq,100,$ofhb100,$ofhe100);
	    PrintBandE($ID,$seq,30,$ofhb30,$ofhe30);
	}

}

sub PrintBandE{
	my ($ID,$seq,$N,$ofhb,$ofhe) = @_ ;
	my ($begin,$end) = util_GetTerminalStrings($seq,$N);
	print $ofhb ">$ID\n";
	print $ofhb "$begin\n";
	print $ofhe ">$ID\n";
	print $ofhe "$end\n";
}

print "Uniquifying $infile to $outfile, wrote $TRSID\n";
if($cutofflength){
	print "Ignore $shortlen due to shortlen\n";
}

system("sortOnLast.pl -rev -in $outfile.length");
system("wc -l $outfile.*");


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
