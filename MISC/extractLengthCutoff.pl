#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($cutoff,$all,$infile,$outfile,$or,$silent,$groupinfo);
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
            "cutoff=i"=>\$cutoff ,
            "or=i"=>\$or ,
            "type=s"=>\@types,
            "protein=s"=>\$protein,
            "motif=s"=>\@motifs,
            "outfile=s"=>\$outfile 
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
usage( "Need to give a cutoff file name => option -cutoff ") if(!defined $cutoff);
#usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
$outfile = "$infile.cutoff.$cutoff" if(!defined $outfile);
my $ofh = util_write("$outfile");
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT) = util_SetEnvVars();
print STDERR "Info: parsing file $infile - might take some time\n";
my ($info,$infoSeq2PDB,$mapChainedName2Name) = util_parsePDBSEQRESNEW($infile,0,$writedata,$reverse,$FASTADIR);

my $mapfullline = util_mapID2fullStringInFastaFile($infile);
my $done = {};
my $nALL = 0 ;
my $nREM = 0 ;
foreach my $name (sort keys %{$info}){
	$nALL++;
	my $seq = $info->{$name};
	next if($done->{$name});
	$done->{$seq} = 1 ;
	my $len = length($seq);
	next if($len < $cutoff);
	die if(! exists $mapfullline->{$name});
	my $mapfulllinechar = $mapfullline->{$name};
	$nREM++;
	print $ofh ">$mapfulllinechar\n";
	print $ofh "$seq\n";
}

print "There were $nALL, out of which wrote $nREM\n";

system("wc -l $outfile*");


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
