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
usage( "Need to give a list file name => option -list ") if(!defined $listfile);
$outfile = "$listfile.fa" if(!defined $outfile);
my $ofh = util_write("$outfile");
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT) = util_SetEnvVars();
#print STDERR "Info: parsing file $infile - might take some time\n";
#my ($info,$infoSeq2PDB,$mapChainedName2Name) = util_parsePDBSEQRESNEW($infile,0,$writedata,$reverse,$FASTADIR);

my @list= util_read_list_sentences($listfile);
my $list = {};
map { s/\s*//g ; $list->{$_} = 1 ; } @list ;

my $mapfullline = util_mapID2fullStringInFastaFile($infile);
my $done = {};
my $nALL = 0 ;
my $nREM = 0 ;

my ($NAME,$sequence);
my $seenNAMES = {};
my $ifh = util_read($infile);
while(<$ifh>){
	     next if(/^\s*$/);
         if(/^\s*>/){
	 	    ## process one fasta
	 	    if(defined $NAME){
	           $nALL++;
			   if(! exists $seenNAMES->{$NAME} ){
	               if(exists $list->{$NAME}){
				   	   die "$NAME does not exist in fullline" if(!exists $mapfullline->{$NAME});
				   	   my $fullline = $mapfullline->{$NAME} ;
					   #print $ofh "$fullline\n";
					   #print $ofh "$sequence\n";
				       $nREM++;

				   }
			       $seenNAMES->{$NAME} = $sequence ;
			   }
			   else{
			        warn "Warning: repeat name $NAME\n" ;
			   }
               $sequence = "";
		    }
    
		    s/>//;
		    ($NAME) = split ;
	     }
	     else{
	        die if(/>/);
		    ## remove spaces 
		    s/ //g;
		    chomp;
		    $sequence = $sequence . $_ ;
	     }
}
if(! exists $seenNAMES->{$NAME} ){
	$seenNAMES->{$NAME} = $sequence ;
	      $nALL++;
	      if(exists $list->{$NAME}){
				die if(!exists $mapfullline->{$NAME});
				my $fullline = $mapfullline->{$NAME} ;
				#print $ofh "$fullline\n";
				#print $ofh "$sequence\n";
				$nREM++;

		}
}
close($ifh);

foreach my $i (@list){
	if(! -exists $seenNAMES->{$i}){
		die " $i not found \n";
	}
	else{
				my $fullline = $mapfullline->{$i} ;
				my $sequence = $seenNAMES->{$i} ;
				print $ofh "$fullline\n";
				print $ofh "$sequence\n";
		
	}
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
