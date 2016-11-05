#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use MyConfigs;
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($writedata,$all,$infile,$outfile,$or,$silent,$groupinfo);
my ($DIR,$length,$listfile,$ignorefile);
my $howmany = 600000 ; 
my $cutofflength = 0 ; 
my @types = (); 
my $protein ;
my @motifs = (); 
GetOptions(
            "all"=>\$all ,
            "groupinfo"=>\$groupinfo ,
            "silent"=>\$silent ,
            "infile=s"=>\$infile ,
            "length=s"=>\$length ,
            "dir=s"=>\$DIR ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            "howmany=i"=>\$howmany ,
            "writedata=i"=>\$writedata ,
            "or=i"=>\$or ,
            "cutofflength=i"=>\$cutofflength ,
            "type=s"=>\@types,
            "protein=s"=>\$protein,
            "motif=s"=>\@motifs,
            "outfile=s"=>\$outfile 
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_open_or_append($outfile);
my $ofhmapping = util_write("$outfile.mapping");

print STDERR "Info: parsing file $infile - might take some time\n";
my ($info,$infoSeq2PDB,$mapChainedName2Name) = util_parsePDBSEQRESNEW($infile,0,$writedata);


my $ignoretable = {};
if(defined $ignorefile){
   my @lll= util_read_list_sentences($ignorefile);
   map { my @l = split ; $ignoretable->{$l[0]} = 1 ; } @lll ;
}


my ($tmp) = Config_getCodonTable();
my @AA = (values %{$tmp});
my $AA = util_make_table(\@AA);
my $NN = (keys %{$AA});

my $done = {};
my $total = 0 ;
my $final = 0 ;
my $ofhlength ;
$ofhlength = util_open_or_append($length) if(defined $length);
my $ign = 0 ;
my $processed = {};
foreach my $seq (keys %{$infoSeq2PDB}){
	my $len = length($seq);

	die if(exists $done->{$seq});
	my @l = @{$infoSeq2PDB->{$seq}};
	die "Expect no duplicates" if(@l > 1);
	my $trs = $l[0];
	if(exists $ignoretable->{$trs}){
		$ign++;
		next ;
	}

	my $CNTPERCENT = 0;
	foreach my $i (keys %{$AA}){
		my (@l) = ($seq =~ /$i/g);
		my $N = @l ;
		my $percent = int(100*$N/$len);
		if($percent > 20){
			print $ofh "$trs PERCENTSINGLE $i $percent\n";
			exit ;
		}
		elsif($percent > 15){
			$CNTPERCENT++;
			if($CNTPERCENT > 2){
			    print $ofh "$trs PERCENTMULT $i $percent\n";
			    exit ;
			}
		}
		#elsif($i eq "C" && $percent > 10){
			    #print $ofh "$trs PERCENTCYS $i $percent\n";
			    #exit ;
		#}
	}

	foreach my $i (keys %{$AA}){
		my $str = "";
		foreach my $j (6..10){
		    foreach my $idx (1..$j){
			    $str = $str . $i ;
		    }
		    if($seq =~ /$str/){
			    $processed->{$trs} = {} if(! defined $processed->{$trs});
			    $processed->{$trs}->{$i} = $str ;
		    }
		 }
	}
}

foreach my $seq (keys %{$infoSeq2PDB}){
	my $len = length($seq);

	die if(exists $done->{$seq});
	my @l = @{$infoSeq2PDB->{$seq}};
	die "Expect no duplicates" if(@l > 1);
	my $trs = $l[0];
	if(exists $ignoretable->{$trs}){
		$ign++;
		next ;
	}

	foreach my $i (keys %{$AA}){
		if(exists $processed->{$trs}){
			if( exists $processed->{$trs}->{$i}){
			     my $str = $processed->{$trs}->{$i} ;
			     $seq =~ s/$str/X/g;
			 
			}
		}
	}
	foreach my $i (keys %{$AA}){

		my $CNT = 0 ;
		my $JOINED = "";
		my $str = "";
		foreach my $j (3..5){
		    foreach my $idx (1..$j){
			    $str = $str . $i ;
		    }
		    if($seq =~ /$str/){
			   $CNT++;
			   $JOINED = $JOINED . $str;
		    }
		 }
		 if($CNT > 1){
		     $processed->{$trs} = {} if(! defined $processed->{$trs});
		     $processed->{$trs}->{$JOINED} = $JOINED ;
		}
	}
}

my $found = (keys %{$processed});
foreach my $trs (keys %{$processed}){
   my $tab = $processed->{$trs};
   print $ofh "$trs STRETCH ";
   foreach my $k (keys %{$tab}){
   	my $v = $tab->{$k};
	print $ofh " $v ";
   }
   print $ofh "\n";
}

$processed = (keys %{$infoSeq2PDB})  - $ign;
print "ign = $ign, processed = $processed, found =$found \n";

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
