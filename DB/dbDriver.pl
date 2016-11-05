#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($all,$infile,$outfile,$or,$silent,$groupinfo);
my ($DIR);
my $howmany = 600000 ; 
my $cutofflength = 0 ; 
my @types = (); 
my @ntypes = (); 
my @motifs = (); 
GetOptions(
            "all"=>\$all ,
            "groupinfo"=>\$groupinfo ,
            "silent"=>\$silent ,
            "infile=s"=>\$infile ,
            "dir=s"=>\$DIR ,
            "howmany=i"=>\$howmany ,
            "or=i"=>\$or ,
            "cutofflength=i"=>\$cutofflength ,
            "type=s"=>\@types,
            "ntype=s"=>\@ntypes,
            "motif=s"=>\@motifs,
            "outfile=s"=>\$outfile 
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
usage( "Need to give a dir => option -dir ") if(!defined $DIR);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);

my $CNT = 0 ; 

my ($RESULTDIR,$PDBDIR) = util_SetEnvVars();
system ("mkdir -p PDBINFO");
print "Info: $PDBDIR is PDBDIR\n";


my $seqlist = {};

my $chain5 = {};
my $info = {};
my $chain5protein = {};
my $ifh = util_read($infile);

## certain PDBs cant be extracted
my $ofhproblems = util_open_or_append("$outfile.problems");
my ($problemtable) = util_maketablefromfile("$outfile.problems");

print "Reading the $infile\n";
while(<$ifh>){
    next if(/^\s*$/);
	if(/^\s*>/){
	        my ($nm,$origchain,$type,$len,$fullnm) = parseSingleLine_pdb_seqres($_);
		    die if(!defined $nm); 
		    $nm = uc($nm);

			# TODO - make function
			if(! exists $problemtable->{$nm}){
			    print "READING PDB to find problem for $nm\n";
			    my $PDBhasproblem = 0 ;
	            if( -e "$PDBDIR/$nm.pdb"){
	                my $IFH = util_read("$PDBDIR/$nm.pdb");
				    while(<$IFH>){
					    if(/^HEADER/ || /^ATOM/ || /^REMARK/) {
						    last ;
					    }
					    elsif (/but the requested file is not available/){
						$PDBhasproblem = 1 ;
						    last ;
					    }
					    
				    }
				    close($IFH);
	            }
			    print $ofhproblems "$nm $PDBhasproblem\n";
			    $problemtable->{$nm} = $PDBhasproblem ;
			}
            next if(exists $problemtable->{$nm} && $problemtable->{$nm} eq 1);



			my $chain = ($origchain);
		    my $nmwithchain = $nm . $chain ;
		    next if(exists $chain5->{$nmwithchain}); 

	        my $seq = <$ifh> ;
		    chomp $seq ; 

		    $chain5->{$nmwithchain} = $seq ;
			if(($type =~ /protein/i)){
				$chain5protein->{$nmwithchain} =$seq ;
			}

			$CNT++ ; 
			if(! defined $info->{$nm}){
		        $info->{$nm} = {}  ;
		        $info->{$nm}->{CHAINS} = [] ;
		        $info->{$nm}->{FULLNM} = $fullnm ;
			}
			push @{$info->{$nm}->{CHAINS}}, $nmwithchain ;
    }
}


my $fiveletterN = keys %{$chain5};
my $fourletterN = keys %{$info};
print "Info: There are fiveletterN $fiveletterN and fourletterN $fourletterN\n";

my $ofhgetpdbs = util_write("$outfile.step0.getpdbs.csh");
my $ofhsplit = util_write("$outfile.step1.split.csh");
my $ofhahbs = util_write("$outfile.step2.ahbs.csh");
my $ofhfinddist = util_write("$outfile.step3.finddist.csh");
my $ofhhetatms = util_write("$outfile.step4.hetatm.csh");

util_writelist2file("$outfile.list5",sort keys %{$chain5});
util_writelist2file("$outfile.list4",sort keys %{$info});

my $ofhlist5protein = util_write("$outfile.list5.protein");
my $ofhlist5proteinnotthere = util_write("$outfile.list5.protein.notthere");
my $ofhlist5proteinseq = util_write("$outfile.list5.protein.seq.uniq");

#####
### Write the "uniq" sequences - also generate the mapping
#####
my $UNIQSEQ = {};
foreach my $i (sort keys %{$chain5protein}){
	print $ofhlist5protein "$i\n";
	if(! -e "$PDBDIR/$i.pdb"){
	    print $ofhlist5proteinnotthere "$i\n";
	}

	## Write the sequences
	my $seq = $chain5protein->{$i};
	if(!exists $UNIQSEQ->{$seq}){
		$UNIQSEQ->{$seq} = [];
	    print $ofhlist5proteinseq ">$i\n";
	    print $ofhlist5proteinseq "$seq\n";
	}
	push @{$UNIQSEQ->{$seq}},$i ; 
}

## Exact mapping
my $ofhlist5proteinseqmapping = util_write("$outfile.list5.protein.seq.uniq.mapping");
foreach my $seq (keys %{$UNIQSEQ}){
     my @l = @{$UNIQSEQ->{$seq}};
	 next if(@l eq 1);
	 print $ofhlist5proteinseqmapping "@l \n";
}



### Generate the commands
print $ofhgetpdbs  "cd $PDBDIR\n";
foreach my $nm (sort keys %{$info}){
    if(! -e "$PDBDIR/$nm.pdb"){
	   print $ofhgetpdbs  "wget   http://www.rcsb.org/pdb/files/$nm.pdb \n";
	}
	print $ofhsplit  "pdbgetlistAndSplit.csh $nm \n";
	print $ofhahbs  "pdbAPBS.csh $nm \n";


	my @chains = @{$info->{$nm}->{CHAINS}};
	my $ofhlist = util_write("$PDBDIR/$nm.list");
	foreach my $i (@chains){
	     print $ofhlist "$i\n";
	}

	if(@chains > 1){
	   print $ofhfinddist  "pdbgetPWDistance.csh $nm\n";
	}
	my $anno = $info->{$nm}->{FULLNM} ;
	$, = "  " ;
	#print $ofh "@chains $anno\n";
}

print $ofhhetatms "hetatmALLPDBS.pl -lis $outfile.list5.protein\n";

system("wc -l $outfile.*");

system("makepdbfourletterfromfive.pl -inf $outfile.list5.protein.notthere");
system("UNIQ $outfile.list5.protein.notthere.four");

sub parseSingleLine_pdb_seqres{
	my ($line) = @_ ; 
	my ($nm,$chain,$type,$len,$fullnm) = ($line =~ /^.(....).(.*)\s\s*mol:(\w+)\s*length:(\d+)\s*(.*)/);
	
	return ($nm,$chain,$type,$len,$fullnm) ;

}
sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
