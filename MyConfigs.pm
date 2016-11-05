package MyConfigs;
use Carp ;
use POSIX ;
require Exporter;
use MyUtils;
no warnings 'redefine';
my $EPSILON = 0.01;

local $SIG{__WARN__} = sub {};

my $RAY = 1 ; 

@ISA = qw(Exporter);
@EXPORT = qw( 
	Config_GetALPHA 
	Config_Helix
	Config_GetUncharStrings
	Config_getCodonTable
	Config_tableforRevComl
	Config_getCountStr
	Config_getTissueList
	Config_NTFastaUnknown
	Config_AACodes
	Config_HetIgnore
	    );

use strict ;
use FileHandle ;
use Getopt::Long;

sub Config_GetALPHA{
   my @ALPHA = qw( A B C D E F G H I J K L M N O P Q R S T U V W X Y Z );
   return @ALPHA ;

}

sub Config_Helix{
   my $tableATOM = {};
   $tableATOM->{LYS} = "NZ";
   $tableATOM->{ARG} = "NH1";
   $tableATOM->{HIS} = "ND1";
   $tableATOM->{SER} = "OG";
   $tableATOM->{THR} = "OG1";
   $tableATOM->{CYS} = "SG";
   $tableATOM->{ASN} = "OD1";
   $tableATOM->{GLN} = "NE2";
   $tableATOM->{TYR} = "CB";
   $tableATOM->{PHE} = "CB";
   $tableATOM->{TRP} = "CB";
   $tableATOM->{ASP} = "OD1";
   $tableATOM->{GLU} = "OE1";
   $tableATOM->{GLY} = "CA";
   $tableATOM->{VAL} = "CB";
   $tableATOM->{ALA} = "CB";
   $tableATOM->{LEU} = "CB";
   $tableATOM->{MET} = "CB";
   $tableATOM->{ILE} = "CB";
   $tableATOM->{PRO} = "CB";

    my $HYDROVAL = {};
    ## from "Computer programs to identify and classify amphipathic a helical domains"
	#http://blanco.biomol.uci.edu/hydrophobicity_scales.html
	#$HYDROVAL->{ALA} = -0.17 ;
	#$HYDROVAL->{ARG} = -0.81 ;
	#$HYDROVAL->{ASN} = -0.42 ;
	#$HYDROVAL->{ASP} = -1.23 ;
	#$HYDROVAL->{CYS} = 0.24 ;
	#$HYDROVAL->{GLN} = -0.58 ;
	#$HYDROVAL->{GLU} = -2.02 ;
	#$HYDROVAL->{GLY} = -0.01 ;
	#$HYDROVAL->{HIS} = -0.96 ;
	#$HYDROVAL->{ILE} = -0.31 ;
	#$HYDROVAL->{LEU} = 0.56 ;
	#$HYDROVAL->{LYS} = -0.99 ;
	#$HYDROVAL->{MET} = 0.23 ;
	#$HYDROVAL->{PHE} = 1.13 ;
	#$HYDROVAL->{PRO} = -0.45 ;
	#$HYDROVAL->{SER} = -0.13 ;
	#$HYDROVAL->{THR} = -0.14 ;
	#$HYDROVAL->{TRP} = 1.85 ;
	#$HYDROVAL->{TYR} = 0.94 ;
	#$HYDROVAL->{VAL} = -0.07 ;


	## OLD VALUES 
    $HYDROVAL->{MET} = 0.975;
    $HYDROVAL->{ILE} = 0.913;
    $HYDROVAL->{LEU} = 0.852;
    $HYDROVAL->{VAL} = 0.811;
    $HYDROVAL->{CYS} = 0.689;
    $HYDROVAL->{ALA} = 0.607;
    $HYDROVAL->{THR} = 0.525;
    $HYDROVAL->{GLY} = 0.484;
    $HYDROVAL->{SER} = 0.402;
    $HYDROVAL->{HIS} = 0.333;
    $HYDROVAL->{PRO} = 0.239;
    $HYDROVAL->{PHE} = 1.036;
    $HYDROVAL->{TRP} = 0.668;
    $HYDROVAL->{TYR} = 0.137;
    $HYDROVAL->{GLN} = -0.558;
    $HYDROVAL->{ASN} = -0.701;
    $HYDROVAL->{GLU} = -1.396;
    $HYDROVAL->{LYS} = -1.518;
    $HYDROVAL->{ASP} = -1.600;
    $HYDROVAL->{ARG} = -2.233;

    my $colortable = {};
    my $value = {};
    my $chargedtable = {};
    my @NP = qw(G P A V L I M );
    foreach my $k (@NP){
    	$colortable->{$k} = "red";
    	$value->{$k} = 100 ;
    }
    
    my @pos = qw (H K R HIS LYS ARG);
    my @pos3 = qw (HIS LYS ARG);
    foreach my $k (@pos){
    	$chargedtable->{$k} = 1 ;
    	$colortable->{$k} = "blue";
    	$value->{$k} = 100 ;
    }
    
    my @neg = qw (E D ASP GLU); 
    foreach my $k (@neg){
    	$chargedtable->{$k} = -1 ;
    	$colortable->{$k} = "blue";
    	$value->{$k} = 60 ;
    }
    my @amide = qw(Q N );
    foreach my $k (@amide){
    	$colortable->{$k} = "blue";
    	$value->{$k} = 10 ;
    }
    
    my @misc = qw(F Y W C S T);
    foreach my $k (@misc){
    	$colortable->{$k} = "red";
    	$value->{$k} = 75 ;
    }

    my @aromatic = qw(F Y W TYR PHE TRP);
    my $aromatic = {};
    foreach my $k (@aromatic){
    	$aromatic->{$k} = "red";
    }

	return ($tableATOM,$HYDROVAL,$colortable,$value,$chargedtable,$aromatic);
}


sub Config_GetUncharStrings{
   #my $UNCHARSTRINGS = "clone:|bacterium|Uncultured bacterium|putative|unknown|chromosome|hypothe|unnamed|uncharacterized|Predicted protein";
   my $UNCHARSTRINGS = "variable|clone|bacterium|Uncultured bacterium|putative|unknown|chromosome|hypothe|unnamed|uncharacterized|Predicted protein";
   return $UNCHARSTRINGS ;
}

sub Config_getCodonTable{
my %DNA_code = (
'GCT' => 'A', 'GCC' => 'A', 'GCA' => 'A', 'GCG' => 'A', 'TTA' => 'L',
'TTG' => 'L', 'CTT' => 'L', 'CTC' => 'L', 'CTA' => 'L', 'CTG' => 'L',
'CGT' => 'R', 'CGC' => 'R', 'CGA' => 'R', 'CGG' => 'R', 'AGA' => 'R',
'AGG' => 'R', 'AAA' => 'K', 'AAG' => 'K', 'AAT' => 'N', 'AAC' => 'N',
'ATG' => 'M', 'GAT' => 'D', 'GAC' => 'D', 'TTT' => 'F', 'TTC' => 'F',
'TGT' => 'C', 'TGC' => 'C', 'CCT' => 'P', 'CCC' => 'P', 'CCA' => 'P',
'CCG' => 'P', 'CAA' => 'Q', 'CAG' => 'Q', 'TCT' => 'S', 'TCC' => 'S',
'TCA' => 'S', 'TCG' => 'S', 'AGT' => 'S', 'AGC' => 'S', 'GAA' => 'E',
'GAG' => 'E', 'ACT' => 'T', 'ACC' => 'T', 'ACA' => 'T', 'ACG' => 'T',
'GGT' => 'G', 'GGC' => 'G', 'GGA' => 'G', 'GGG' => 'G', 'TGG' => 'W',
'CAT' => 'H', 'CAC' => 'H', 'TAT' => 'Y', 'TAC' => 'Y', 'ATT' => 'I',
'ATC' => 'I', 'ATA' => 'I', 'GTT' => 'V', 'GTC' => 'V', 'GTA' => 'V',
'GTG' => 'V',);

return \%DNA_code ;

}

sub Config_getCountStr{
    return " TRS   CE   CI   CK   EM   FL   HC   HL   HP   HU   IF   LE   LM   LY   PK   PL   PT   RT   SE   TZ   VB";
}

sub Config_getTissueList{
    my @l = qw (CE   CI   CK   EM   FL   HC   HL   HP   HU   IF   LE   LM   LY   PK   PL   PT   RT   SE   TZ   VB);
    return @l ;
}


sub Config_NTFastaUnknown{
    return "R|Y|K|M|S|W|B|D|H|V|N";
}


sub Config_tableforRevComl{
    my $tableforRevComl = {};
    $tableforRevComl->{"A"} = "T";
    $tableforRevComl->{"T"} = "A";
    $tableforRevComl->{"G"} = "C";
    $tableforRevComl->{"C"} = "G";
    return $tableforRevComl ;
}




sub Config_AACodes{
   my ($tableThree2One,$tableOne2Three) ;

   $tableOne2Three->{"G"}  = "GLY";
   $tableOne2Three->{"P"}  = "PRO";
   $tableOne2Three->{"A"}  = "ALA";
   $tableOne2Three->{"V"}  = "VAL";
   $tableOne2Three->{"L"}  = "LEU";
   $tableOne2Three->{"I"}  = "ILE";
   $tableOne2Three->{"M"}  = "MET";
   $tableOne2Three->{"F"}  = "PHE";
   $tableOne2Three->{"Y"}  = "TYR";
   $tableOne2Three->{"W"}  = "TRP";
   $tableOne2Three->{"H"}  = "HIS";
   $tableOne2Three->{"K"}  = "LYS";
   $tableOne2Three->{"R"}  = "ARG";
   $tableOne2Three->{"Q"}  = "GLN";
   $tableOne2Three->{"N"}  = "ASN";
   $tableOne2Three->{"E"}  = "GLU";
   $tableOne2Three->{"D"}  = "ASP";
   $tableOne2Three->{"C"}  = "CYS";
   $tableOne2Three->{"S"}  = "SER";
   $tableOne2Three->{"T"}  = "THR";
   
   
   
   foreach my $k (keys %{$tableOne2Three}){
   	  my $v = $tableOne2Three->{$k};
   	  $tableThree2One->{$v} = $k ;
   }

   my @sortedSingle = (sort keys %{$tableOne2Three}) ;

   return ($tableThree2One,$tableOne2Three,@sortedSingle) ;

}


sub Config_HetIgnore{
   my @l = qw (HOH CE FE CD O PO4 FE2 CU NA MN ACN MG HG ZN CA BME ACT CL SO4 GOL MSE);
   my ($t) = MyUtils::util_make_table(\@l);
   return $t ;
}
