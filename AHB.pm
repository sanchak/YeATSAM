package AHB;
use Carp ;
use POSIX ;
require Exporter;
use MyUtils;
use MyGeom;
no warnings 'redefine';
my $EPSILON = 0.01;
use Algorithm::Combinatorics qw(combinations permutations) ;
local $SIG{__WARN__} = sub {};

my $RAY = 1 ; 

@ISA = qw(Exporter);
@EXPORT = qw( 
	
		ahb_GetPermuations
		ahb_printPairWiseDistance
		ahb_PrintData
		ahb_printCrossWiseDistance
		ahb_ParseValues
	    );
sub ahb_GetPermuations{
	my ($what,$list,$info,$num) = @_ ;
    my $iter = combinations($list, $num);
    my $cnt = 0 ;
	my @ret ;
    while (my $c = $iter->next) {
          my @combo = @{$c} ; 
          my $iter1 = permutations(\@combo, $num);
          while (my $c1 = $iter1->next) {
             my @perm = @{$c1} ; 
		     my $idx = 1 ;
			 my @lll ;
			 my @COORDS ;
		     foreach my $a (@perm){
		        my $nm = "$what". $idx++;
		 	    my $str = "$nm $info->{$a}";
			    my @COORD = split " ", $info->{$a} ;
				push @COORDS, \@COORD;
				push @lll, $str ;
		       }
			 my $toofar = 0 ;
			 while(@COORDS){
			 	my $A = shift @COORDS;
				foreach my $B (@COORDS){
		           my $dist =  geom_Distance(@{$A},@{$B});
				   if($dist > 20){
				      $toofar = 1 ;
				   }
				}
			 }

			  if(!$toofar){
			     push @ret, \@lll ;
			  }

	      }
    }
	return @ret ;
}

sub ahb_printPairWiseDistance{
   my ($protein,$pdb1,$RESTABLE,@retB) = @_ ; 
   my $info  = {};
   while (@retB){
	   my ($B1,$x1,$y1,$z1) = split " ", shift @retB;

	   my $L1 = $RESTABLE->{$B1};
	   foreach my $XXX (@retB){
	       my ($B2,$x2,$y2,$z2) = split " ", $XXX;
		   my $dist =  geom_Distance($x1,$y1,$z1,$x2,$y2,$z2) ;
	       my $L2 = $RESTABLE->{$B2};
		   my ($sorted,$minvalue,$A,$B,$AB) = 
		               util_GetDistancesBetween2SetsOfResidues($protein,$protein,$pdb1,$pdb1,$L1,$L2,10);
		   if($minvalue < 4){
		   	   my $str = "DIST $B1 $B2 $dist $minvalue";
		   	   my $key = $B1.$B2 ;
		       $info->{$key} = $str ;
		   }
	   }
   }
   return $info ;
}

sub ahb_printCrossWiseDistance{
   my ($protein,$pdb1,$RESTABLE,$retB,$retA) = @_ ; 
   my @retA = @{$retA};
   my @retB = @{$retB};
   my $info  = {};
   while (@retB){
	   my ($B1,$x1,$y1,$z1) = split " ", shift @retB;
	   my $L1 = $RESTABLE->{$B1};
	   foreach my $XXX (@retA){
	       my ($B2,$x2,$y2,$z2) = split " ", $XXX;
		   my $dist =  geom_Distance($x1,$y1,$z1,$x2,$y2,$z2) ;
	       my $L2 = $RESTABLE->{$B2};
		   my ($sorted,$minvalue,$A,$B,$AB) = 
		               util_GetDistancesBetween2SetsOfResidues($protein,$protein,$pdb1,$pdb1,$L1,$L2,10);
		   if($minvalue < 4){
		   	   my $str = "DIST $B1 $B2 $dist $minvalue";
		   	   my $key = $B1.$B2 ;
		       $info->{$key} = $str ;
		   }
	   }
   }
   return $info ;
}



sub ahb_ParseValues{
	my ($protein,$pdb1,$values,$isCAOnly) = @_ ;
   my $RESTABLE = {};
   my @AHList ;
   my @BList ;
   my $DISUPLH = {};
   my $AHinfo = {};
   my $Binfo = {};
   my $ifh = util_read($values);
   my $HM = {};
   while(<$ifh>){
		my ($p,$id,$start,$end,$len,$hm) = split ;
		next if($p ne $protein);

		my ($tableofres,@listofres) = util_GetSequentialSetResidues($pdb1,$start,$end) ;
		foreach my $i ($start..$end){
			$DISUPLH->{$i} = $id ;
		}
        my ($meanX,$meanY,$meanZ) = util_GetCentreOfMassFromSet($pdb1,$isCAOnly,@listofres);
		$RESTABLE->{$id} = \@listofres;
		#print "$start $end \n";
		#print "llll $id $meanX,$meanY,$meanZ\n";
		my $str = " $id $meanX $meanY $meanZ ";
		if($id =~ /HELIX/){
			$HM->{$id} = $hm ;
			#print $ofh "HM $id $hm \n";
			push @AHList,$str ;
			$AHinfo->{$id} = " $meanX $meanY $meanZ ";
		}
		else{
			push @BList,$str ;
			$Binfo->{$id} = " $meanX $meanY $meanZ ";
		}
   }
   return ($RESTABLE,$AHinfo,\@AHList,$Binfo,\@BList,$DISUPLH,$HM);
}

sub ahb_PrintData{
	my ($info,$ofh) = @_ ;
	foreach my $k (keys %{$info}){
		my $v = $info->{$k} ;
		print $ofh "$v\n";
	}
}
