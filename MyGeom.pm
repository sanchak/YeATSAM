
package MyGeom;
use MyUtils;
use Atom;
use Residue;
use Carp ;
use POSIX ;
use Algorithm::Combinatorics qw(combinations) ;
require Exporter;
use Math::Geometry ; 
#use Math::Trig;
#use Math::Trig ':radial';
@ISA = qw(Exporter);
@EXPORT = qw( 
	geom_R geom_DirCosines geom_Distance geom_AngleBetween2Points
	geom_AngleBetweenThreePoints geom_TransformAxis
	geom_IsZero geom_AlignXaxis2Vector
	geom_MidPointPoints geom_MidPointAtoms
	geom_PointsOnACircle geom_Makepoint
	geom_Closest2
	geom_RotateWithXConstantAndZwillbeZero
	geom_Align3PointsToXYPlane
	geom_Distance_2D geom_Distance_3D
	geom_GetABForEllipse
	geom_GetPointsBetween
	geom_ExtendTwoPointsDouble
	geom_GetCircleAroundPoint
	geom_GetPointsAroundOnePoint 
	geom_MakeKeyFromCoord
	geom_GetPointsTopDownNorthSouthEastWest
	geom_MakeCoordFromKey

	MakeVectorFrom2Points_Atoms
	MakeVectorFrom2Points


	geom_AngleBetweenTwoVectors
	geom_AngleBetweenTwoVectors_4Atoms


	geom_GetVectorAtAngle
	geom_GetSetOfVectorsForFullCircle
	    );

use strict ;
use FileHandle ;
use Getopt::Long;

my $EPSILON = 0.01;


sub geom_R{
	my ($x,$y,$z) = @_; 
	my $s = $x*$x + $y*$y + $z*$z ; 
	return sqrt($s) ; 
}

sub geom_MidPointAtoms{
	my ($a,$b) = @_; 
	return geom_MidPointPoints($a->Coords(),$b->Coords());
}


sub geom_MidPointPoints{
	my ($x,$y,$z,$a,$b,$c) = @_; 
	my $p = ($x+$a)/2 ;
	my $q = ($y+$b)/2 ;
	my $r = ($z+$c)/2 ;
	return ($p,$q,$r);
}



sub geom_DirCosines{
	my ($x,$y,$z) = @_; 

	my $d = geom_R($x,$y,$z);
	return ($x/$d,$y/$d,$z/$d);

}

sub geom_AngleBetween2Points{
	my ($x1,$y1,$z1,$x2,$y2,$z2) = @_; 

	my ($l1,$m1,$n1) = geom_DirCosines($x1,$y1,$z1);
	my ($l2,$m2,$n2) = geom_DirCosines($x2,$y2,$z2);

	return rad2deg(acos($l1*$l2 + $m1*$m2 + $n1*$n2)) ; 

}
sub geom_Distance_3D{
	my ($x1,$y1,$z1,$x2,$y2,$z2) = @_ ; 
	return geom_Distance($x1,$y1,$z1,$x2,$y2,$z2) ;
}

sub geom_Distance{
	my ($x1,$y1,$z1,$x2,$y2,$z2) = @_ ; 
	carp " kkk " if(!defined $x2);
	

	my $dx = $x1 - $x2 ; 
	my $dy = $y1 - $y2 ; 
	my $dz = $z1 - $z2 ; 

	my $dx2 = $dx*$dx ; 
	my $dy2 = $dy*$dy ; 
	my $dz2 = $dz*$dz ; 

	my $s = $dx2 + $dy2 + $dz2 ; 
	return util_format_float(sqrt($s),3) ; 
}

sub geom_TransformAxis{
	my ($x1,$y1,$z1,$x2,$y2,$z2) = @_ ; 
	

	my $dx = $x1 - $x2 ; 
	my $dy = $y1 - $y2 ; 
	my $dz = $z1 - $z2 ; 

	return ($dx,$dy,$dz);
}

sub geom_AngleBetweenThreePoints{
	my ($x1,$y1,$z1,$x2,$y2,$z2,$x3,$y3,$z3) = @_ ; 
	

	my ($d1,$d2,$d3) = geom_TransformAxis($x1,$y1,$z1,$x2,$y2,$z2) ;
	my ($e1,$e2,$e3) = geom_TransformAxis($x3,$y3,$z3,$x2,$y2,$z2) ;

	return geom_AngleBetween2Points($d1,$d2,$d3,$e1,$e2,$e3);
}

#debug function
sub geom_CrossProduct{
   my ($a,$b,$str) = @_ ; 
   my $out = $a  x $b;
   #print "Crossproduct $out  = $str \n";
	
}

sub geom_DotProduct{
   my ($a,$b,$str) = @_ ; 
   my $out = $a  * $b;
   #print "dotproduct $out  = $str \n";
	
}

########  The coordinate axes have been transformed to make the first atom at (0,0,0)
########  The first thing to do is to align two atoms such that both lie on the X axis.
########  So, first rotate the X axis, keeping Z const - geom_RotateWithXAligned 
########  This gives a rotation matrix R1 - note rotation matrices are symmetric 3x3.
#######  Then rotate keeping Y const - geom_RotateWithZAligned. Very similar
########  Function naming is not good though

sub geom_AlignXaxis2Vector{
   my ($point) = @_ ; 
   my ($XXX,$newY,$ZZZ,$vec,$rotMatrix1) = geom_RotateWithXAligned($point);
   $point = vector( $vec->x(), $vec->y(), $point->z() );
   my ($newX,,$newZ,$rotMatrix2) ;
   ($newX,$newY,$newZ,$vec,$rotMatrix2) = geom_RotateWithZAligned($point,$newY);

   geom_CrossProduct($newX,$newY,"xy");
   geom_CrossProduct($newY,$newZ,"yz");
   geom_CrossProduct($newZ,$newX,"zx");

   my $R = $rotMatrix1 * $rotMatrix2 ;

   return ($vec,$R,$rotMatrix1,$rotMatrix2,$newX,$newY,$newZ) ;
}

sub geom_RotateWithXAligned{
    use Math::VectorReal qw(:all);  # Include O X Y Z axis constant vectors

    my ($orig) = @_ ; 
    my $a = vector( $orig->x(), $orig->y(), 0 );
    
    my $nx = $a->norm ;
    my $ny = $nx x Z ;
    my $R = vector_matrix( $nx, $ny, Z );   # make the rotation matrix

    my $rotatedWithXAligned = $a * (  $R ) ;
    my $newX = $rotatedWithXAligned->norm ;
    my $newY =    $newX x Z ;
    
    #print '$a * $R (vector -> vector)',"\n", $rotatedWithXAligned, "\n";
    die "Expected y to be zero " if(! geom_IsZero($rotatedWithXAligned->y()));
    die "Expected z to be zero " if(! geom_IsZero($rotatedWithXAligned->z()));
    
    carp "Expected length to be the same " if(! geom_IsLengthSame($a,$rotatedWithXAligned));
    
    return ($newX,$newY,Z,$rotatedWithXAligned,$R);
}


sub geom_RotateWithZAligned{
    use Math::VectorReal qw(:all);  # Include O X Y Z axis constant vectors

    my ($orig,$newY) = @_ ; 
    my $a = vector( $orig->x(), 0,  $orig->z() );
    
    my $nx = $a->norm ;
    my $nz = $nx x $newY;
    my $R = vector_matrix( $nx, $newY , $nz );   # make the rotation matrix

    my $rotatedWithZAligned = $a * $R ;
    my $newX = $rotatedWithZAligned->norm ;
    my $newZ = $newX x $newY;
    

    die "Expected y to be zero " if(! geom_IsZero($rotatedWithZAligned->y()));
    die "Expected z to be zero " if(! geom_IsZero($rotatedWithZAligned->z()));
    die "Expected length to be the same " if(! geom_IsLengthSame($a,$rotatedWithZAligned));
    
    
    return ($newX,$newY,$newZ,$rotatedWithZAligned,$R);
}

sub geom_IsZero{
	my ($num) = @_ ; 
	if(abs($num) < $EPSILON){
		return 1 ; 
	}
	else{
		return 0 ; 
	}

}

sub geom_IsLengthSame{
	my ($a,$b) = @_ ; 
    my $len1 =  $a->length;
    my $len2 =  $b->length;
	my $lendiff = $len1 - $len2 ; 

	#print "$len1 $len2 \n";
    return ( geom_IsZero($lendiff));
}

sub geom_Makepoint{
	my (@p) = @_ ; 
	$, = " " ;
	#print @p , "\n";

	return \@p ;
}

sub geom_PointsOnACircle{
    my ($radius,$x,$y,$howmany) = @_ ;
	$howmany = $howmany /4 ;
	my $radius2 = $radius * $radius ;
	my @p ;
	my $delta = ($radius/($howmany));
	#print "delta = $delta \n";

	push @p, geom_Makepoint($radius,0); 
	push @p, geom_Makepoint(-$radius,0); 
	push @p, geom_Makepoint(0,$radius); 
	push @p, geom_Makepoint(0,-$radius); 

	my $initial = $radius -$delta ;
    while($initial > 0 ){
		#print "$initial = init \n";
		my $y = sqrt($radius2 - $initial*$initial);
	    push @p, geom_Makepoint($initial,$y); 
	    push @p, geom_Makepoint($initial,-$y); 
	    push @p, geom_Makepoint(-$initial,$y); 
	    push @p, geom_Makepoint(-$initial,-$y); 
		$initial = $initial - $delta ;
	}
    my $len = @p ;
	#print STDERR "There were $len points \n";
	return \@p ;
}

sub geom_Closest2{
    my ($x,$y,$z,$atoms) = @_;
    my ($mindist,$atom); 
    $mindist = 10000 ;
        foreach my $atom1 (@{$atoms}){
		    my ($p,$q,$r) = $atom1->Coords();
			if(!defined $p){
				$atom1->Print();
			}

            my $d = geom_Distance($x,$y,$z,$p,$q,$r);
            if($d < $mindist){
			     $mindist = $d ;
				 $atom = $atom1 ;
             }
        }
    return ($mindist,$atom);
}

sub geom_RotateWithXConstantAndZwillbeZero{

    my ($orig) = @_ ; 
    my $a = vector(0, $orig->y(), $orig->z());
    
    my $ny = $a->norm ;
    my $nz = $ny x X ;


    my $R = vector_matrix( X , $ny, $nz );   # make the rotation matrix


    my $rotatedWithXAligned = $a * (  $R ) ;

    my $newY = $rotatedWithXAligned->norm ;
    my $newZ =    $newY x X ;
    
    die "Expected z to be zero " if(! geom_IsZero($rotatedWithXAligned->z()));
    carp "Expected length to be the same " if(! geom_IsLengthSame($a,$rotatedWithXAligned));
    
    return (X,$newY,$newZ,$rotatedWithXAligned,$R);
}

sub geom_Align3PointsToXYPlane{
	my ($pdb1,$atoms1,$verbose,$flip) = @_ ; 
	my @atoms = @{$atoms1};
	my $natoms = @atoms ; 
	#die "Need only 3 atoms got $natoms" if(@atoms != 3);
    my $a0  = shift @atoms or die ;
    my $a1  = shift @atoms  or die;
    my $a2  = shift @atoms  or die;
	my @done ;
	push @done, $a0 ;
	push @done, $a1 ;
	push @done, $a2 ;
    
	if($verbose){
	print STDERR "before\n";
        $a0->Print();
        $a1->Print();
        $a2->Print();
	}

    my ($newX,$newY,$newZ) = $pdb1->AlignXto2Atoms($a0,$a1);

    my $p2 = vector( $a2->x(), $a2->y(), $a2->z() );
	if(defined $flip){
	    print "FLIPPING\n";
	    $p2 = -$p2 ;
	}

    my ($XXX,$newY1,$ZZZ,$vec,$R) = geom_RotateWithXConstantAndZwillbeZero($p2);
	my @allatoms = $pdb1->GetAtoms();
    $pdb1->ApplyRotationMatrix($R,\@allatoms);

	if($verbose){
	print STDERR "after \n";
        $a0->Print();
        $a1->Print();
        $a2->Print();
	}
	return (\@done,\@atoms) ;
}


sub geom_Distance_2D{
	my ($x1,$y1,$x2,$y2) = @_ ; 
	my $dx = $x1 - $x2 ; 
	my $dy = $y1 - $y2 ; 

	my $dx2 = $dx*$dx ; 
	my $dy2 = $dy*$dy ; 

	my $s = $dx2 + $dy2 ;
	return util_format_float(sqrt($s),3) ; 
}

sub geom_GetPointWithRequiredAngle{
    my ($midX,$midY,$maxX,$maxY,$reqangle,$allpoints,$thesh) = @_ ;
	die if(!defined $thesh);
	my @l = @{$allpoints} ;
	my $min = 1000 ; 
	my $mx; 
	my $my; 
	while(@l){
		my $x = shift @l;
		my $y = shift @l;
	    my $angle = geom_AngleBetweenThreePoints($x,$y,0,$midX,$midY,0,$maxX,$maxY,0);
		my $diff = abs ($angle -$reqangle);
		if($diff < $min){
			$min = $diff ; 
			$mx = $x ;
			$my = $y ;
		}
	}
	die "Expect min = $min less than $thesh " if($min > $thesh);
	return ($mx,$my) ;
}


sub geom_GetABForEllipse{
    my ($midX,$midY,$allpoints) = @_ ;
	my ($max,$min);
	$max = 0;
	$min = 10000000;
	my ($maxX,$maxY);
	my ($minX,$minY);
	my @l = @{$allpoints} ;
	while(@l){
		my $x = shift @l;
		my $y = shift @l;
		my $dist = geom_Distance_2D($midX,$midY,$x,$y);
		if($dist > $max){
		    $max = $dist ;
			$maxX = $x ; 
			$maxY = $y ; 
		}
	}

	# adjust midpoint
    my ($x,$y) = geom_GetPointWithRequiredAngle($midX,$midY,$maxX,$maxY,180,$allpoints,3);
	die "Expected a point with 180" if(!defined $x);
	$midX = ($x + $maxX)/2;
	$midY = ($y + $maxY)/2;

    ($x,$y) = geom_GetPointWithRequiredAngle($midX,$midY,$maxX,$maxY,90,$allpoints,3);
	die "Expected a point with 90" if(!defined $x);
	my $dist = geom_Distance_2D($midX,$midY,$x,$y);
	$min = $dist ;
	$minX = $x ;
	$minY = $y ;
	return ($max,$min,$maxX,$maxY,$minX,$minY,$midX,$midY);
}


sub AddVectorsAboutaCentre{
	my ($centre,$otherpoints,$potdiffs) = @_ ;
	my @otherpoints = @{$otherpoints};
	my @potdiffs = @{$potdiffs};
	my $firspoint = shift @otherpoints;
	my $fmagnitude = shift @potdiffs;
	my $vector = MakeVectorFrom2Points($centre,$firspoint);
	while (@otherpoints){
		my $p = shift @otherpoints ;
		my $magnitude = shift @potdiffs ;
	    my $v = MakeVectorFrom2Points($centre,$p);
		$vector = Add2Vectors($vector,$v,$fmagnitude,$magnitude);
		$fmagnitude = GetMagnitude($vector);
	}
	return $vector;
}

sub MakeVectorFrom2Points_Atoms{
	my ($a1,$a2) = @_;
	my (@p1) = $a1->Coords();
	my (@p2) = $a2->Coords();
	return MakeVectorFrom2Points(\@p1,\@p2);
}
sub MakeVectorFrom2Points{
	my ($p1,$p2) = @_;
	my ($x,$y,$z) = @{$p1};
	my ($p,$q,$r) = @{$p2};
    my $P1 = vector( $x,$y,$z);
	my $P2 = vector( $p,$q,$r);
	my $P = $P1 - $P2 ;
	return $P ;
}

sub geom_GetPointsBetween{
	my ($x1,$y1,$x2,$y2) = @_ ; 
	my @ret ; 
	return @ret if(abs($x1 -$x2) <=1 && abs($y1 -$y2) <=1 );
	my $midx = ($x1 + $x2)/2 ; 
	my $midy = ($y1 + $y2)/2 ; 
	push @ret, $midx ;
	push @ret, $midy ;

	my @ret1 = geom_GetPointsBetween($x1,$y1,$midx,$midy);
	push @ret1, $x1 ;
	push @ret1, $y1 ;
	my @ret2 = geom_GetPointsBetween($midx,$midy,$x2,$y2);
	push @ret2, $x2 ;
	push @ret2, $y2 ;
	
	push @ret, @ret1 ;
	push @ret, @ret2 ;
	return @ret ;
}

sub geom_ExtendTwoPointsDouble{
	my ($x1,$y1,$midX,$midY) = @_ ; 

	my $x2 = 2*$midX - $x1;
	my $y2 = 2*$midY - $y1;
	return ($x2,$y2);
}

sub geom_GetCircleAroundPoint{
	my ($x1,$y1,$radius) = @_ ; 

	my $maxX = $x1 + $radius;
	my $maxY = $y1 + $radius;

	my $minX = $x1 - $radius;
	my $minY = $y1 - $radius;

	my $midX = ($maxX + $minX)/2;

	my $x = $minX ;
	my @points ; 
	while($x < $midX){
		foreach my $y ($minY...$maxY){
	         my $dist = geom_Distance_2D($x1,$y1,$x,$y);
			 if(abs($dist - $radius) < 0.1){
			   push @points, $x ;
			   push @points, $y ;
			   last ;
			 }
		}
	    $x++ ; 
	}
	$x = $maxX ;
	while($x > $midX){
		foreach my $y ($minY...$maxY){
	         my $dist = geom_Distance_2D($x1,$y1,$x,$y);
			 if(abs($dist - $radius) < 0.1){
			   push @points, $x ;
			   push @points, $y ;
			   last ;
			 }
		}
	    $x-- ; 
	}
	return @points ;
}

sub geom_GetPointsTopDownNorthSouthEastWest{
	    my ($k,$distsurface) = @_ ;
	    my ($x,$y,$z) = geom_MakeCoordFromKey($k);
		  my @l ; 
	      my $newkey = geom_MakeKeyFromCoord($x-1,$y,$z);
		  push @l, $newkey ;
	      $newkey = geom_MakeKeyFromCoord($x+1,$y,$z);
		  push @l, $newkey ;
	      $newkey = geom_MakeKeyFromCoord($x,$y-1,$z);
		  push @l, $newkey ;
	      $newkey = geom_MakeKeyFromCoord($x,$y+1,$z);
		  push @l, $newkey ;
	      $newkey = geom_MakeKeyFromCoord($x,$y,$z-1);
		  push @l, $newkey ;
	      $newkey = geom_MakeKeyFromCoord($x,$y,$z+1);
		  push @l, $newkey ;
		  return @l ;
}

sub geom_GetPointsAroundOnePoint{
	    my ($k,$distsurface) = @_ ;
	    my ($x,$y,$z) = geom_MakeCoordFromKey($k);
		my @l ; 
	    my $Sx  = $x- $distsurface ; my $Ex  = $x + $distsurface ; 
		my $Sy  = $y- $distsurface ; my $Ey  = $y + $distsurface ; 
		my $Sz  = $z- $distsurface ; my $Ez  = $z + $distsurface ; 
        foreach my $p ($Sx..$Ex){
            foreach my $q ($Sy..$Ey){
                foreach my $r ($Sz..$Ez){
	                  push @l, geom_MakeKeyFromCoord($p,$q,$r);
			    }
			}
		}
		return @l ;
}

sub geom_MakeKeyFromCoord{
	my ($X,$Y,$Z) = @_ ;
	$Z = "0" if(!defined $Z);

	## create new grp 
	my $str = "$X#$Y#$Z";
	return $str ;
}
sub geom_MakeCoordFromKey{
	my ($str) = @_ ;
	my ($x,$y,$z) = split "#", $str ;
	return ($x,$y,$z);
}


sub geom_AngleBetweenTwoVectors_4Atoms{
	my ($a1,$a2,$b1,$b2) = @_ ;
	my $v1 = MakeVectorFrom2Points_Atoms($a1,$a2);
	my $v2 = MakeVectorFrom2Points_Atoms($b1,$b2);
	return geom_AngleBetweenTwoVectors($v1,$v2);
}

sub geom_AngleBetweenTwoVectors{
	my ($V1,$V2) = @_ ;
	my $a = $V1->norm;
	my $b = $V2->norm;
	my $angle = util_format_float(rad2deg(acos($a . $b)),1)  ;
	return $angle ;
}


sub geom_GetVectorAtAngle{
     my ($angledegree) = @_ ;
     my $angle = util_format_float(((3.14/180)*$angledegree));
	 my $y = cos($angle);
	 my $z = sin($angle);
     return ($y,$z);
}

sub geom_GetSetOfVectorsForFullCircle{
     my ($total,$numberofpieces) = @_ ;
     my $pieceangle = $total/$numberofpieces ; 
     my @l ; 
     foreach my $val (1..$numberofpieces){
           my ($y,$z) = geom_GetVectorAtAngle($val*$pieceangle);
	       push @l, $y;
	       push @l, $z;
     }
     return @l ;
}
