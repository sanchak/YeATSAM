
sub util_wait_on_lockfiles{
    my ($ofh,$lockfiles) = @_ ;
    my @lockfiles = @{$lockfiles};
    print $ofh "echo Waiting for lock files to disappear\n";
    print $ofh "while(1) \n";
    print $ofh "sleep 5\n";
    my $str = "";
    map { $str.= " -e $_ || " ; } (@lockfiles);
    $str .= " -e $lockfiles[0]" ;
    print $ofh "if( $str)  then \n";
    print $ofh "else \n";
    print $ofh "break \n";
    print $ofh "endif \n";
    print $ofh "end \n";
}

sub util_printHtmlHeader{
    my ($ofh,$header1,$header2) = @_ ;
    print $ofh "<html> \n";
    print $ofh "<h1>$header1</h1> \n";
    print $ofh "<body> \n";

    if(defined $header2){
    print $ofh "<html> \n";
    print $ofh "<h2>$header2</h2> \n";
    print $ofh "<body> \n";
    }
}



sub util_printHtmlEnd{
    my ($ofh) = @_ ;
    print $ofh "</body> \n";
    print $ofh "</html> \n";
}

sub util_HtmlizeLine{
    my ($line) = @_ ;
    chomp $line ; 
    #return $line  . "<br />" ;
    return $line ;

}

sub util_HtmlTableHead{
    my ($ofh,@headers) = @_ ;
    print $ofh "<table border=\"1\" cellpadding=\"5\" cellspacing=\"5\" width=\"100%\">\n";
    print $ofh "<tr>\n";
    foreach my $h (@headers){
        print $ofh "<th>$h</th>";
    }
    print $ofh "\n";
    print $ofh "</tr>\n";

}


sub util_HtmlTableEnd{
    my ($ofh) = @_ ;
    print $ofh "</table>\n";
}



sub util_HtmlTableCell{
    my ($str) = @_ ;
    return "<td>$str</td>";
}





sub util_MakeLink{
    my ($nm,$link) = @_ ;
    return  "<a href=\"$link\"> $nm</a>";
}

sub util_EC_CorrelatePDBS{
    my ($info,$a,$b,$VALUE) = @_ ; 
    print "util_EC_CorrelatePDBS for $a $b\n";
    my $x = util_getECfromPDB($info,$a);
    my $y = util_getECfromPDB($info,$b);
    if(!defined $x || !defined $y){
        print "Could not get EC for $a \n" if(!defined $x);
        print "Could not get EC for $b \n" if(!defined $y);;
        return undef;
    }

    my @ec1 = @{$x};
    my @ec2 = @{$y};
    my $ec1 = $ec1[0];
    my $ec2 = $ec2[0];
    if(!defined $ec1 || !defined $ec2){
        print "Could not get EC for $a \n" if(!defined $ec1);
        print "Could not get EC for $b \n" if(!defined $ec2);;
        return undef;
    }
    print "util_EC_CorrelatePDBS for $ec1 $ec2\n";

    $ec1 =~ s/\./YYY/g;
    $ec2 =~ s/\./YYY/g;
    my @l1 = split "YYY", $ec1 ;
    my @l2 = split "YYY", $ec2 ;
   
   my $N1 = @l1 -1 ;
   my $N2 = @l2 -1 ;
   my $score = abs($N1 - $N2);

   foreach my $n (0..$N1){
            my $v1 = $l1[$n] ; 
            my $v2 = $l2[$n] ; 
            if($v1 != $v2){
                $score = $score + $VALUE/($n+1) ; 
                print " Found mismatch $n $score $v1 $v2\n";
                last ;
            }
   }
   return $score ; 
}

sub util_EC_AddPDB{
    my ($info,$pdb) = @_ ; 
    my $ec = util_getECfromPDB($info,$pdb);
    return $ec if(!defined $ec);
    my @l = split "\.", $ec ;
    
   my $N = @l -1 ;
   my $obj = {} ;
   foreach my $n (0..$N){
            my $v = $l[$n] ; 
            $obj->{$n}->{$v} = [] if(!exists $obj->{$n}->{$v}) ;
            push @{$obj->{$n}->{$v}}, $pdb ; 
   }
}


sub util_EC_CreateLevels{
    my ($ecdone) = @_ ; 
    my $MAXLEVEL = 0 ;
    my $obj = {} ; 
    foreach my $ec (sort keys %{$ecdone}){
        my @l = split "\.", $ec ;
        $MAXLEVEL = @l if(@l > $MAXLEVEL);
    
        my $N = @l -1 ;
        foreach my $n (0..$N){
            my $obj->{$n} = [];
        }
    }
    return ($obj,$MAXLEVEL) ;
}
