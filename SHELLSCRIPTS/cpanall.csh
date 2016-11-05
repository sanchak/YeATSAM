# run this first:
% cpan
cpan> o conf makepl_arg INSTALL_BASE=/home/sandeepc/perl5
cpan> o conf commit

OR 

% cpan
 o conf mbuildpl_arg '--install_base /home/sandeepc/perl5'
 o conf commit

#cpan> o conf prerequisites_policy follow 

cpan -i Algorithm/Combinatorics.pm
cpan -i Math/NumberCruncher.pm
cpan -i Math/MatrixReal.pm
cpan -i Math/Geometry.pm
cpan -i Math/VectorReal.pm
cpan -i Statistics/Descriptive.pm
cpan -i Math/Geometry/Planar.pm
cpan -i Math/Combinatorics.pm
cpan -i Math/Geometry.pm 

