#!/bin/csh  -f

if($#argv != 3  ) then 
    echo "<name> <image  or pdf> <filetype>"
	exit 
endif 

set nm=$1 
set file=$2
set filetype=$3

cp fig.sample.tex fig.$nm.tex 
replacestring.pl -in fig.$nm.tex -wh figlabel -with "fig$nm" -same -outf kkkkk
replacestring.pl -in fig.$nm.tex -wh microbio -with $nm -same -outf kkkkk
replacestring.pl -in fig.$nm.tex -wh TYPE -with "$filetype" -same -outf kkkkk

echo "\input{fig.$nm.tex}" >> fig.tex


\cp -f $2 $nm.$filetype 

 

