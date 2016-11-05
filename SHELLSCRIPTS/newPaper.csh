#!/bin/csh -f

if($#argv != 1  ) then 
  echo "Usage <name of paper>"
  exit 
endif 

set what=$1

if(! -e /home/sandeepc/Bio/Data/Paper/$what/) then 
   cp -r /home/sandeepc/Bio/Data/Paper/tmpl /home/sandeepc/Bio/Data/Paper/$what/ 
   echo "alias $what 'cd /home/sandeepc/Bio/Data/Paper/$what/ ; ./do.sh  '"
   echo "alias $what 'cd /home/sandeepc/Bio/Data/Paper/$what/ ; ./do.sh  '" >> ~/.alias
   echo "$what" > ! /home/sandeepc/Bio/Data/Paper/$what/title.tex
endif 

source  ~/.alias

$1
