#!/bin/csh -f

if($#argv != 3) then 
  echo "Usage : <id> <num> <tag> "
  exit 
endif 

set i = $1
set num = $2
set tag = $3


ann2simpleinput.pl -out ANNOTATE.$num/$i.in -in ANNOTATE.$num/$i.outconf.annotated

$SRC/MISC/createPremoninput.pl -out $i.$num.$tag.premon.in -con $CONFIGGRP -li ANNOTATE.$num/$i.in -pr $i

