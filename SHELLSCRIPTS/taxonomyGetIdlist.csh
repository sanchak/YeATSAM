#!/bin/csh  -f

if($#argv != 2  ) then 
    echo "Usage: <list> <DIR> " ; exit 
endif 

set list = $1 
set DIR = $2 

foreach i (`cat $list`)
	if(! -e $DIR/$i.txt) then 
     wget "http://www.ncbi.nlm.nih.gov/taxonomy/?term=$i&report=info" 
	 \mv -f "index.html?term=$i&report=info" $DIR/$i.html
	 html2text $DIR/$i.html > ! $DIR/$i.txt
	 else
	 	echo Already done $DIR/$i.txt
	 endif
end 

 

