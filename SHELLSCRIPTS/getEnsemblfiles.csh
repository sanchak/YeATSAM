#!/bin/csh  -f

if($#argv != 1  ) then 
    echo "Wrong args, required 1 " ; exit 
endif 

set list = $1 

foreach i (`cat $list`)
	wget -r --no-parent --reject "index.html*"  ftp://ftp.ensemblgenomes.org/pub/plants/current/fasta/$i/pep/
end 

 

