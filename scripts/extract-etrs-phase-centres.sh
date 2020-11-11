#!/usr/bin/bash

pushd data

# Remove the temp file if it already exists
rm -f ../share/lofarantpos/etrs-phase-centres-tmp.csv

echo "STATION,FIELD,ETRS-X,ETRS-Y,ETRS-Z" > ../share/lofarantpos/etrs-phase-centres-tmp.csv
for fn in `ls *antenna-posi*.csv`; do 
    STATION=`echo $fn|sed -e 's/\(.*\)-ant.*/\1/'`; 
    cat $fn | egrep '^C' \
    		| sed -e "s/C/${STATION^^?},/" ;
done	| tr -d ' ' \
		| sort -k 1.3,1.9 \
		| grep -v FI \
		| sed -e's/^\(.*,.*,.*,.*\),.*,.*,.*,.*,.*/\1/' >> ../share/lofarantpos/etrs-phase-centres-tmp.csv

popd
