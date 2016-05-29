#!/usr/bin/bash

pushd data
for fn in `ls *antenna-posi*.csv`; do 
    STATION=`echo $fn|sed -e 's/\(.*\)-ant.*/\1/'`; cat $fn|egrep '^C'|sed -e "s/C/${STATION^^?}/" ;
done|sort -k 1.3,1.9|grep -v FI |sed -e's/^\(.*,.*,.*,.*\),.*,.*,.*,.*,.*/\1/'> etrs-phase-centres-tmp.csv
popd
