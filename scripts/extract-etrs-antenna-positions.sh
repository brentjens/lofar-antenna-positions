#!/usr/bin/bash

pushd data
for fn in `ls *antenna-posi*.csv`; do      STATION=`echo $fn|sed -e 's/\(.*\)-ant.*/\1/'`; cat $fn|egrep -v '^C'|sed -e "s/^/${STATION^^?},/" ; done|grep -v FI |sed -e's/^\(.*,.*,.*,.*\),.*,.*,.*\(,.*,.*\)/\1\2/'|sed -e's/,L/,LBA,/'|sed -e's/,H/,HBA,/' |grep -v NAME |grep BA > ../etrs-antenna-positions-tmp.csv
popd
