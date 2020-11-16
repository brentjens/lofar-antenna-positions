#!/usr/bin/bash

# Save working directory
pushd data


# Remove the previous version of the file
rm -f ../share/lofarantpos/etrs-antenna-positions-tmp.csv


# Parse the files to extract the antenna positions
echo "STATION,ANTENNA-TYPE,ANTENNA-ID,ETRS-X,ETRS-Y,ETRS-Z,RCU-X,RCU-Y" > ../share/lofarantpos/etrs-antenna-positions-tmp.csv
for fn in `ls *antenna-posi*.csv`; do
      STATION=`echo $fn|sed -e 's/\(.*\)-ant.*/\1/'`; 
      cat $fn|egrep -v '^C'|sed -e "s/^/${STATION^^?},/" ; 
done 	| grep -v FI |sed -e's/^\(.*,.*,.*,.*\),.*,.*,.*\(,.*,.*\)/\1\2/' \
		| sed -e's/,L/,LBA,/'|sed -e's/,H/,HBA,/'  \
		| grep -v NAME  \
		| grep BA >> ../share/lofarantpos/etrs-antenna-positions-tmp.csv

# Return to original working directory
popd
