#!/bin/bash

# Download files if they do not exist
STATIONS=( CS001 CS002 CS003 CS004 CS005 CS006 CS007 CS011 CS013 CS017 CS021 CS024 CS026 CS028 CS030 CS031 CS032 CS101 CS103 CS201 CS301 CS302 CS401 CS501 RS106 RS205 RS208 RS210 RS305 RS306 RS307 RS310 RS406 RS407 RS409 RS503 RS508 RS509 DE601 DE602 DE603 DE604 DE605 FR606 SE607 UK608 DE609 PL610 PL611 PL612 IE613 LV614 FI901 )

# Chose a branch of the LOFAR repo to scrape
branch="master"

# Get antenna CSV information
for stn in "${STATIONS[@]}"; do
	echo "Downloading antenna CSV for $stn"
	lowerstn=$(echo ${stn} | awk '{print tolower($0)}')
	curl --silent "https://git.astron.nl/ro/lofar/-/raw/${branch}/MAC/Deployment/data/Coordinates/ETRF_FILES/${stn}/${lowerstn}-antenna-positions-etrs.csv" --output "${lowerstn}-antenna-positions-etrs.csv";
done


# Get the StationInfo.dat file
if [ ! -f "StationInfo.dat" ]; then
	curl --silent "https://git.astron.nl/ro/lofar/-/raw/${branch}/MAC/Deployment/data/StaticMetaData/StationInfo.dat" --output "StationInfo.dat"
fi


# Run the extraction scripts
rm -f ../share/lofarantpos/stationinfo-tmp.csv
./extract-cabinet-positions.py

./extract-etrs-antenna-positions.sh
./extract-etrs-phase-centres.sh


# Cleanup working directory if NO_CLEANUP is not defined
if [ -z ${NO_CLEANUP+x} ]; then
	rm -f "StationInfo.dat"
	for stn in "${STATIONS[@]}"; do
		lowerstn=$(echo ${stn} | awk '{print tolower($0)}')
		rm ${lowerstn}-antenna-positions-etrs.csv
	done
fi