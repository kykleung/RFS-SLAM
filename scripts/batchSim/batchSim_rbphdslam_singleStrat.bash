#!/bin/bash

LANG=en_US # so that seq doesn't use commas as decimal separator
LC_NUMERIC=en_US.UTF-8 # so that seq doesn't use commas as decimal separator

cd ~/Projects/phdFilter
mkdir -p data/batch/rbphdslam_singleStrat

for Pd in $(seq 1.0 -0.1 0.1)
do
    if (( $(bc <<< "$Pd == 1.0") == 1 ))
    then
	Pd=0.99  
    fi

    for c in 0.0001 0.001 0.01 0.1 0.28 0.46 0.64 0.82 1
    do
	
	for trialseed in 1 2 3 4 5 6 7 8 9 10
	do

	    echo "Pd = $Pd   c = $c   seed = $trialseed"

	    # Edit xml config files
	    sed -e "s/<probDetection>.*<\/probDetection>/<probDetection>$Pd<\/probDetection>/" -e "s/<clutterIntensity>.*<\/clutterIntensity>/<clutterIntensity>$c<\/clutterIntensity>/" -e "s/<logDirPrefix>.*<\/logDirPrefix>/<logDirPrefix>data\/batch\/rbphdslam_singleStrat\/<\/logDirPrefix>/" -e "s/<nEvalPt>.*<\/nEvalPt>/<nEvalPt>1<\/nEvalPt>/" cfg/rbphdslam2dSim.xml > data/batch/rbphdslam_singleStrat/rbphdslam2dSim_singleStrat.xml

	    # Run the simulator
	    build/bin/rbphdslam2dSim 28 data/batch/rbphdslam_singleStrat/rbphdslam2dSim_singleStrat.xml

	    # Analyze results
	    build/bin/analysis2dSim data/batch/rbphdslam_singleStrat/

	    # Get the position and map errors at the very end
	    posError=$(tail -n 1 data/batch/rbphdslam_singleStrat/poseEstError.dat | tr -s ' ' | cut -d ' ' -f5)
	    mapError=$(tail -n 1 data/batch/rbphdslam_singleStrat/landmarkEstError.dat | tr -s ' ' | cut -d ' ' -f4)

	    # Write results to file
	    echo "$Pd   $c   $posError   $mapError" >> data/batch/batch_results_rbphdslam_singleStrat.dat

	done

    done

done
