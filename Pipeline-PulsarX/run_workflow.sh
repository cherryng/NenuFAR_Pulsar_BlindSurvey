#!/bin/bash
set -e  # stop script if any command fails

# Check argument
if [ $# -lt 1 ]; then
    echo "Usage: $0 <list_file>"
    exit 1
fi

#[To be edited] The path where all the scripts and depedencies have been downloaded to, links in this doc:
#https://docs.google.com/document/d/1fg8rR74tqlfrpkYrHIXKpuB_H2Uf3CI5rhhII6GiWLg/edit?usp=sharing
SCRIPTPATH=""


#This is a text file that continue the full paths of all the .fbk files to be processed
LISTFILE="$1"

# Check file exists
if [ ! -f "$LISTFILE" ]; then
    echo "Error: '$LISTFILE' not found"
    exit 1
fi

# Loop through each line
while IFS= read -r filename; do
    [[ -z "$filename" ]] && continue      # skip empty lines
    [[ "$filename" =~ ^# ]] && continue   # skip comment lines

    echo "Processing: $filename"
    INF=$(echo "$filename" | awk -F'/' '{print $NF}')
    ROOT=$(echo "$INF" | awk -F'.spectra' '{print $1}')
    FIL=$ROOT".flat.fil"
    echo file root $ROOT, fil file $FIL
    mkdir -p $ROOT
    cd $ROOT/
    echo Working dir is $PWD
    
    #Collect required files
    chmod 777 $PWD
    cp $SCRIPTPATH/ddplan.txt .
    cp $SCRIPTPATH/fold.template .
    cp $SCRIPTPATH/ACCEL_sift_pulsarx.py .
    cp $SCRIPTPATH/FlatTime.tar .
    chmod a+x ACCEL_sift_pulsarx.py
    cp $SCRIPTPATH/presto5.sif .
    cp $SCRIPTPATH/pulsarx.sif .
    
    #Step 1 ----------------------
    if [ -e $FIL ]; then
	echo  $FIL exists, skipping step 1
    else
	echo "1) Running Flattening ------------------"
	cp $filename .
	tar -xvf FlatTime.tar
	#docker run --rm -v "$PWD":/work -w /work alex88ridolfi/presto5 python FlatTimeSeries.py  -c 25 -t 3 -o $ROOT $INF
	apptainer exec --bind "$PWD":/work --pwd /work  presto5.sif  python FlatTimeSeries.py -c 25 -t 3 -o $ROOT $INF
	mv $ROOT.flat.fbk $FIL
    fi    

    #Step 2 --------------------
    if [ -e dedisp.done ]; then
	echo "Already dedispersed"
    else
	echo "2) Running dedisperse_all_fil -----------------"
	#docker run --rm -v "$PWD":/work -w /work ypmen/pulsarx bash -c "pwd && dedisperse_all_fil -v --ddplan  ddplan.txt -z zdot --format presto -f $FIL"
	ulimit -n 10000
	apptainer exec --bind "$PWD":/work --pwd /work pulsarx.sif bash -c "pwd && dedisperse_all_fil -v --ddplan  ddplan.txt -z zdot --format presto -f $FIL"
	touch dedisp.done
    fi

    #Step 3-----------------
    echo "3) Running realfft, searching and sifting ------"
    if [ -e sift.done ]; then
	echo "Already sifted"
    else
	#docker run --rm -v "$PWD":/work -w /work alex88ridolfi/presto5 bash -c 'realfft *.dat && find . -name "*.fft" -print0 | xargs -0 -n 1 -P 16 sh -c '\''accelsearch -zmax 0 -numharm 32 "$1" || echo "FAILED: $1"'\'' _ && python ACCEL_sift_pulsarx.py'
	apptainer exec --bind "$PWD":/work --pwd /work presto5.sif bash -c 'realfft *.dat && find . -name "*.fft" -print0 | xargs -0 -n 1 -P 16 sh -c '\''accelsearch -zmax 0 -numharm 32 "$1" || echo "FAILED: $1"'\'' _ && python ACCEL_sift_pulsarx.py'
	touch sift.done
    fi

    #Step 4 --------------
    if [ -e cands.txt ]; then
	echo "4) folding ---------------------------------"
	#docker run --rm -v "$PWD":/work -w /work ypmen/pulsarx psrfold_fil2 -v --template fold.template --candfile cands.txt --plotx -n 64 -b 64 --clfd 2 -t 16 --noarch -z zdot -f $FIL
	ulimit -n 10000
	apptainer exec --bind "$PWD":/work --pwd /work pulsarx.sif psrfold_fil2 -v --template fold.template --candfile cands.txt --plotx -n 64 -b 64 --clfd 2 -t 16 --noarch -z zdot -f $FIL
    else
	touch NoCands.txt
    fi
    
    #Tidy up--------------
    rm -f *.dat
    rm -f *.fft
    rm -f *.inf
    rm -f *.spectra.fbk fold.template ddplan.txt FlatTimeSeries.py ACCEL_sift_pulsarx.py FlatTime.tar *.py
    rm -rf *pycache*
    mkdir CANDS/
    mv *_ACCEL_0* CANDS/

    #Select only plots with SNR>10 to tar
    candfile=`ls J*.cands`
    candroot=`echo $candfile | awk -F'.cands' '{print $1"_"}'`
    awk -v root="$candroot" '!/^#/ { if ($16>10.0) printf "%s%05d.png\n", root, $1}' $candfile | tar -cvf PNG-selected-$ROOT.tar -T -

    tar -cvf PNG-$ROOT.tar *.png
    rm -f *.png
    rm -f *.fil
    touch All.done
    cd -
    
done < "$LISTFILE"
