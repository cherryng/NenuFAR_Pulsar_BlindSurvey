#!/bin/bash 

#-------------------------------
# Update Extended_grid_[date].csv in case the previously planned observations were not actually carried out
# 
# INF = (a copy) of the last good csv grid (the version before the month with issue
# DIR = the .parset_user that was actually observed in the month with issue
# 
# The script will look in the $DIR and change the status of those observed lines in the last-good csv
# A new file named Extended_grid_[date].csv_patched will be written out
#-------------------------------


export INF=$1
export OUT=$1.patched
cp $INF $OUT
export DIR=$2

for F in $DIR/*user; do
    echo "Working on $F ----------------------------------------"

    # Extract ObsID
    ObsID=$(echo "$F" | awk -F'_' '{print $2}')

    # Beam0 ID
    B0ID=$(grep "Beam\[0\].target" "$F" | awk -F'_' '{print $2}')
    OLDLINE=$(grep "^$B0ID" "$OUT")
    NEWLINE=$(echo "$OLDLINE $ObsID" | awk '{print $1,$2,$3,$4,$5,"1",$7,$9}')

    # Escape variables for sed
    ESC_OLDLINE=$(printf '%s\n' "$OLDLINE" | sed 's/[\/&]/\\&/g')
    ESC_NEWLINE=$(printf '%s\n' "$NEWLINE" | sed 's/[\/&]/\\&/g')
    echo "ESC_NEWLINE: $ESC_NEWLINE"

    if [[ -n "$ESC_OLDLINE" && -n "$ESC_NEWLINE" ]]; then
        sed -i '' "s|$ESC_OLDLINE|$ESC_NEWLINE|g" "$OUT"  # Use '' if on macOS
    else
        echo "Error: ESC_OLDLINE or ESC_NEWLINE is empty"
    fi

    # Beam1 ID
    B1ID=$(grep "Beam\[1\].target" "$F" | awk -F'_' '{print $2}')
    OLDLINE=$(grep "^$B1ID" "$OUT")
    NEWLINE=$(echo "$OLDLINE $ObsID" | awk '{print $1,$2,$3,$4,$5,"1",$7,$9}')

    # Escape variables for sed
    ESC_OLDLINE=$(printf '%s\n' "$OLDLINE" | sed 's/[\/&]/\\&/g')
    ESC_NEWLINE=$(printf '%s\n' "$NEWLINE" | sed 's/[\/&]/\\&/g')
    echo "ESC_NEWLINE: $ESC_NEWLINE"    

    if [[ -n "$ESC_OLDLINE" && -n "$ESC_NEWLINE" ]]; then
        sed -i '' "s|$ESC_OLDLINE|$ESC_NEWLINE|g" "$OUT"  # Use '' if on macOS
    else
        echo "Error: ESC_OLDLINE or ESC_NEWLINE is empty"
    fi

    # Beam2 ID
    B2ID=$(grep "Beam\[2\].target" "$F" | awk -F'_' '{print $2}')
    OLDLINE=$(grep "^$B2ID" "$OUT")
    NEWLINE=$(echo "$OLDLINE $ObsID" | awk '{print $1,$2,$3,$4,$5,"1",$7,$9}')

    # Escape variables for sed
    ESC_OLDLINE=$(printf '%s\n' "$OLDLINE" | sed 's/[\/&]/\\&/g')
    ESC_NEWLINE=$(printf '%s\n' "$NEWLINE" | sed 's/[\/&]/\\&/g')
    echo "ESC_NEWLINE: $ESC_NEWLINE"

    if [[ -n "$ESC_OLDLINE" && -n "$ESC_NEWLINE" ]]; then
        sed -i '' "s|$ESC_OLDLINE|$ESC_NEWLINE|g" "$OUT"  # Use '' if on macOS
    else
        echo "Error: ESC_OLDLINE or ESC_NEWLINE is empty"
    fi


    # Beam3 ID
    B3ID=$(grep "Beam\[3\].target" "$F" | awk -F'_' '{print $2}')
    OLDLINE=$(grep "^$B3ID" "$OUT")
    NEWLINE=$(echo "$OLDLINE $ObsID" | awk '{print $1,$2,$3,$4,$5,"1",$7,$9}')

    # Escape variables for sed
    ESC_OLDLINE=$(printf '%s\n' "$OLDLINE" | sed 's/[\/&]/\\&/g')
    ESC_NEWLINE=$(printf '%s\n' "$NEWLINE" | sed 's/[\/&]/\\&/g')
    echo "ESC_NEWLINE: $ESC_NEWLINE"    

    if [[ -n "$ESC_OLDLINE" && -n "$ESC_NEWLINE" ]]; then
        sed -i '' "s|$ESC_OLDLINE|$ESC_NEWLINE|g" "$OUT"  # Use '' if on macOS
    else
        echo "Error: ESC_OLDLINE or ESC_NEWLINE is empty"
    fi    
done
