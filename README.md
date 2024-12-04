# Scheduling observations

- To plan upcoming observations:

`python SelectPointings.py -f Dates-to-observe.dat -g Extended_grid_[lastdate].csv -write -flag` 

- To restore the flags of previously scheduled observations that somehow were not executed:

`Patch_Schedules_In_GridFile.sh Extended_grid_[lastgood].csv [DIR-of-executed-parset_user]`
