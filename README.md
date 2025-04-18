# Scheduling observations

- To plan upcoming observations:

`python SelectPointings.py -f Dates-to-observe.dat -g Extended_grid_[lastdate].csv -write -flag` 

- To restore the flags of previously scheduled observations that somehow were not executed:

`Patch_Schedules_In_GridFile.sh Extended_grid_[lastgood].csv [DIR-of-executed-parset_user]`

# Search pipeline
It is based on the public [presto](https://github.com/scottransom/presto) suite, with specific elements adapted to NenuFAR configurations.
![BS flowchart](BS-flowchart.png?raw=true)

# Viewing candidates
Extract results by a month at a time (mandatory input flag `-month YYYYMM`), filter by significance range (using optional input flag `-s X Y`), by period range (using optional input flag `-p X Y`), by DM range (using optional input flag `-d X Y`) and by maximum harmonic association (by using optional flag `-m X`). A cross-matching with PSRCAT known pulsars is conducted to reduce the number of candidates. Finally, two tarballs `Plots_CAND_YYYYMM.tar.gz` and  `Plots_PSR_YYYYMM.tar.gz` will be generated with presto plots of events that are indentified as candidates vs known pulsars.

`python CollectPlotsForViewing.py -month 20YYMM`
