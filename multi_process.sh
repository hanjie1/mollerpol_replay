#!/bin/tcsh

###########################################################################################################
# Optimized Parallel Replay for ifarm
###########################################################################################################

if ( $#argv < 5 ) then
    echo "Usage: ./multi_process.sh <runnum> <first_seg> <max_seg> <first_event> <nevents>"
    exit 1
endif

set runnum = $1
set firstsegment = $2
set maxsegments = $3
set firstevent = $4
set nevents = $5

# 1. Environment Setup: Setup as per your analyser build, input and output directories
source /w/halla-scshelf2102/moller12gev/asar/Analyzer/local/analyzer-1.7.18/bin/setup.csh
setenv MOLLERPOL_HOME /w/halla-scshelf2102/moller12gev/asar/Moller_Polarimeter/PolGEM/mollerpol_analyzer
setenv MOLLERPOL_BUILD ${MOLLERPOL_HOME}/build
setenv LD_LIBRARY_PATH ${MOLLERPOL_BUILD}:${LD_LIBRARY_PATH}

setenv MOLLER_REPLAY /w/halla-scshelf2102/moller12gev/asar/Moller_Polarimeter/Replay/mollerpol_replay
setenv DB_DIR ${MOLLER_REPLAY}/DB
setenv OUT_DIR /volatile/halla/moller12gev/asar/GEMRootfiles
setenv LOG_DIR ${OUT_DIR}/logs

mkdir -p ${LOG_DIR}

# Set the number of simultaneous cores to use (ifarm nodes usually have 32-64,
# but keep it reasonable, e.g., 8-16)
set MAX_CORES = 12

# 2. Parallel Loop
echo "Starting parallel analyzer processes (Max Cores: ${MAX_CORES})..."

set seg = $firstsegment
@ lastseg_limit = $firstsegment + $maxsegments

while ( $seg < $lastseg_limit )
    # --- Data Discovery per Segment ---
    set found_dir = ""
    foreach d ( /adaqeb1/data1 /adaqeb2/data1 /adaqeb3/data1 /cache/halla/sbs/GEp/raw /volatile/halla/moller12gev/asar/TestLab_Data )
        if ( -f "${d}/moller_ssp_${runnum}.evio.${seg}" ) then
            set found_dir = "$d"
            break
        endif
    end

    if ( "$found_dir" == "" ) then
        echo "WARNING: Segment ${seg} not found. Skipping."
        @ seg++
        continue
    endif

    setenv DATA_DIR "$found_dir"

    # --- Job Throttling ---
    # Count current background jobs; wait if we are at MAX_CORES
    set njobs = `jobs | wc -l`
    while ( $njobs >= $MAX_CORES )
        sleep 2
        set njobs = `jobs | wc -l`
    end

    # --- Launch ---
    echo "Launching Segment ${seg} from ${DATA_DIR}"
    echo "Launching file moller_ssp_${runnum}.evio.${seg}"
    analyzer -b -q "${MOLLER_REPLAY}/replay_moller_gem_multi.C(${runnum},${seg},${firstevent},${nevents})" >& ${LOG_DIR}/replay_${runnum}_seg${seg}.log &

    @ seg++
end

echo "Waiting for remaining segments to finish..."
wait

# 3. Merging
set O_FNAME = "${OUT_DIR}/moller_uvagem_replayed_${runnum}_merged.root"

# We look specifically for the files created by THIS run number
set ofiles = `ls ${OUT_DIR}/moller_polgem_run${runnum}_seg*.root`

if ( "$ofiles" != "" ) then
    echo "Merging `echo $ofiles | wc -w` files into $O_FNAME"

    # Run hadd
    hadd -k -f $O_FNAME $ofiles

    # --- SAFETY CHECK ---
    # Only delete if hadd returned a success code (0)
    # and the merged file actually exists.
    if ( $status == 0 && -f "$O_FNAME" ) then
        echo "Merge successful. Deleting individual segment files for Run ${runnum}..."
        rm -f $ofiles
    else
        echo "ERROR: hadd failed or merged file not found. Keeping segment files for safety."
    endif
else
    echo "ERROR: No root files found to merge for Run ${runnum}."
endif
rm *.log
echo "Done."
