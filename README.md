# mollerpol_replay

# Set up the directory
1. Set up the correct enviroment variables, such as:
```
export LD_LIBRARY_PATH=/home/hamoller/hanjie/moller_analyzer/build:$LD_LIBRARY_PATH
```
The above is adding the moller customized library "libMollerPol.so" to the path.
```
export DB_DIR=/home/hamoller/hanjie/moller_analyzer/replay/DB
```
This points to the replay data base.

2. Set up the data file directory. The default data searching path is data/
3. Set up the Rootfiles directory. The root files shouldn't be saved under the /home directory. A better place is /adaq1/data1/hamoller/<username>

# Run the replay
```
analyzer "replay.C(run_number[,number of events])"
```
# To analyze GEM data

You must set up the environment variables to point to the correct libraries and database.

Example:

```bash
export MOLLERPOL_HOME=/path/to/your/mollerpol_analyzer
export MOLLERPOL_BUILD=$MOLLERPOL_HOME/build
export LD_LIBRARY_PATH=$MOLLERPOL_BUILD:$LD_LIBRARY_PATH
export MOLLER_REPLAY=/path/to/mollerpol_replay
export DB_DIR=$MOLLER_REPLAY/DB
```

---

# Parallel GEM Replay

For GEM data containing multiple EVIO segments, the replay can be run using:

```bash
./multi_process.sh run_number first_segment number_of_segments first_event number_of_events
```

For example:

```bash
./multi_process.sh 2315 0 10 0 -1
```

This processes run `2315`, starting from segment `0`, for 10 segments, and processes all events in each segment.

The arguments are:

```text
run_number
first_segment
number_of_segments
first_event
number_of_events
```

Use `-1` for `number_of_events` to process all events.

The script runs multiple Analyzer processes in parallel and searches for raw data files with names of the form:

```text
moller_ssp_<run_number>.evio.<segment>
```

The number of simultaneous Analyzer jobs is controlled inside `multi_process.sh` by:

```tcsh
set MAX_CORES = 12
```

Each segment is processed using:

```text
replay_moller_gem_multi.C
```

After all segments are processed, the script merges the individual ROOT files using `hadd`.

The merged ROOT file has the form:

```text
moller_uvagem_replayed_<run_number>_merged.root
```

Individual segment ROOT files are removed only after a successful merge.

Replay logs are stored under the output directory in:

```text
logs/
```

---

# Alignment Script

The script `GEM_align.C`, located in:

```text
scripts/gems/
```

is used to align the GEM modules. It recalculates module positions from reconstructed tracks by iterating until the chi-square is minimized.

To run the alignment:

```bash
root -l -q 'scripts/gems/GEM_align.C("input_alignment.txt")'
```

The script generates output files containing optimized geometry parameters.

Update the GEM database module positions using the resulting alignment parameters to improve tracking resolution.
