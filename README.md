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
