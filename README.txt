To run pythODE++
1) Ensure that pythODE++ is checked out under your /home/user directory. If it isn't, the scripts will not execute
2) In this same directory, build the program by executing the build script: ./build
3) Navigate back to your /home/user directory. Create a file called hostfile. The format for the file is as follows:
hostname slots=(Insert number of processes here) max_slots=(Insert number of processes here)
This file tells MPI which nodes to run on, with how many processes on each node. For example
moler slots=12 max_slots=12
would run the scripts on the moler head node with 12 processes.
4) navigate to the pythODE++/scripts directory. To run a script, execute the run-experiment.sh file, and enter the name
of the script you want to run *without* the .py extension:
./run-experiment.sh brusselator-experiment

For more information on the format of the experiment scripts, see the tutorial under the scripts directory.
