
This directory contains the following files:  

(1) **locate_all_vs_all_variableslovecs_DPRK.sh**  

A script which takes the input locations and slowness vectors, randomizes them within the specified limits, and iterates towards
a new solution which fits the set of differential times better - moving both event locations and solution vectors. The script
takes an identifier runID which makes the output files unique

(2) **run_all.sh**  

A (trivial) script which just calls locate_all_vs_all_variableslovecs_DPRK.sh repeatedly for a set of different runID.
The outcome of the runs will vary due to the randomization of the intial values.  

(3) **new_rellocs_init.txt**  

A set of starting location estimates. All flags are set to "S" (solve) except for DPRK2 which
is set to "F" (fixed).  

(4) **new_slovecs_init.txt**  

A set of starting slowness vectors.  

