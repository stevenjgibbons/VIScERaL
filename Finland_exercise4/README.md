
In this example we solve for both the location vectors and the slowness vectors
simultaneously. The files new_rellocs_init.txt and new_slovecs_init.txt contain a
set of starting estimates for the locations and slowness vectors respectively
(these were the median location vectors from a previous iteration of this procedure:
you could equally well start with the purely theoretical slowness vectors and
then any initial set of location vectors which performs well for a good set of
fixed slowness vectors).

The script **run_all.sh** runs the script **locate_all_vs_all_variableslovecs.sh**
repeatedly (100 times in this example) perturbing the starting values of locations
and slowness vectors randomly each time within the specified limits.
