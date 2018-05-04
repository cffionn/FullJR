To run:

source setFullJREnv.sh

This will setup some necessary environmental variables for each directory to build. Note that directories, while interdependent in code, are built independently. So if one wants only to run or modify a specific subset of the code, this is done by working within the desired directory.

If source setFullJREnv.sh gives empty paths (likely if not working on CERN computing resources), fill in the appropriate path and source.

To build a directory, cd into directory and run 'make'
To clean, cd into directory and run 'make clean'
