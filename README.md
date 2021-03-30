# Proseeker

In silico directed evolution

This is the research version of the Proseeker program written in MATLAB. It will operate in a parallel manner managed by MATLAB. Upon upload of the Python3 based distribution version an updated README will be provided.

It is set up for runs with k=100, g=20 and using the optmised PPR library for ease of replicating submitted results. An example blocking scheme is included but both blocking schemes used in the paper are listed in the file BSCHEMES.txt.

Once the run is set up in script (i.e. file paths specified) the script should be started in the working directory where BLOCK.txt, DHR8.txt and the bres files are located (they should all be in the same directory) using the command 

proseeker 

in the MATLAB command line.

Included are MATLAB scripts to allow creation of new BRES files and pre-formatted datasets for functional protein evolution.
