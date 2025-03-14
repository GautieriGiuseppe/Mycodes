This folder contains two versions of a Multiple Sequence Aligner script.
The standard one makes usage of both ClustalW and Muscle through the subprocess module
while the built_in performs the alignment in local resembling the ClustalW algorithm
and storing everything using NumPy arrays and Pandas dataframes.

The v2 of the python script uses multiprocessing to run in parallel ClustalW and Muscle, enhancing the computation time.
Now is also present in the Java version.
The Java version is slightly faster than the v2, while the fastest is the built-in.
