Step 1: 'make' inside the directory named 'demo'. Please not that the following dependencies need to be installed in the 		 system for the 'make' command to run succesfully:
		1. g++ (4.9.2) [sudo apt install build-essential]
		2. Boost C++ library (1.58.0) [sudo apt-get install libboost-all-dev]
		3. libnuma (2.0.9) [sudo apt-get install -y libnuma-dev]
		4. libtcmalloc_minimal in google-perftools (2.1) [sudo apt install -y google-perftools libgoogle-perftools-dev]

Step 2: Prepare the input graph file such that it is in the form of an edgelist (similar format as the datasets on SNAP). 		  Please note that the indexing of the nodes in the edgelist must be continuous and start from 0, otherwise the 			program may consume more memory and show lower performance (the python script for this transformation has been 			included with Gorder.

Step 3: The command to run Rabbit-Order is: ./reorder /path/to/input/file >> /path/to/output/vertex/file
		The output will be written to the text file specified in the form of a list of vertex ids (1 per line) where the i-th line will have the new label of the vertex whose original id was 'i'.
