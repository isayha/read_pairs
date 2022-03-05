# String Reconstruction from Read-Pairs Problem Implementation in Python

A Python implementation of the String Reconstruction from Read-Pairs Problem as described in **CPSC 450 - Bioinformatics**.

## Instructions

### How to run

- To run the program, simply run `main.py`, either using the command `python` on the command line (while in the same directory as the source files), or the IDE of your choice
	- For reference, this program was written in Python 3.8.6

- This program requires the provision of 1 (one) command line argument:
	- The command line argument should be the name of the desired input data (text) file (including the `.txt` extension)
    	- e.g.: `example_input_1.txt`
		- The desired input data file should be placed in the root of this project
    		- If the file is placed in a project *subdirectory*, the file can *still* be specified via relative pathing
            	- e.g.: If you were to place `example_input_1.txt` in a new subdirectory called `data`, `.\data\example_input_1.txt` would be the correct argument to provide
      	- The desired input data file should be formatted in the same manner specified in the Assignment Details

- Failure to follow these instructions will result in the program writing a brief error message to console and immediately halting

### Results

- After (successfully) executing the algorithm, the program will:
	- Write the results to console
	- Output the results to `x_output.txt`, where `x` is the name of the specified input data file
		- The output file is found in the same directory as the specified input data file (either the project root or a project subdirectory)
- Note that generating new results for a given data file will overwrite any existing results for the same data file
