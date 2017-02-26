# Simple NW Alignment Tool
A simple C++ implementation of the Needleman-Wunsch algorithm.
## Usage
./Alignment input_sequence1 input_sequence2 output_result output_profile_1 output_profile_2 step_1 step_2 window_length filter_threshold

 - **input_sequence1**, *file*, input, containing the first sequence;
 - **input_sequence2**, *file*, input, containing the second sequence;
 - **output_result**, *file*, main output, has 3 columns;
 - **output_profile_1**, *file*, output, has 1 column;
 - **output_profile_2**, *file*, output, has 1 column;
 - **step_1**, *int*, input, step size 1;
 - **step_2**, *int*, input, step size 2;
 - **window_length**, *int*, input, window size;
 - **filter_threshold**, *int*, input, filter threshold;
 
## Example
With input data from the `data` folder, run `./Alignment data/input1 data/input2 output1 output2 output3 1 1 100 0`.
Check that the output files (`output1`, `output2` and `output3`) match those in the `data` folder.

## Previous Version
The four integer parameters were part of the first input file as shown in the `data_ott` folder
