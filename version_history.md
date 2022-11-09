# Version History

## 0.41

- Replaced the python test file with a inbuilt tester function in the C program.
- Deleted old result csv files.

## 0.33

Will now give each thread a group of elements to compute the averages on instead of a group of rows. This is done to reduce the difference between time taken for the fast thread and the slowest thread to compute the average of four values.

## 0.32

Updated Comments and documentaion.

## 0.31

- Changed code to repeatedly calculate the average of the four neighbors until the diffrence between the previous average and current average is less than a given precision level.
- Removed the function to calculate the average of four neighbours and included the code where the function is called.
- Updated Python batch test to include precision level.
- Removed decimal place support.
- Updated function to calculate expected output arrayto support the precision level.

## 0.26

- Added decimals place support.
- Added optional arguments.

## 0.25

- Commenting added.
- Deleted a series of result csv files.

## 0.24

First implementation to pass all 100 tests for each 44 thread producing the correct output.

## 0.23

Code now has a method to compare calculated results with expected results. The outcome of this comparison is printed to the console.

## 0.22

Improved By batch testing.

## 0.21

Code will now distributes rows to each thread instead of cells.

## 0.12

Updated the results logging file. Will now no longer overwrite the previous results file.

## 0.11

- Included timing functions in .c
- Fixed "pthread" include error.

## 0.10

Implimented a Python batch test script

## 0.02

First successfully execution of the program.

## 0.01

Init Repo
