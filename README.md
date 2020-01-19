# Mittag-Leffler for Stata
 Mittag-Leffler function in Stata/Mata using Garrappa (2015) inverse Laplace transform algorithm
 
This is a direct port of Roberto Garrappa's MATLAB implementation of his inverse-Laplace-transform method for computing the one-, two-, or three-parameter Mittag-Leffler function. The package provides Mata functions, mlf() and LTInversion(), with identical syntax and behavior as the MATLAB original.
 
Source:
R. Garrappa, Numerical evaluation of two and three parameter
  Mittag-Leffler functions, SIAM Journal of Numerical Analysis, 2015,
  53(3), 1350-69.

Documentation for the MATLAB original:
  mathworks.com/matlabcentral/fileexchange/48154-the-mittag-leffler-function

To use this file in Stata, more precisely in Mata, run MattagLeffler.mata one time as if it were a do file. 
Running this file  will create the compiled module lMattagLeffer.mlib, and make the Mata function mlf() 
permanently available.

No additional documentation is provided.

The library includes a function verify_mlf(), which validates the program against Tran Quoc Viet's excellent collection
of test values (https://github.com/tranqv/Mittag-Leffler-function-and-its-derivative/tree/master/tcases). This function
takes a string argument, which is the path to the directory containing the test files. At least on one Intel-based Windows PC,
all results match the target values up to a relative difference of at most 1.75e-12.
