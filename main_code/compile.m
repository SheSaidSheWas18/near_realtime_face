mex -O dramanan_files/resize.cc
mex -O dramanan_files/reduce.cc
mex -O dramanan_files/shiftdt.cc
mex -O dramanan_files/features.cc

% use one of the following depending on your setup.
% 1 is fastest, 3 is slowest.
% If you are using a Windows machine, please use 3. 

% 1) multithreaded convolution using blas
% mex -O dramanan_files/fconvblas.cc -lmwblas -o fconv
% 2) mulththreaded convolution without blas
mex -O dramanan_files/fconvMT.cc -o fconv
% 3) basic convolution, very compatible
% mex -O dramanan_files/fconv.cc

