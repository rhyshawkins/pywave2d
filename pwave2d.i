%module pywave2d
%{
  #define SWIG_FILE_WITH_INIT
  #include "pywave2d.hpp"
%}

%include "numpy.i"
%init %{
import_array();
%}

%apply (int DIM1, double* INPLACE_ARRAY1) {(int n0, double *a0)};
%apply (int DIM1, double* IN_ARRAY1) {(int n, double *a), (int m, double *b)};
%apply (int DIM1, double* ARGOUT_ARRAY1) {(int size, double *arr)};

%include "pywave2d.hpp"
