%module pywave2d
%{
  #define SWIG_FILE_WITH_INIT
  #include "pywave2d.hpp"
%}

%include "numpy.i"
%init %{
import_array();
%}

//
// For counts
//
%apply (int *DIM1, int** ARGOUTVIEW_ARRAY1) {(int *n, int **counts)};
%apply (int *DIM1, int** ARGOUTVIEW_ARRAY1) {(int *n, int **indices)};

//
// For get image
//
%apply (int *DIM1, int *DIM2, double** ARGOUTVIEW_ARRAY2) {(int *w, int *h, double **image)};

//
// For get model values
//
%apply (int *DIM1, double** ARGOUTVIEW_ARRAY1) {(int *nv, double **values)};

//
// For set model indices
//
%apply (int DIM1, int* IN_ARRAY1) {(int n, int *indices)};

//
// For set model values
//
%apply (int DIM1, double *IN_ARRAY1) {(int nv, double *values)};

%include "pywave2d.hpp"
