# Introduction

This Python module provides a wrapper around C/C++ code for trans-dimensional
trees with a wavelet parameterization described in

Hawkins R. and Sambridge M., "Geophysical imaging using trans-dimensional trees",
Geophysical Journal International, 2015.

The Numpy swig interface file, numpy.i, is included in this source bundle for
convenience. See numpy.i for their copyright and license.

Released Opensource under GPL v3 License. Copyright (c) 2014 - 2018 Rhys Hawkins.
See LICENSE for details.

# Prequisites

For building the module:

* gcc
* g++ (Tested on 6.3 and 7.2)
* Python 2.7
* SWIG
* Numpy
* gmp (Gnu Arbitrary precision math library)
* gsl (Gnu scientific library)
* Trans-dimensional Tree base libraries (see TDTbase repo)


For plotting:

* Matplotlib

# Building the module

It is assumed that the TDTbase repository is extracted under the same parent
folder as pywave2d, e.g.

* Parent Directory
..* TDTbase
..* pywave2d

If this is not the case, then setup.py needs to be modified to point to the
correct location of the TDTbase repository (see the setting of the TDTBASE
variable at the beginning of setup.py).

```
tar -xzf pywave2d.tar.gz
cd pywave2d
python2 setup.py build
```

This module is still experimental and is probably safest to run in the
build directory for now. See the imports section of example/regression.py
for how to import the module from a known subdirectory.

# Running the example

The example is a simple 2D regression problem in which the true model is
a Gaussian hill. The synthetic data with independent Gaussian noise is created
with the mkdata.py script.

The inversion is performed by regression.py and includes hierarchical error
scale estimation and hierarchical prior steps.

The regression.py file outputs the entire ensemble to a python pickle file
and various statistics from the ensemble can be computed. This is done
in the postprocess.py script which computes the ensemble mean and
some other diagnostics.

Various simple plotting scripts show how to plot some of the outputted files.

To run the entire example, starting from the the pywave2d directory, do:

```
cd example
python2 mkdata.py

python2 regression.py
python2 postprocess.py

python2 plot_mean.py
python2 plot_prior.py
python2 plot_like.py
python2 plot_coeff.py
```

# Notes

A work in progress so most of the module documenation is in comments in the
regression.py script.






