from distutils.core import setup, Extension
import numpy
import os

TDTBASE='../TDTbase'

os.environ['CC'] = 'g++';
setup(name='pyWave2D',
      version='1.0',
      ext_modules = [Extension('_pywave2d',
                               ['pywave2d.cpp', 'pywave2dexception.cpp', 'pywave2d.i'],
                               include_dirs = [os.path.join(TDTBASE, 'hnk'),
                                               os.path.join(TDTBASE, 'log'),
                                               os.path.join(TDTBASE, 'oset'),
                                               os.path.join(TDTBASE, 'sphericalwavelet'),
                                               os.path.join(TDTBASE, 'tracking'),
                                               os.path.join(TDTBASE, 'wavelet'),
                                               os.path.join(TDTBASE, 'wavetree')],
                               swig_opts=['-c++'],
                               libraries=['gsl', 'gslcblas', 'gmp', 'm'],
                               extra_objects=[os.path.join(TDTBASE, 'log/liblog.a'),
                                              os.path.join(TDTBASE, 'hnk/libhnk.a'),
                                              os.path.join(TDTBASE, 'sphericalwavelet/libsphericalwavelet.a'),
                                              os.path.join(TDTBASE, 'tracking/libtracking.a'),
                                              os.path.join(TDTBASE, 'wavelet/libwavelet.a'),
                                              os.path.join(TDTBASE, 'wavetree/libwavetree.a'),
                                              os.path.join(TDTBASE, 'oset/liboset.a')])])
