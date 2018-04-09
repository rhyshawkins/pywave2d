
TDTBASE=../TDTbase
INCLUDES = \
	-I/usr/include/python2.7 \
	-I$(TDTBASE)/log \
	-I$(TDTBASE)/oset \
	-I$(TDTBASE)/hnk \
	-I$(TDTBASE)/sphericalwavelet \
	-I$(TDTBASE)/tracking \
	-I$(TDTBASE)/wavelet \
	-I$(TDTBASE)/wavetree

CXX = g++
CXXFLAGS = -c -g -Wall $(INCLUDES) -fPIC

TARGETS = 

OBJS = pywave2d.o \
	pywave2dexception.o

STATICLIBS = $(TDTBASE)/log/liblog.a \
	$(TDTBASE)/hnk/libhnk.a \
	$(TDTBASE)/sphericalwavelet/libsphericalwavelet.a \
	$(TDTBASE)/tracking/libtracking.a \
	$(TDTBASE)/wavelet/libwavelet.a \
	$(TDTBASE)/wavetree/libwavetree.a \
	$(TDTBASE)/oset/liboset.a

LIBS = $(shell gsl-config --libs) \
	-lgmp

all : $(TARGETS)

SOURCE = Makefile \
	README.md \
	LICENSE \
	numpy.i \
	pywave2d.cpp \
	pywave2d.hpp \
	pywave2d.i \
	pywave2dexception.cpp \
	pywave2dexception.hpp \
	setup.py \
	tests/pywave2d_instantiate.py \
	example/mkdata.py \
	example/priorproposal.txt \
	example/regression.py \
	example/postprocess.py \
	example/chainhistory.py \
	example/plot_mean.py \
	example/plot_prior.py \
	example/plot_like.py \
	example/plot_coeff.py \
	example/plot_khist.py 

EXTRADIST = 

INSTALL = install
INSTALLFLAGS = -D

DATE = $(shell date +"%Y%m%d%H%M")
DIR = pywave2d
TGZ = $(DIR).tar.gz

dist :
	mkdir -p $(DIR)
	echo $(DATE) > $(DIR)/Version
	for f in $(SOURCE) $(EXTRADIST); do \
	    $(INSTALL) $(INSTALLFLAGS) $$f $(DIR)/$$f ; \
	done
	tar -czf $(TGZ) $(DIR)/*
	rm -rf $(DIR)


%.o : %.cpp
	$(CXX) $(CXXFLAGS) -o $*.o -c $*.cpp

clean :
	rm -f $(TARGETS) *.o

