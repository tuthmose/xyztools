#Makefile for xdrlib with python SWIG interface
#Giordano Mancini Sept 2012

# Set compiler flags 
CC     = gcc
CFLAGS  = -O3 -ffast-math -fomit-frame-pointer -pipe -fPIC
F77 = gfortran
FFLAGS = $(CFLAGS)

SHELL	= /bin/sh

OBJECTS = xdrfile.o xdrfile_c_test.o xdrfile_xtc.o xdrfile_trr.o 
INCLUDE = xdrfile.h xdrfile_trr.h xdrfile_xtc.h 

#CPPFLAGS = -I/home/g.mancini/pkg/python27/include/python2.7 -I/home/g.mancini/pkg/python27/lib/python2.7/site-packages/numpy/core/include/
CPPFLAGS = -I/usr/include/python2.7/ -I/usr/lib64/python2.7/site-packages/numpy/core/include/

SWIG     = swig 
SWIGOPT  = -python -Wall
SWIGOBJS = xtcpy_wrap.o

all: clean xtcpy

check: clean xtcpy test

test:
	python test_xtc.PY 

xtcpy: $(OBJECTS) $(SWIGOBJS)
	$(CC) -shared -o _xtcpy.so $(OBJECTS) $(SWIGOBJS)
	
staticlib:  trr2xtc.o $(OBJECTS)
	ar rcs libxtc.a $(INCLUDE) $(OBJECTS)
	
xtcpy_wrap.c: xtcpy.i
	$(SWIG) $(SWIGOPT) xtcpy.i
	
xtcpy_wrap.o: xtcpy_wrap.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -c xdrfile.c xtcpy_wrap.c

clean:
	rm -f *.so *.o *.py *.pyc *.a *.gch xtcpy_wrap.c test.xtc
	
.SUFFIXES : .c .h .o

.c.o:
	$(CC) $(INCLUDE) $(CFLAGS) -c $*.c
