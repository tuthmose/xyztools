/* 
SWIG interface file for xdrlib
Giordano Mancini Sept 2012
*/

// start with SWIG module
%module xtcpy
%{
#define SWIG_FILE_WITH_INIT
// XTC headers
#include "xdrfile.h"
#include "xdrfile_trr.h"
#include "xdrfile_xtc.h"
%}

// include interface for Numpy arrays 
%include "numpy.i"

%init %{
import_array();
%}

/* 
import stuff only for open/close and high level functions
*/
enum { exdrOK, exdrHEADER, exdrCLOSE, exdrMAGIC,
        exdrNOMEM, exdrENDOFFILE, exdrFILENOTFOUND};

extern XDRFILE* xdrfile_open(const char *path, const char *mode);
extern int xdrfile_close(XDRFILE *fp);

#define DIM 3
typedef float matrix[DIM][DIM];
typedef float rvec[DIM];

//map input arrays and wrap write_xtc
%apply (float IN_ARRAY2[ANY][ANY]) {(matrix box)};
%apply (int DIM1, int DIM2, float *IN_ARRAY2) {(int natoms,int dim,float *XIN)};

%inline %{
        int xtc_writeframe(XDRFILE *xd,int natoms, int dim,float *XIN,int step,
                float time, float prec,matrix box){

                int out;
                out = write_xtc(xd,natoms,step,time,box,(rvec *)XIN,prec);
                if (out == exdrOK){
                        return 0;
                } else {
                        return -1;
                }
        }
%}

// wrap read_xtc_natoms
%inline{
        int read_xtcnatoms(char *filename){
 
                int out, natoms;
                out = read_xtc_natoms(filename,&natoms);
                if (out == exdrOK){
                        return natoms;
                } else {
                        return -1;
                }
        }
}

// map output arrays 
%apply (float *INPLACE_ARRAY2, int DIM1, int DIM2) {(float *XOUT,int natoms,int dim)};
%inline %{
        int xtc_readframe(XDRFILE *xd,float *XOUT,int natoms, int dim,int step,
                float time, float prec,matrix box){

                int out;
                out = read_xtc(xd,natoms,&step,&time,box,(rvec *)XOUT,&prec);
                if (out == exdrOK){
                        return 0;
                } else {
                        return -1;
                }
        }
%}

//clear mappings
%clear (matrix box);
%clear (int natoms,int dim, float *XIN);
%clear (float *XOUT,int natoms,int dim);

