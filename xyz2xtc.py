from itertools  import islice
from math import sqrt
from optparse   import OptionParser
from subprocess import check_output
from sys        import argv
import numpy as np
import scipy.constants as const
import xtcpy

parser=OptionParser(usage="python %prog -i file.xtc -o outfile -R <Rsph> -d <delta t>\
-n <natoms> -N <Nframes>",version="%prog 0.23")
parser.add_option("-d","--delta-t",dest="deltat",default=1.0,action="store",\
help="time step (ps); default=1.0")
parser.add_option("-i","--in",dest="infile",action="store",help="input xyz file")
parser.add_option("-o","--out",dest="outfile",action="store",\
help="output xtc file",default=False)
parser.add_option("-q","--no-charge",dest="charges",default=True,action="store_false",\
help="if present, do not save charges in a separate file")
parser.add_option("-p","--prec",dest="prec",default=1000.0,action="store",\
help="XTC file precision; default=1000.0")
parser.add_option("-n","--natoms",dest="natoms",default=None,action="store",\
help="number of atoms in system")
parser.add_option("-N","--nframes",dest="nframes",default=False,action="store",\
help="number of frames in system;if not set count them")
parser.add_option("-B","--box-type",dest="boxtype",default="sphere",action="store",\
help="box type: cube or sphere")
parser.add_option("-R","--rad",dest="sphere",default=False,action="store",help="radius or \
box side in nm; default; for spheres write a cubic box with [2R,2R,2R]; if not set then \
R=0.5*(maximum interatomic distance) for each frame (is very slow)")
parser.add_option("-u","--units",dest="unitofL",default="ang",action="store",\
help="unit of lenght to convert from xyz to xtc (bohr or angstrom); default is angstrom to nm")
parser.add_option("-v","--verbose",dest="vb",default=1000.0,action="store",\
help="print conversion progress every x frames; default=1000")
parser.add_option("-s","--shuffle",dest="shuffle",default=None,action="store",\
help="shuffle according to shuffle provided file")
options,args = parser.parse_args(argv[1:])

###### INPUT
if options.outfile is False:
    print "ERROR: missing output file name"
    quit()
else:
    outname = options.outfile+".xtc"
    
Qsaved = False    
 
if options.sphere:
    R = float(options.sphere)

if options.natoms:
    natoms = int(options.natoms)
else:
    print "ERROR: missing number of atoms"
    quit()
    
if options.nframes:
    nframes = int(options.nframes)
else:
    #print "error: missing number of frames"
    #quit()
    print "Counting frames"
    nframes = 0

deltat = float(options.deltat)
prec   = float(options.prec)
vb     = float(options.vb)

###### CONSTANTS
PI = const.pi
DIM = 3
C = (4.0/3.0)*PI
Bohr = const.value("atomic unit of length")
Ang2nm = 0.1
Angstrom = const.angstrom
Bohr2Ang = Bohr/Angstrom
#Bohr2Ang = 0.529177208590
Bohr2nm  = Bohr2Ang/10.0

if options.unitofL == "ang":
    unitofL = Ang2nm
elif options.unitofL == "au":
    unitofL = Bohr2nm
else:
    print "ERROR: unit of lenght not set"
    quit()

##### OBJECTS AND FUNCTIONS
def gen_shuffle(sfile):
    s = open(sfile,"r")
    lines = s.readlines()
    s.close()
    lines = map(lambda x: x.split(),lines)
    slist = list()
    for l in lines:
        slist.append(int(l[2]))
    return len(slist), slist

def shuffle_atoms(natoms,slist,coords):
    newc = np.empty((natoms,3))
    for i in range(natoms):
        newc[slist[i]-1] = coords[i]
    coords[:natoms] = newc[:natoms]
    return coords

def count_lines(F,n):
    """
    count number of lines in file
    """
    a = check_output(["wc","-l",F])
    a = int(a.split()[0])
    b = a/(n+2)
    print "--- file contains ",a," lines, i.e. ",b," frames"
    return b

def calc_box(X,N):
    """
    calculate current box and use it
    """
    r = 0.0
    if options.sphere and options.boxtype=="sphere":
        r = R
    elif options.sphere and options.boxtype=="cube":
        r = 0.5*R
    else:
        for j in xrange(N):
            for k in xrange(j+1,N):
                DX = X[j,0]-X[k,0]
                DY = X[j,1]-X[k,1]
                DZ = X[j,2]-X[k,2] 
                d = 0.5*sqrt(DX*DX+DY*DY+DZ*DZ)
                if d > r:
                    r = d
            r = r + 0.2
    l = 2.0*r
    box = np.array([[l,.0,.0],[.0,l,.0],[.0,.0,l]],dtype=np.float32)
    return box

def ParseXyz(infile,nframes,natoms,Qsaved,nshuffle,slist):
    """
    open xtc parse xyz file frame-wise 
    and write in xtc then close it
    """
    try:
        xyzfile = open(infile,"r")
    except:
        print "ERROR opening xyz file"
        quit()
    xdr1 = xtcpy.xdrfile_open(outname, "w")      
    step = 0
    time = 0.0
    #nrecords = nframes*natoms+2
    lrec = natoms+2
    while True:
# read lines in blocks of natoms + header        
        lines = islice(xyzfile,lrec)
        if lines:
#discard header and check dimension            
            lines = list(lines)[2:]
            if len(lines) != natoms:
                break
#check for charges and convert to numpy arrays        
            if options.charges is True:
                if Qsaved is False:
                    Qfile  = open(options.outfile+"_Q.dat","a")
                    Qsaved = True
                coords  = []
                charges = []
                records = map(str.split,lines)
                for rec in records:
                    a=rec[1:4]
                    coords.append(rec[1:4])
                    charges.append((rec[0],rec[4]))
            elif options.charges is False:
                coords = [i.split()[1:4] for i in lines]
            coords = np.asarray(coords,dtype=np.float32)
            coords.shape = (natoms,DIM)
            coords = unitofL*coords
            box = calc_box(coords,natoms)
            step += 1
            time = step*deltat
#if selected reshuffle atoms
            if nshuffle is not None:
                coords = shuffle_atoms(nshuffle,slist,coords)
            xtcpy.xtc_writeframe(xdr1,coords,step,time,prec,box) 
#read/write charges if present 
            if options.charges:
                Qfile.write("%10d\n Frame %12.6f%12.6f\n" % (natoms,step,time))
                charges = np.asarray(charges,dtype=np.float32)
                charges.shape = (natoms,2)
                np.savetxt(Qfile,charges,fmt="%15.6f")
            if (step % vb)==0 and vb!=0:
                print "--- xyz2xtc: writing frame ",step,"time ",time," ps"
        else:
            print "ERROR writing frame ", step+1
            break
    if step == nframes:
        xtcpy.xdrfile_close(xdr1)
        print "--- xyz2xtc: XTC file closed"    
        if Qsaved:
            Qfile.close()
            print "--- xyz2xtc: charge file closed"
    return None

##### execute
if nframes==0:
    nframes = count_lines(options.infile,natoms)
    
if options.shuffle is not None:
    nshuffle, shuffle_list = gen_shuffle(options.shuffle)
else:
    nshuffle = None
    shuffle_list = None
    
ParseXyz(options.infile,nframes,natoms,Qsaved,nshuffle,shuffle_list)

quit()
