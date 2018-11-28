#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 10 17:08:21 2012

@author: gmancini
"""

# reads Gaussian formatted checkpoint and writes xtc trajectory
# NB a proper coordinate file (.pdb, .gro) may be needed to
# read the xtc

#check path for xtcpy module
import socket,sys
host = socket.gethostname()
if host=="wintermute" or host=="magnusthered":
  sys.path.append("/home/gmancini/SNS/SnSDev/script/libxtcpy")
elif host=="neuromancer":
  sys.path.append("/home/gmancini/script/libxtcpy")
else: 
  sys.path.append("/home/gmancini/Work/psdev/devel/intel/script/libxtcpy")

#import & help
from csv import reader
from math import sqrt
from optparse import OptionParser
from subprocess import call
from sys import argv
import numpy as np
import re
import xtcpy

parser=OptionParser(usage="python %prog <options> <input>.<[f]chk/xyz>",version="%prog 0.23")
parser.add_option("-b","--box-type",dest="boxtype",default=True,action="store_false",\
help="calculate box at each step")
#parser.add_option("-B","--Box",dest="Box",default=False,action="store",\
#nargs=3,help="use these values for box")
parser.add_option("-d","--delta-t",dest="dt",default=0.002,action="store",\
help="time step; default=0.002")
parser.add_option("-p","--prec",dest="prec",default=1000.0,action="store",\
help="XTC file precision; default=1000000.0")
parser.add_option("-s","--skip",dest="skip",default=0,action="store",\
help="write one frame every x frames; default=0")
parser.add_option("-S","--symbols",dest="symbols",help="2 column file with atomic number and atom names",\
default=False)
parser.add_option("-t","--type",dest="filetype",default="fchk",action="store",\
help="file type:chk/fchk/xyz; default fchk")
#parser.add_option("-u","--input-um",dest="um",default='bohr',action="store",\
#help="input file units: ang/bohr")
#parser.add_option("-U","--output-um",dest="UM",default='nm',action="store",\
#help="output file units: ang,nm; chk is in bohr")
parser.add_option("-v","--verbose",dest="vb",default=1000.0,action="store",\
help="print conversion progress every x frames; default=1000")
parser.add_option("-x","--xyz",dest="xyz",default=False,action="store_true",\
help="print 1st frame in xyz format")
options,args = parser.parse_args(argv[1:])

#process options
#if options.Box is not False:
#    Box = map(float,options.Box)
#    a   = Box[0]
#    b   = Box[1]
#    c   = Box[2]
#    Box = np.array([[a,.0,.0],[.0,b,.0],[.0,.0,c]],dtype=np.float32)
#else:
#    Box = False
Box = False
dt   = float(options.dt)        
prec = float(options.prec)
skip = int(options.skip)
if skip==0:
    skip = False
#um   = options.um
#UM   = options.UM
vb   = float(options.vb)
input_file = args[0]
ft = options.filetype
if options.symbols is not False:
    fsymb = open(options.symbols,"r")
    AtomNames = dict()
    for i in fsymb.readlines():
        line = i.split()
        AtomNames[int(line[0])] = line[1]
#    with fsymb as f:
#        AtomNames = dict(reader(f, delimiter=' '))

#conversion
Bohr2Ang = 0.529177208590
Bohr2nm  = Bohr2Ang/10.0
nm2bohr  = 1.0/Bohr2nm
ang2bohr = 1.0/Bohr2Ang
nm2ang   = 10.0
ang2nm   = 0.1
dim      = 3

class Trajectory:
    """
    Basic container class for data
    Checkpoint files are read-only!
    """
    def __init__(self,filename,filetype,box,bfixed):
        """
        Create trajectory object and parse imput file
        G09 util formchk must be in the path to convert
        binary checkpoint files on the fly
        Box for xtc file is either given in input or 
        generated using maximum distance between atoms
        """
        self.filename = filename
        if filetype == "xyz":
            nf,nat,nc = self.XYZParse(filename,filetype)
            self.Coordinates = ang2nm*self.Coordinates
        else:
            nf,nat,nc = self.ChkParse(filename,filetype)
            self.Coordinates = Bohr2nm*self.Coordinates
        self.NAtoms  = nat
        self.NCoords = nc
        self.NFrames = nf
        if box is not False:
            self.box = box
        elif bfixed is True:
            f = self.Coordinates[0,:]
            f.shape = (self.NAtoms,3)
            self.box = self.genbox(f)
        else:
            self.box = None
            
    def XYZParse(self,name,ftype):
        """
        read and parse frames.xyz
        """
        rg = r'(.*)\.(xyz)'
        self.chkprefix = re.search(rg,name).group(1)
        print "--- chk2xtc: reading input file"        
        xyzfile = open(name,"r")
        records = xyzfile.readlines()        
        nat = int(records[0])
        ncoords = nat*4 
        records = map(lambda x:x.split(),records)
        nframes = len(filter(lambda x: x[0]=='Frame',records))
#        print nframes,ncoords
        print "--- chk2xtc: getting coordinates"
        records = filter(lambda x:(len(x)>1 and x[0]!='Frame'),records)
#        print records[0]
        coords  = np.array(records,dtype=np.float32)
#       print coords[0]
#       print coords.shape
        coords.shape = (nframes,ncoords)
#        print coords[0,:]        
        atoms = coords[0,::4].astype('int')
#       print atoms
        self.Atoms = dict(zip(xrange(nat),atoms))
        self.Coordinates = np.delete(coords, np.s_[::4],1)
#       print self.Coordinates
        ncoords = nat*3
#        print ncoords        
#        quit()
        print "--- chk2xtc: parsing input file done"
        print "--- chk2xtc: checkpoint file contains:"
        print "                                     ",nat," atoms"
        print "                                     ",nframes," frames"
        return nframes,nat,ncoords
                    
    def ChkParse(self,chkname,ftype):
        """
        read .chk file and filt 
        coordinates and mol spec
        """
        rg = r'(.*)\.(f?chk)'
        self.chkprefix = re.search(rg,chkname).group(1)
        if ftype is "chk":
            print "--- chk2xtc: convert binary checkpoint to ASCII"
            call(["formchk",chkname])
            chkfile = open(self.chkprefix+".fchk","r")
        else:
            chkfile = open(chkname,"r")
        print "--- chk2xtc: reading input file"
        chk = chkfile.read()
# find atomic numbers      
        print "--- chk2xtc: getting coordinates"        
        regexp = r'Atomic numbers\s+I\s+N=\s+(\d+)\n(.*?)Nuclear'
        atoms  = re.search(regexp,chk,re.S)
        NAtoms  = int(atoms.group(1))
        NCoords = NAtoms*3
        #atoms = map(int,(atoms.group(2)).split())
        atoms = (atoms.group(2)).split()
        if len(atoms) != NAtoms:
            print "error parsing atomic numbers"
            quit()
        self.Atoms   = dict(zip(xrange(NAtoms),atoms))
# number of frames
        regexp  = r'Trajectory MaxStp\s+I\s+(\d+)'
        NFrames = int(re.search(regexp,chk).group(1))-1
        print "--- chk2xtc: parsing input file done"
        print "--- chk2xtc: checkpoint file contains:"
        print "                                     ",NAtoms," atoms"
        print "                                     ",NFrames," frames"
#read frames in numpy array  
        regexp  = r'Traj num\s+1\s+Geometries\s+R\s+N=\s+\d+\n(.*?)Traj num'
        frames = (re.search(regexp,chk,re.S).group(1)).split()
        self.Coordinates = np.array(frames,dtype=np.float32)
        self.Coordinates.shape = (NFrames,NCoords)
        return NFrames,NAtoms,NCoords
        
    def genbox(self,f):
        """
        calculate maximum distance between atoms and
        generate box to write xtc file
        """
        r = 0.0
        for j in xrange(self.NAtoms):
            for k in xrange(j+1,self.NAtoms):
                DX = f[j,0]-f[k,0]
                DY = f[j,1]-f[k,1]
                DZ = f[j,2]-f[k,2] 
                d = 0.5*sqrt(DX*DX+DY*DY+DZ*DZ)
                if d > r:
                    r = d
#  cubic  box
        l = 2.5*r
        box = np.array([[l,.0,.0],[.0,l,.0],[.0,.0,l]],dtype=np.float32)
        return box
        
    def writexyz(self,frame):
        """
        write given frame in xyz format
        in a new file
        """
        coords = self.Coordinates[frame,:]*10.0
        coords.shape = (self.NAtoms,dim)
        #coords = coords+0.4*np.diag(self.box)
        xyzname = "out_" + self.chkprefix+".xyz"
        xyz = open(xyzname,"w")
        xyz.write("%d\n created by chk2xtc\n" % self.NAtoms)
#better would be to create a list and use writelines(?)
        for i in xrange(self.NAtoms):
            atom = AtomNames[self.Atoms[i]]
            x = coords[i,0] 
            y = coords[i,1] 
            z = coords[i,2] 
            xyz.write("%s %16.6f %12.6f %12.6f\n" % (atom,x,y,z))
        xyz.close()
        print "--- chk2xtc: first frame saved in xyz format"
        return None
        
    def writextc(self,dt,prec,skip,vb,xyz):
        """
        writes xtc file; first frame saved in xyz
        """
        xtctraj = self.chkprefix+".xtc"    
        xdr1 = xtcpy.xdrfile_open(xtctraj, "w")
        print "--- chk2xtc: XTC file open"
        if xyz is True: 
            self.writexyz(0)
        for frame in xrange(self.NFrames):
            if skip!=0 and frame%skip==0:
                continue
            if frame%vb==0:
                print "--- chk2xtc: writing frame ",frame
            time = float(frame)*dt
            this_frame = self.Coordinates[frame,:]
            this_frame.shape = (self.NAtoms,dim)
            if self.box is None:
                box = self.genbox(this_frame)
            else:
                box = self.box
            this_frame = this_frame #+ 0.4*np.diag(box)
            #if frame%vb==0:
            #    print "--- chk2xtc: coordinates shifted by ",0.5*np.diag(box)
            xtcpy.xtc_writeframe(xdr1,this_frame,frame,time,prec,box) 
        xtcpy.xdrfile_close(xdr1)
        print "--- chk2xtc: XTC file closed"
        return None
          
trj = Trajectory(input_file,ft,Box,options.boxtype)
trj.writextc(dt,prec,skip,vb,options.xyz)

quit()
    
    
        
