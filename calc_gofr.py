#!/usr/bin/env python

import math
import argparse

def nrst_long(argm):
    if(argm>0):
        return int(argm+0.5)
    else:
        return int(argm-0.5)

parser = argparse.ArgumentParser(description = 'calculate radial distribution function from XYZ file')
parser.add_argument('-i', help='xyz file name',type=str, default='WATER-pos-1.xyz')
parser.add_argument('-atom1', help='atom type to compute g(r), options are atomic symbol or all',type=str,default='all')
parser.add_argument('-atom2', help='atom type to compute g(r), options are atomic symbol or all',type=str,default='all')
parser.add_argument('-box', help='size of box in units of length',type=float,default=12.42)
parser.add_argument('-o', help='output file name for g(r)',type=str,default='gr.dat')

args = parser.parse_args()

xyzfile = args.i
atomtype1 = args.atom1
atomtype2 = args.atom2
boxs = args.box
grdatafile = args.o

print('Coordinate file name: ',xyzfile)
print('will read atom type 1: ',atomtype1)
print('will read atom type 2: ',atomtype2)
print('simulation box size: ',boxs)
print('output file name: ',grdatafile)

f = open(xyzfile,'r')


delr = 0.05   # histogram bin width
scal = 0.5   # box dimension scaling factor
#boxs = 12.42  # box size

nbin = int(boxs*scal/delr)

hist = []
# initialize HISTOGRAM
for i in range(nbin):
    hist.append(0.0)

print('number of histogram bins', nbin)

invd = 1.0/delr
invb = 1.0/boxs
rmax = float(nbin*delr)

## first check is atom1 is same as atom2 ###

sametype = False
difftype = False
alltype = False

if(atomtype1 == 'all'):
    alltype = True
    if(atomtype2 != 'all'):
        print('WARNING: Atom type 1 is set to all atoms, so atom type 2 must also be set to all. Changing atomtype 2 to all')
        atomtype2 = 'all'

if(atomtype2 == 'all'):
    if(atomtype1 != 'all'):
        print('WARNING: Atom type 2 is set to all atoms, but atom type 1 type: ', atomtype1)
        print('Changing atomtype 2 to ',atomtype1)
        print('To compute radial distribution between all pairs of atoms, use -atom1 all')
        atomtype2 = atomtype1

if(not alltype):
    if(atomtype1 == atomtype2):
        sametype = True

    if(atomtype1 != atomtype2):
        difftype = True

if(sametype):
    print('Computing radial distribution between same atom types', atomtype1, 'and', atomtype2)

if(difftype):
    print('Computing radial distribution between different atom types', atomtype1, 'and', atomtype2)

if(alltype):
    print('Computing radial distribution between all atom types')

nframes = 0
if(alltype):
    coordx = []
    coordy = []
    coordz = []
    for line in f:
        line=line.strip()
        columns=line.split()
        s=columns[0]
        if(s.isdigit()):
            nframes += 1
        else:
            if(s != 'i'):
                coordx.append(float(columns[1]))
                coordy.append(float(columns[2]))
                coordz.append(float(columns[3]))

    print(nframes)
    print(len(coordx))
    natoms = int(len(coordx)/nframes)
    print(natoms, ' atoms per frame')
    for t in range(nframes):
        print('cycle ', t+1, ' of ', nframes)
        for i in range(natoms-1):
            m=natoms*t+i
            for j in range(i+1,natoms):
                n=natoms*t+j
                xdis = coordx[m] - coordx[n]
                ydis = coordy[m] - coordy[n]
                zdis = coordz[m] - coordz[n]
                xdis = xdis - boxs*(float(nrst_long(xdis*invb)))
                ydis = ydis - boxs*(float(nrst_long(ydis*invb)))
                zdis = zdis - boxs*(float(nrst_long(zdis*invb)))
                dist = math.sqrt(xdis*xdis + ydis*ydis + zdis*zdis)
                if(dist<rmax):
                    ibin = int(dist*invd)
                    hist[ibin] = hist[ibin] + 2.0

if(sametype):
    coordx = []
    coordy = []
    coordz = []
    for line in f:
        line=line.strip()
        columns=line.split()
        s=columns[0]
        if(s.isdigit()):
            nframes += 1
        if(s==atomtype1):
            coordx.append(float(columns[1]))
            coordy.append(float(columns[2]))
            coordz.append(float(columns[3]))

    print(nframes)
    print(len(coordx))
    natoms = int(len(coordx)/nframes)
    print(natoms, ' atoms per frame')
    for t in range(nframes):
        print('cycle ', t+1, ' of ', nframes)
        for i in range(natoms-1):
            m=natoms*t+i
            for j in range(i+1,natoms):
                n=natoms*t+j
                xdis = coordx[m] - coordx[n]
                ydis = coordy[m] - coordy[n]
                zdis = coordz[m] - coordz[n]
                xdis = xdis - boxs*(float(nrst_long(xdis*invb)))
                ydis = ydis - boxs*(float(nrst_long(ydis*invb)))
                zdis = zdis - boxs*(float(nrst_long(zdis*invb)))
                dist = math.sqrt(xdis*xdis + ydis*ydis + zdis*zdis)
                if(dist<rmax):
                    ibin = int(dist*invd)
                    hist[ibin] = hist[ibin] + 2.0

if(difftype):
    coord1x = []
    coord1y = []
    coord1z = []
    coord2x = []
    coord2y = []
    coord2z = []
    for line in f:
        line=line.strip()
        columns=line.split()
        s=columns[0]
        if(s.isdigit()):
            nframes += 1
        if(s==atomtype1):
            coord1x.append(float(columns[1]))
            coord1y.append(float(columns[2]))
            coord1z.append(float(columns[3]))
        if(s==atomtype2):
            coord2x.append(float(columns[1]))
            coord2y.append(float(columns[2]))
            coord2z.append(float(columns[3]))

    print(nframes)
    print(len(coord1x))
    print(len(coord2x))
    natoms1 = int(len(coord1x)/nframes)
    natoms2 = int(len(coord2x)/nframes)
    print(natoms1, ' type 1 atoms per frame')
    print(natoms2, ' type 2 atoms per frame')
    for t in range(nframes):
        print('cycle ', t+1, ' of ', nframes)
        for i in range(natoms1):
            m=natoms1*t+i
            for j in range(natoms2):
                n=natoms2*t+j
                xdis = coord1x[m] - coord2x[n]
                ydis = coord1y[m] - coord2y[n]
                zdis = coord1z[m] - coord2z[n]
                xdis = xdis - boxs*(float(nrst_long(xdis*invb)))
                ydis = ydis - boxs*(float(nrst_long(ydis*invb)))
                zdis = zdis - boxs*(float(nrst_long(zdis*invb)))
                dist = math.sqrt(xdis*xdis + ydis*ydis + zdis*zdis)
                if(dist<rmax):
                    ibin = int(dist*invd)
                    hist[ibin] = hist[ibin] + 1.0

## Normalize  and print histogram

if(difftype):
    natoms2 = (natoms1*natoms2)
else:
    natoms2 = natoms*natoms

fout = open(grdatafile,'w')
for i in range(nbin):
    shex = (i+1)*(i+1)*(i+1)
    shin = i*i*i
    volu = 4.0/3.0*math.pi*(shex-shin)*delr*delr*delr
    hist[i] = hist[i]/(float(natoms2*nframes*volu*invb*invb*invb))
    rval = float((i+0.5)*delr)
    txt = str(rval) + ' ' + str(hist[i])+'\n'
    fout.write(txt)

f.close()
fout.close()
