#! /usr/bin/env python

import os, sys
import numpy as np
from math import *
from scipy import interpolate
from scipy.interpolate import interp1d, interp2d, RegularGridInterpolator
from collections import OrderedDict
from re import findall 

pref = 'LO-13'

grids = ['CCsame', 'NNsame']
grids = ['C2+C1-', 'C2-C1+', 'NN', 'NC+', 'NC-']
scales = ['mu1','mu2','mu05']

for grid in grids:
    for sc in scales:

        if 'same' in grid:
            dim = 2
        else:
            dim = 3

        outdir = 'lookups'
        tag0 = '{pref}_{grid}'.format(pref=pref, grid=grid)
        tag = '{pref}_{grid}_{sc}'.format(pref=pref, grid=grid, sc=sc)        
        filename = '{tag}.table'.format(tag=tag)

        if dim == 2:
            dir_path = os.path.dirname(os.path.realpath(__file__))
            fpath = os.path.join(dir_path, filename)
            data = np.loadtxt(fpath)
            m1_ar = data[:,0]
            m2_ar = data[:,1]
            m1 = np.array(sorted(list(set(m1_ar))))
            m2 = np.array(sorted(list(set(m2_ar))))            
            m1.dump('{outdir}/{tag}.m1'.format(outdir=outdir, tag=tag0))
            m2.dump('{outdir}/{tag}.m2'.format(outdir=outdir, tag=tag0))
            F_ar = []
            sizeF = np.shape(data)[1] - dim - 1            
            for ii in xrange(sizeF): 
                Fdm = data[:,ii+dim]
                Far = Fdm.reshape(len(m1), len(m2))
                Far.dump('{outdir}/{tag}.F{i}'.format(outdir=outdir, tag=tag, i=ii+1))

        if dim == 3:
            dir_path = os.path.dirname(os.path.realpath(__file__))
            fpath = os.path.join(dir_path, filename)
            data = np.loadtxt(fpath)
            m1_ar = data[:,0]
            m2_ar = data[:,1]
            m3_ar = data[:,2]            
            m1 = np.array(sorted(list(set(m1_ar))))
            m2 = np.array(sorted(list(set(m2_ar))))            
            m3 = np.array(sorted(list(set(m3_ar))))                        
            m1.dump('{outdir}/{tag}.m1'.format(outdir=outdir, tag=tag0))
            m2.dump('{outdir}/{tag}.m2'.format(outdir=outdir, tag=tag0))            
            m3.dump('{outdir}/{tag}.m3'.format(outdir=outdir, tag=tag0))                        
            F_ar = []
            sizeF = np.shape(data)[1] - dim - 1            
            # print filename
            # print len(m1), len(m2), len(m3), np.shape(data[:,1+dim])
            for ii in xrange(sizeF): 
                Fdm = data[:,ii+dim]
                Far = Fdm.reshape(len(m1), len(m2), len(m3))
                F_ar.append( RegularGridInterpolator( (m1, m2, m3), Far) )
                Far.dump('{outdir}/{tag}.F{i}'.format(outdir=outdir, tag=tag, i=ii+1))

