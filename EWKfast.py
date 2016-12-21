#! /usr/bin/env python

import os, sys
import numpy as np
from math import *
from pyslha import * 
from functions import *
from scipy import interpolate
from scipy.interpolate import interp1d, interp2d, RegularGridInterpolator

if __name__ == "__main__":

    Tabs = Get_Tables()
    Tabs.add( key=('LO','13','CCsame'), filename='LO-13_CCsame.table', dim=2, lenF=6 )
    Tabs.add( key=('LO','13','NNsame'), filename='LO-13_NNsame.table', dim=2, lenF=5 )
    Tabs.add( key=('NLO','13','CCsame'), filename='NLO-13_CCsame.table', dim=2, lenF=6 )
    Tabs.add( key=('NLO','13','NNsame'), filename='NLO-13_NNsame.table', dim=2, lenF=5 )

    try:
        SLHAfile = sys.argv[1]
        mode = sys.argv[2]
        rs = sys.argv[3]
        order = sys.argv[4]        
    except:
        print 'Input [SLHAfile] [mode (c1c1, n1n1, all)] [rootS] [NLO, LO]'
        exit()

    dir_path = os.getcwd()
    #dir_path = os.path.dirname( os.path.realpath(__file__) )
    SLHApath = os.path.join(dir_path, SLHAfile)
    params = get_params(SLHApath)
    if params['mQ'] >= 6000: params['mQ'] = 5999.

    if mode == 'c1c1': runmode, mass, indices = 'CCsame', params['mCha'][0], [0,0]
    if mode == 'n1n1': runmode, mass, indices = 'NNsame', params['mNeu'][0], [0,0]

    res = get_xsec(params, mode=runmode, indices=indices, table=Tabs.tables[(order,rs,runmode)])
    print '{mode} {rs} {order} xsec: '.format(mode=mode, rs=rs, order=order), res




