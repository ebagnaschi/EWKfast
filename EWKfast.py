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
        input_file = sys.argv[1] 
        slha_file = sys.argv[2]
    except:
        print '[input file] [SLHA file]'
        exit()

    dir_path = os.getcwd()
    slha_path = os.path.join(dir_path, slha_file)
    params = get_params(slha_path)
    params, range_memo = ranges(params)

    input_path = os.path.join(dir_path, input_file)
    input_list = process_input(input_path)
    for li in input_list: 
        rs, order, mode, grid, indices = li['rs'], li['order'], li['mode'], li['grid'], li['indices']
        res, warn = get_xsec(params, mode=grid, indices=indices, table=Tabs.tables[(order,rs,grid)])
        print rs, order, mode, res, warn

    if len(range_memo) > 0:
        print ''
        for memo in range_memo: print memo

    exit()

    if mode == 'c1c1': runmode, mass, indices = 'CCsame', params['mC'][0], [0,0]
    if mode == 'n1n1': runmode, mass, indices = 'NNsame', params['mN'][0], [0,0]

    res = get_xsec(params, mode=runmode, indices=indices, table=Tabs.tables[(order,rs,runmode)])
    print '{mode} {rs} {order} xsec: '.format(mode=mode, rs=rs, order=order), res




