#! /usr/bin/env python

import os, sys
import numpy as np
from math import *
from pyslha import * 
from functions import *
from scipy import interpolate
from scipy.interpolate import interp1d, interp2d, RegularGridInterpolator

if __name__ == "__main__":

    #Tabs = Get_Tables()
    # Tabs.add( key=('LO','13','CCsame'), filename='LO-13_CCsame.table', dim=2, lenF=6 )
    # Tabs.add( key=('LO','13','NNsame'), filename='LO-13_NNsame.table', dim=2, lenF=5 )
    # Tabs.add( key=('NLO','13','CCsame'), filename='NLO-13_CCsame.table', dim=2, lenF=6 )
    # Tabs.add( key=('NLO','13','NNsame'), filename='NLO-13_NNsame.table', dim=2, lenF=5 )

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
    input_list, options = process_input(input_path)

    Tabs = load_tables(input_list, options) 

    #show_Fs('LO','13','NC+','mu1', [0,591,391], Tabs)
    #show_Fs('LO','13','C2+C1-','mu1', [591,  -146,  100], Tabs)    
    #exit()

    method = 'log'

    for data in input_list: 
        order, rs, mode, grid, indices = data['order'], data['rs'], data['mode'], data['grid'], data['indices']
        xsecs, warn = get_xsec(params, data, options, Tabs, method)
        if options['scale_var'] == 'ON':            
            x1, x2, x05 = str(xsecs['mu1']).ljust(18), str(xsecs['mu2']).ljust(18), str(xsecs['mu05']).ljust(18)
            result = '{rs}  {order}  {mode}   {x1}   {x2}  {x05}'.format(rs=rs, order=order.ljust(3), mode=mode.ljust(7), x1=x1, x2=x2, x05=x05)
        else:
            x1 = str(xsecs['mu1']).ljust(18)
            result = '{rs}  {order}  {mode}   {x1}'.format(rs=rs, order=order.ljust(3), mode=mode.ljust(7), x1=x1)
        print result, warn
    if len(range_memo) > 0:
        print ''
        for memo in range_memo: print memo

    exit()





