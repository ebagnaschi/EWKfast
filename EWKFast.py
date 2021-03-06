#! /usr/bin/env python

import os, sys
import numpy as np
from math import *
from pyslha import * 
from functions import *
from scipy import interpolate
from scipy.interpolate import interp1d, interp2d, RegularGridInterpolator

if __name__ == "__main__":

    ninp = len(sys.argv)
    if ninp not in [2, 3]:
        print '[SLHA file]'
        print '[SLHA file] [input file]'        
        exit()

    slha_file = sys.argv[1] 
    input_file = 'input_default.dat'
    if ninp == 3: input_file = sys.argv[2]

    dir_path = os.getcwd()
    slha_path = os.path.join(dir_path, slha_file)
    input_path = os.path.join(dir_path, input_file)

    if not os.path.join(slha_path):
        print 'SLHA file not found:', slha_path
    if not os.path.join(input_path):
        print 'Input file not found:', input_path


    params = get_params(slha_path)
    params, range_memo = ranges(params)

    input_list, options = process_input(input_path)

    Tabs = load_tables(input_list, options) 

    #show_Fs('LO','13','NC+','mu1', [0,591,391], Tabs)
    #show_Fs('LO','13','C2+C1-','mu1', [591,  -146,  100], Tabs)    
    #exit()

    results = get_results(input_list, params, options, Tabs, method='linear')
    if options['Sort'] in ['ON', 'on', 'On']: results = sort_results(results)

    fout = open('output.ewk', 'w')

    extra = ''
    out_string = 'Energy  Precision  PDF  Process' + extra + '\n'  

    for key, data in results.items():
        order, rs, pdf, mode = key
        xsec, warn = data['xsec'], data['warn']
        x1 = str(xsec).ljust(18)
        result = '{rs}      {order}   {pdf}        {mode}    {x1}'.format(rs=rs, order=order.ljust(3), pdf=pdf.ljust(10), mode=mode.ljust(7), x1=x1)
        out_string += result +' '+ warn + '\n'
    out_string += '\n'

    out_string += 'SLHA file: ' + slha_path + '\n'
    out_string += 'Input file: ' + input_path + '\n'    
    for key, item in options.items():
        if item != '':
            out_string += '{key}: {item} \n'.format(key=key, item=item)

    if len(range_memo) > 0:
        out_string += '\n'
        for memo in range_memo: out_string += memo + '\n'


    print out_string
    fout.write(out_string)

    exit()





