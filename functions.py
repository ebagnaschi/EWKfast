#! /usr/bin/env python

import os, sys
import numpy as np
from math import *
from scipy import interpolate
from pyslha import * 
from scipy.interpolate import interp1d, interp2d, RegularGridInterpolator
from collections import OrderedDict
from re import findall 

def get_line(List):
    line = str(List[0])
    for i in range(1, len(List)):
        line += '  ' + str(List[i])
    return line

class Get_Tables:

    def __init__(self):

        #xar = np.load( os.path.join(dirpath, grid_file + '_x') )
        #yar = np.load( os.path.join(dirpath, grid_file + '_y') )
        #zar = np.load( os.path.join(dirpath, grid_file + '_z') )
        #self.xslim_interp = RegularGridInterpolator((xar, yar), zar)                
        self.tables = {}

    def add(self, key):

        order, rs, grid, sc = key

        dir_path = os.path.dirname(os.path.realpath(__file__))

        gridname = key[2]
        if 'same' in gridname:
            dim = 2
        else:
            dim = 3

        nF = {}
        nF['CCsame'] = 6
        nF['NNsame'] = 5
        nF['NN'] = 5
        nF['NC+'] = 6
        nF['NC-'] = 6
        nF['C2+C1-'] = 3
        nF['C2-C1+'] = 3

        basetag = '{order}-{rs}_{grid}'.format(order=order, rs=rs, grid=grid)            
        tag = '{order}-{rs}_{grid}_{sc}'.format(order=order, rs=rs, grid=grid, sc=sc)            
        if dim == 2:
            fpath_mass = os.path.join(dir_path, 'lookups', basetag)
            fpath = os.path.join(dir_path, 'lookups', tag)                
            m1 = np.load( fpath_mass + '.m1' )
            m2 = np.load( fpath_mass + '.m2' )
            F_ar = []
            for ii in xrange(nF[grid]): 
                i_F = ii+1
                Far = np.load( fpath + '.F{i}'.format(i = i_F))
                F_ar.append( RegularGridInterpolator( (m1, m2), Far) )
            self.tables[key] = F_ar

        if dim == 3:
            fpath_mass = os.path.join(dir_path, 'lookups', basetag)
            fpath = os.path.join(dir_path, 'lookups', tag)                
            m1 = np.load( fpath_mass + '.m1' )
            m2 = np.load( fpath_mass + '.m2' )
            m3 = np.load( fpath_mass + '.m3' )                
            F_ar = []
            for ii in xrange(nF[grid]): 
                i_F = ii+1
                Far = np.load( fpath + '.F{i}'.format(i = i_F) )
                F_ar.append( RegularGridInterpolator( (m1, m2, m3), Far) )
            self.tables[key] = F_ar


def load_tables(input_list, options):

    if options['scale_var'] in ['ON', 'On', 'on']:
        scales = ['mu1', 'mu2', 'mu05']
    else:
        scales = ['mu1']

    Tabs = Get_Tables()

    for data in input_list: 
        rs, order, grid, = data['rs'], data['order'], data['grid']        
        for sc in scales:
            key=(order, rs, grid, sc)
            if key not in Tabs.tables.keys(): 
                #print data['mode'], data['grid']
                Tabs.add(key)
    return Tabs


def show_Fs(order,rs,grid,scale, masses, Tabs):
    table = Tabs.tables[(order,rs,grid,scale)]
    Fs = []
    for tab in table: Fs.append( tab(masses)[0] )  
    line = get_line(Fs)
    print grid, masses, line


def process_input(input_path):
    input_list = []
    options = {'scale_var': ''}
    for line in open(input_path):
        line0 = line.split('#')[0]
        ops = line0.split(':')
        if len(ops) > 1:
            if ops[0] == 'scale variation': options['scale_var'] = ops[1].strip()
            continue
        elems = line.split('#')[0].split()            
        if len(elems) != 3: continue
        data = OrderedDict()
        rs, order, mode = elems
        data['rs'] = rs 
        data['order'] = order
        data['mode'] = mode
        indices = findall(r'[1-4]+', mode)
        indices = np.array(map(int, indices))
        indices = indices - np.array([1, 1])         
        data['indices'] = indices
        data['grid'] = 'empty'
        if mode.count('N') == 2:
            if mode in ['N1N1', 'N2N2', 'N3N3', 'N4N4']: 
                data['grid'] = 'NNsame'
            else: 
                data['grid'] = 'NN'
        if mode in ['C1C1', 'C2C2']: 
            data['grid'] = 'CCsame'
        if mode == 'C1+C2-': data['grid'] = 'C2-C1+' 
        if mode == 'C1-C2+': data['grid'] = 'C2+C1-' 
        if (mode.count('N'), mode.count('C'), mode.count('+')) == (1,1,1): data['grid'] = 'NC+'
        if (mode.count('N'), mode.count('C'), mode.count('-')) == (1,1,1): data['grid'] = 'NC-'
        if data['grid'] == 'empty':
            print 'grid is empty'
            print mode, 'is not defined'
            exit()
        input_list.append( data )
    return input_list, options


def get_xsec(params, data, options, table_dic):

    rs, order, grid, indices = data['rs'], data['order'], data['grid'], data['indices']

    if options['scale_var'] in ['ON', 'On', 'on']:
        scales = ['mu1', 'mu2', 'mu05']
    else:
        scales = ['mu1']

    m1, m2 = 0, 0

    i1, i2 = indices
    mQ = params['mQL']    
    if grid == 'CCsame':
        m1 = abs(params['mC'][i1])
    if grid == 'NNsame':
        m1 = abs(params['mN'][i1])
    if grid in ['C2-C1+', 'C2+C1-']:
        m1 = abs(params['mC'][i1])
        m2 = abs(params['mC'][i2])
        m1, m2 = max(m1, m2), min(m1, m2)
    if grid in ['NN']:
        m1 = params['mN'][i1]
        m2 = params['mN'][i2]
        sgn = np.sign(m1*m2)
        m1, m2 = max(abs(m1), abs(m2)), sgn*min(abs(m1), abs(m2))
    if grid in ['NC+', 'NC-']:
        m2 = params['mC'][i2]
        m1 = params['mN'][i1] 
        m2 = np.sign(m1*m2)*m2
        m1 = abs(m1)
        print 'here', m1, m2

    warn = ''
    if m1 == 1500 - 1: warn = '# mass is shifted in the grid'
    if m2 == 1500 - 1: warn = '# mass is shifted in the grid'

    vec = get_vec(grid, params, i1, i2)

    xsecs = OrderedDict()
    for sc in scales:
        key=(order, rs, grid, sc)
        table = table_dic[key]
        Fs = []
        if grid in ['CCsame', 'NNsame']:
            for tab in table: Fs.append( tab([m1, mQ])[0] )  
        else:
            for tab in table: Fs.append( tab([m1, m2, mQ])[0] )  

        Fs = np.array(Fs)
        xsecs[sc] = np.dot(vec, Fs)
    return xsecs, warn

def get_params(SLHAfile):

    blocks, decays = readSLHAFile(SLHAfile)

    mN, mC, N, V, U = [], [], [], [], []

    mN.append( blocks["MASS"].entries[1000022] ) 
    mN.append( blocks["MASS"].entries[1000023] ) 
    mN.append( blocks["MASS"].entries[1000025] ) 
    mN.append( blocks["MASS"].entries[1000035] ) 

    mC.append( blocks["MASS"].entries[1000024] )
    mC.append( blocks["MASS"].entries[1000037] )

    mQL = 0
    mQL += blocks["MASS"].entries[1000001]
    mQL += blocks["MASS"].entries[1000002]
    mQL += blocks["MASS"].entries[1000003]
    mQL += blocks["MASS"].entries[1000004]
    mQL = mQL/4.

    mQR = 0
    mQR += blocks["MASS"].entries[2000001]
    mQR += blocks["MASS"].entries[2000002]
    mQR += blocks["MASS"].entries[2000003]
    mQR += blocks["MASS"].entries[2000004]
    mQR = mQR/4.

    for ii in xrange(4):
        i = ii + 1
        elem1 = blocks['NMIX'].entries[i][1]
        elem2 = blocks['NMIX'].entries[i][2]
        elem3 = blocks['NMIX'].entries[i][3]
        elem4 = blocks['NMIX'].entries[i][4]
        N.append([elem1, elem2, elem3, elem4])

    for ii in xrange(2):
        i = ii + 1        
        v1 = blocks['VMIX'].entries[i][1]
        v2 = blocks['VMIX'].entries[i][2]
        u1 = blocks['UMIX'].entries[i][1]
        u2 = blocks['UMIX'].entries[i][2]
        V.append([v1, v2])
        U.append([u1, u2])

    params = {}
    params['mQL'] = mQL    
    params['mQR'] = mQR    
    params['mN'] = mN
    params['mC'] = mC
    params['N'] = np.array(N)
    params['V'] = np.array(V)
    params['U'] = np.array(U)

    return params


def ranges(params):
    range_memo = []
    bound_dic = OrderedDict([('mQR', 6000), ('mQL', 6000)])
    for p, bound in bound_dic.items():
        if params[p] >= bound: 
            memo = '{p} = {val} is out of range. {p} = {bound} was used.'.format(p=p,val=params[p], bound=bound)
            range_memo.append(memo)        
            params[p] = bound - 1

    bound = 1500
    p = 'mC'
    for i in [0, 1]:
        val = params[p][i]
        if abs(val) >= bound:
            memo = '{p}{i} = {val} is out of range. |{p}{i}| = {bound} was used.'.format(p=p,val=params[p][i], bound=bound, i=i+1)
            range_memo.append(memo)        
            params[p][i] = bound - 1

    bound = 1500
    p = 'mN'
    for i in [0, 1, 2, 3]:
        val = params[p][i]
        if abs(val) >= bound:
            memo = '{p}{i} = {val} is out of range. |{p}{i}| = {bound} was used.'.format(p=p,val=params[p][i], bound=bound, i=i+1)
            range_memo.append(memo)        
            params[p][i] = bound - 1


    return params, range_memo

#=================================================#

sw = sqrt(1. - 80.410003662109375**2 / 91.186996459960938**2)
cw = sqrt(1 - sw**2)
r2 = sqrt(2.)

deno = sw * cw
Lu = ( 1./2. - 2./3. * sw**2 ) / deno
Ld = (-1./2 - (-1./3.) * sw**2 ) / deno
Ru = (-2./3. * sw**2 ) / deno
Rd = ( - (-1./3.) * sw**2 ) / deno

def Lut(i, n): 
    return (n[i, 0]*cw + n[i, 1]*sw) * (2./3.) + (-n[i, 0]*sw + n[i, 1]*cw) * Lu

def Ldt(i, n): 
    return (n[i, 0]*cw + n[i, 1]*sw) * (-1./3.) + (-n[i, 0]*sw + n[i, 1]*cw) * Ld

def Rut(i, n): 
    return (n[i, 0]*cw + n[i, 1]*sw) * (2./3.) + (-n[i, 0]*sw + n[i, 1]*cw) * Ru

def Rdt(i, n): 
    return (n[i, 0]*cw + n[i, 1]*sw) * (-1./3.) + (-n[i, 0]*sw + n[i, 1]*cw) * Rd

def gNC( i, j, n, v, u ):
    elem1 = n[i, 1] * v[j, 0] - n[i, 3] * v[j, 1]/sqrt(2.) 
    elem2 = v[j, 0] * Lut(i, n)
    elem3 = n[i, 1] * u[j, 0] + n[i, 2] * u[j, 1]/sqrt(2.) 
    elem4 = u[j, 0] * Ldt(i, n)
    return [elem1, elem2, elem3, elem4]

def cNC( v, c ):
    if c == 0: return v[0]**2 + v[2]**2
    if c == 1: return v[1]**2 + v[3]**2
    if c == 2: return 2.*( v[0]*v[1] - v[2]*v[3] )
    if c == 3: return 2.*( v[1]*v[2] - v[0]*v[3] )
    if c == 4: return 2.* v[0]*v[2]
    if c == 5: return 2.* v[1]*v[3]

def gNN( i, j, n ):
    elem1 = n[i, 3] * n[j, 3] - n[i, 2] * n[j, 2]
    elem2 = Lut(i, n) * Lut(j, n)
    elem3 = Rut(i, n) * Rut(j, n)
    elem4 = Ldt(i, n) * Ldt(j, n)
    elem5 = Rdt(i, n) * Rdt(j, n)
    return [elem1, elem2, elem3, elem4, elem5]

def cNN( v, c ):
    if c == 0: return v[0]**2
    if c == 1: return v[1]**2 + v[2]**2
    if c == 2: return 2.*( Lu * v[0] * v[1] - Ru * v[0] * v[2] )
    if c == 3: return v[3]**2 + v[4]**2
    if c == 4: return 2.*( Ld * v[0] * v[3] - Rd * v[0] * v[4] )    


def gCC( i, j, u, v ):
    elem1 = u[i, 0] * u[j, 0]
    elem2 = v[i, 0] * v[j, 0]
    elem3 = 0
    if i == j: elem3 = 1
    return [elem1, elem2, elem3]

def cCC( v, c ):
    if c == 0: return v[0]**2
    if c == 1: return 2. * v[0] * v[1]
    if c == 2: return v[1]**2
    if c == 3: return 2. * v[0] * v[2]
    if c == 4: return 2. * v[1] + v[2]
    if c == 5: return v[2]**2


def get_vec(run_mode, params, i1, i2):

    V = params['V']
    U = params['U']
    N = params['N']

    if run_mode in ['CCsame', 'C2+C1-', 'C2-C1+']:
        gcc = gCC(i1, i2, U, V)
        n = 3
        if run_mode == 'CCsame': n=6
        vec = [ cCC(gcc, i) for i in xrange(n) ]

    if run_mode in ['NN', 'NNsame']:
        gnn = gNN(i1, i2, N)
        vec = [ cNN(gnn, i) for i in xrange(5) ]

    if run_mode in ['NC+', 'NC-']:
        gnc = gNC(i1, i2, N, V, U)
        vec = [ cNC(gnc, i) for i in xrange(6) ]

    return vec

#####################################################






# def get_coeffs_2D(run_mode, masses, tableFile):

#     data = np.loadtxt(tableFile)
#     m1_ar = data[:,0]
#     m2_ar = data[:,1]
#     F_ar = []
#     for ii in range(10): 
#         F_ar.append( data[:,ii+2] )

#     coeff_seeds = []
#     for ii in range(10):      
#         tmp_list = []      
#         for j in range(len(list(m1_ar))):
#             m1 = m1_ar[j]
#             m2 = m2_ar[j]
#             val = F_ar[ii][j]    
#             tmp_list.append( [m1, m2, val ] )
#         coeff_seeds.append(tmp_list) 

#     m1_in = masses[0]
#     m2_in = masses[1]

#     Coeffs = []
#     for ii in range(10):
#         v_interp = Interpolate2D(coeff_seeds[ii], [m1_in, m2_in])
#         Coeffs.append(v_interp)

#     return np.array(Coeffs)



# def get_xsec(params, mode, indices, nlo_flag):
#     if indices[0] == indices[1]: 
#         return get_same(params, mode, indices, nlo_flag)
#     else:
#         if 'NN' in mode: return get_NN(params, mode, indices, nlo_flag)



# def get_same(params, mode, indices, nlo_flag):

#     mode = mode + 'same'

#     tableDir = './'
#     table_mid = tableDir + '{order}-{E}_{mode}.table'.format(order=nlo_flag, E=rootS, mode=mode)

#     ii = indices[0]

#     mQ = params['mQ']
#     if 'CC' in mode: 
#         m1_raw = params['mC'][ii]
#         prod = 'C'
#     if 'NN' in mode: 
#         m1_raw = params['mN'][ii]
#         prod = 'N'

#     m1 = abs(m1_raw)        
#     if m1 > 1000.: 
#         print prod + '{ii} is too heavy (|{mass}| >1TeV): skip'.format(ii=str(ii+1), mass=str(m1_raw))
#         return 
#     #Coeff05 = get_coeffs_2D(mode,[m1, mQ], table05)
#     Coeff10 = get_coeffs_2D(mode,[m1, mQ], table10)
#     #Coeff20 = get_coeffs_2D(mode,[m1, mQ], table20)

#     Mline = get_Mline(mode, params, ii, ii)
#     xsec05 = np.dot(Mline, Coeff05)
#     xsec10 = np.dot(Mline, Coeff10)
#     xsec20 = np.dot(Mline, Coeff20)

#     if mode == 'CCsame': print 'C{ii}-C{ii} (m{prod}{ii}, mQ) = ({mX}, {mQ}) {nlo_flag} [pb]'.format(ii=str(ii+1), prod=prod, mX=m1_raw, mQ=mQ, nlo_flag=nlo_flag)
#     if mode == 'NNsame': print 'N{ii}-N{ii} (m{prod}{ii}, mQ) = ({mX}, {mQ}) {nlo_flag} [pb]'.format(ii=str(ii+1), prod=prod, mX=m1_raw, mQ=mQ, nlo_flag=nlo_flag)

#     print 'scale 0.5:  ', xsec05
#     print 'scale 1.0:  ', xsec10
#     print 'scale 2.0:  ', xsec20



# def get_same_old(params, mode, indices, nlo_flag):

#     mode = mode + 'same'

#     if nlo_flag == 'NLO': 
#         tableDir = 'NLO_tables/'
#         table05 = tableDir + 'table_' + mode + '_8_NLO_05.dat'
#         table10 = tableDir + 'table_' + mode + '_8_NLO_10.dat'
#         table20 = tableDir + 'table_' + mode + '_8_NLO_20.dat'
#     if nlo_flag == 'LO': 
#         tableDir = 'table_1/'
#         table05 = tableDir + 'table_' + mode + '_8_LO_05.dat'
#         table10 = tableDir + 'table_' + mode + '_8_LO_10.dat'
#         table20 = tableDir + 'table_' + mode + '_8_LO_20.dat'

#     ii = indices[0]

#     mQ = params['mQ']
#     if 'CC' in mode: 
#         m1_raw = params['mC'][ii]
#         prod = 'C'
#     if 'NN' in mode: 
#         m1_raw = params['mN'][ii]
#         prod = 'N'

#     m1 = abs(m1_raw)        
#     if m1 > 1000.: 
#         print prod + '{ii} is too heavy (|{mass}| >1TeV): skip'.format(ii=str(ii+1), mass=str(m1_raw))
#         return 
#     Coeff05 = get_coeffs_2D(mode,[m1, mQ], table05)
#     Coeff10 = get_coeffs_2D(mode,[m1, mQ], table10)
#     Coeff20 = get_coeffs_2D(mode,[m1, mQ], table20)

#     Mline = get_Mline(mode, params, ii, ii)
#     xsec05 = np.dot(Mline, Coeff05)
#     xsec10 = np.dot(Mline, Coeff10)
#     xsec20 = np.dot(Mline, Coeff20)

#     if mode == 'CCsame': print 'C{ii}-C{ii} (m{prod}{ii}, mQ) = ({mX}, {mQ}) {nlo_flag} [pb]'.format(ii=str(ii+1), prod=prod, mX=m1_raw, mQ=mQ, nlo_flag=nlo_flag)
#     if mode == 'NNsame': print 'N{ii}-N{ii} (m{prod}{ii}, mQ) = ({mX}, {mQ}) {nlo_flag} [pb]'.format(ii=str(ii+1), prod=prod, mX=m1_raw, mQ=mQ, nlo_flag=nlo_flag)

#     print 'scale 0.5:  ', xsec05
#     print 'scale 1.0:  ', xsec10
#     print 'scale 2.0:  ', xsec20



# def get_NN(params, mode, indices, nlo_flag):

#     if nlo_flag == 'NLO': 
#         tableDir = 'NLO_tables/'
#         table05 = tableDir + 'table_' + mode + '_8_NLO_05.dat'
#         table10 = tableDir + 'table_' + mode + '_8_NLO_10.dat'
#         table20 = tableDir + 'table_' + mode + '_8_NLO_20.dat'
#     if nlo_flag == 'LO': 
#         tableDir = 'table_1/'
#         table05 = tableDir + 'table_' + mode + '_8_LO_05.dat'
#         table10 = tableDir + 'table_' + mode + '_8_LO_10.dat'
#         table20 = tableDir + 'table_' + mode + '_8_LO_20.dat'

#     i1, i2 = min(indices), max(indices)

#     mQ = params['mQ']
#     mL_raw = params['mN'][i1]
#     mH_raw = params['mN'][i2]

#     mH, mL = mH_raw, mL_raw
#     if mH_raw < 0 and mL_raw < 0: mH, mL = abs(mH_raw), abs(mL_raw)
#     if mH_raw < 0 and mL_raw > 0: mH, mL = abs(mH_raw), - mL_raw
    
#     if abs(mH) > 1500.: 
#         print prod + '{i2} is too heavy (|{mass}| >1.5TeV): skip'.format(i2=str(i2+1), mass=str(mH_raw))
#         return 
#     Coeff05 = get_coeffs_3D(mode, [mH, mL, mQ], table05)
#     Coeff10 = get_coeffs_3D(mode, [mH, mL, mQ], table10)
#     Coeff20 = get_coeffs_3D(mode, [mH, mL, mQ], table20)

#     Mline = get_Mline(mode, params, i1, i2)
#     xsec05 = np.dot(Mline, Coeff05)
#     xsec10 = np.dot(Mline, Coeff10)
#     xsec20 = np.dot(Mline, Coeff20)

#     print 'N{i1}-N{i2} (mN{i1}, mN{i2}, mQ) = ({m1}, {m2}, {mQ}) {nlo_flag} [pb]'.format(i1=i1+1, i2=i2+1, m1=mL_raw, m2=mH_raw, mQ=mQ, nlo_flag=nlo_flag)

#     print 'scale 0.5:  ', xsec05
#     print 'scale 1.0:  ', xsec10
#     print 'scale 2.0:  ', xsec20
