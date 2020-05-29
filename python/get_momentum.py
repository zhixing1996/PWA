#!/usr/bin/env python
"""
Get four momentum of samples
"""

__author__ = "Maoqiang JING <jingmq@ihep.ac.cn>"
__copyright__ = "Copyright (c) Maoqiang JING"
__created__ = "[2020-05-29 Fri 11:19]"

import sys, os
import logging
from math import *
from ROOT import *
from tools import convert_name
logging.basicConfig(level=logging.DEBUG, format=' %(asctime)s - %(levelname)s- %(message)s')

def usage():
    sys.stdout.write('''
NAME
    get_momentum.py

SYNOPSIS
    ./get_momentum.py [ecms] [D_sample] [D_type] [sample_type]

AUTHOR
    Maoqiang JING <jingmq@ihep.ac.cn>

DATE
    May 2020
\n''')

def extract(path):
    try:
        f = TFile(path[0])
        t = f.Get('save')
        entries = t.GetEntries()
        logging.info('Entries :'+str(entries))
    except:
        logging.error(path[0] +' is invalid!')
        sys.exit()

    f_out = open(path[1], 'w')
    for i in range(entries):
        t.GetEntry(i)
        out = str(t.m_p4_Dbar[0]) + ' ' + str(t.m_p4_Dbar[1]) + ' ' + str(t.m_p4_Dbar[2]) + ' ' + str(t.m_p4_Dbar[3]) + '\n'
        f_out.write(out)
        out = str(t.m_p4_D[0]) + ' ' + str(t.m_p4_D[1]) + ' ' + str(t.m_p4_D[2]) + ' ' + str(t.m_p4_D[3]) + '\n'
        f_out.write(out)
        out = str(t.m_p4_pi[0]) + ' ' + str(t.m_p4_pi[1]) + ' ' + str(t.m_p4_pi[2]) + ' ' + str(t.m_p4_pi[3]) + '\n'
        f_out.write(out)
    f_out.close()

if __name__ == '__main__':
    args = sys.argv[1:]
    if len(args)<4:
        usage()
        sys.exit()
    ecms = int(args[0])
    D_sample = args[1]
    D_type = args[2]
    sample_type = args[3]

    if not os.path.exists('../data/'):
        os.makedirs('../data/')
    path = []
    if sample_type == 'data':
        path.append('/besfs/groups/cal/dedx/jingmq/bes/DstDpi_momentum/python/anaroot/data/data_'+str(ecms)+'_DstDpi_'+D_sample+'_'+convert_name(D_sample, D_type)+'_signal_angle_before.root')
        path.append('../data/PWA_dt.dat')
    if sample_type == 'background':
        path.append('/besfs/groups/cal/dedx/jingmq/bes/DstDpi_momentum/python/anaroot/data/data_'+str(ecms)+'_DstDpi_'+D_sample+'_'+convert_name(D_sample, D_type)+'_background_angle_before.root')
        path.append('../data/bg1.dat')
    if sample_type == 'mc':
        path.append('/besfs/groups/cal/dedx/jingmq/bes/DstDpi_momentum/python/anaroot/mc/mc_'+str(ecms)+'_DstDpi_'+D_sample+'_'+convert_name(D_sample, D_type)+'_angle_before.root')
        path.append('../data/PWA_mc.dat')
    extract(path)
