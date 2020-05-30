#!/usr/bin/env python
"""
Convert pwa results to hist
"""

__author__ = "Maoqiang JING <jingmq@ihep.ac.cn>"
__copyright__ = "Copyright (c) Maoqiang JING"
__created__ = "[2020-05-30 Sat 09:49]"

import sys, os
import logging
from math import *
from ROOT import *
from tools import convert_name
logging.basicConfig(level=logging.DEBUG, format=' %(asctime)s - %(levelname)s- %(message)s')

def usage():
    sys.stdout.write('''
NAME
    pwa_convert.py

SYNOPSIS
    ./pwa_convert.py [ecms] [D_sample] [D_type] [sample_type]

AUTHOR
    Maoqiang JING <jingmq@ihep.ac.cn>

DATE
    May 2020
\n''')

def set_canvas_style(mbc):
    mbc.SetFillColor(0)
    mbc.SetLeftMargin(0.15)
    mbc.SetRightMargin(0.15)
    mbc.SetTopMargin(0.1)
    mbc.SetBottomMargin(0.15)

def convert(path, sample_type):
    try:
        f = TFile(path[0])
        t = f.Get('save')
        entries = t.GetEntries()
        logging.info('Entries :'+str(entries))
    except:
        logging.error(path[0] +' is invalid!')
        sys.exit()

    h_mDpi = TH1F('h_mDpi', 'h_mDpi', 160, 2.0, 2.8)
    h_mDbarpi = TH1F('h_mDbarpi', 'h_mDbarpi', 160, 2.0, 2.8)
    h_mDDbar = TH1F('h_mDDbar', 'h_mDDbar', 180, 3.7, 4.6)
    h_dalitz = TH2F('h_dalitz', 'h_dalitz', 30, 3.5, 8.5, 30, 3.5, 8.5)
    f_wt = open('../output/out/mc.wt', 'r')
    wt_lines = f_wt.readlines()
    for i in range(entries):
        t.GetEntry(i)
        if not sample_type == 'mc':
            wt_sig = 1.
        if sample_type == 'mc':
            temp = wt_lines[i].split()
            wt_sig = float(temp[0])
            wt_bkg = float(temp[1])
        h_mDpi.Fill(t.m_m_Dpi, wt_sig)
        h_mDbarpi.Fill(t.m_rm_D, wt_sig)
        h_mDDbar.Fill(t.m_rm_pi, wt_sig)
        h_dalitz.Fill(t.m_m_Dpi*t.m_m_Dpi, t.m_rm_D*t.m_rm_D, wt_sig)
    f_wt.close()
    f_out = TFile(path[1], 'recreate')
    h_mDpi.Write()
    h_mDbarpi.Write()
    h_mDDbar.Write()
    h_dalitz.Write()
    f_out.Close()

if __name__ == '__main__':
    args = sys.argv[1:]
    if len(args)<4:
        usage()
        sys.exit()
    ecms = int(args[0])
    D_sample = args[1]
    D_type = args[2]
    sample_type = args[3]

    if not os.path.exists('./anaroot/'):
        os.makedirs('./anaroot/')
    path = []
    if sample_type == 'data':
        path.append('/besfs/groups/cal/dedx/jingmq/bes/DstDpi_momentum/python/anaroot/data/data_'+str(ecms)+'_DstDpi_'+D_sample+'_'+convert_name(D_sample, D_type)+'_signal_angle_before.root')
        path.append('./anaroot/PWA_dt.root')
    if sample_type == 'background':
        path.append('/besfs/groups/cal/dedx/jingmq/bes/DstDpi_momentum/python/anaroot/data/data_'+str(ecms)+'_DstDpi_'+D_sample+'_'+convert_name(D_sample, D_type)+'_background_angle_before.root')
        path.append('./anaroot/bg1.root')
    if sample_type == 'mc':
        path.append('/besfs/groups/cal/dedx/jingmq/bes/DstDpi_momentum/python/anaroot/mc/mc_'+str(ecms)+'_DstDpi_'+D_sample+'_'+convert_name(D_sample, D_type)+'_angle_before.root')
        path.append('./anaroot/PWA_mc.root')
    convert(path, sample_type)
