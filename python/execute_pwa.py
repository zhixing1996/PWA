#!/bin/env python
"""
Execute PWA step by step
"""

import os
import sys 
import math
import logging
from pwa_module import *
__author__ = "Maoqiang JING <jingmq@ihep.ac.cn>"
__copyright__ = "Copyright (c) Maoqiang JING"
__created__ = "[2020-05-29 Fri 17:50]"

logging.basicConfig(level=logging.DEBUG, format=' %(asctime)s - %(levelname)s- %(message)s')

def usage():
    sys.stdout.write('''
NAME
    execute_pwa.py

SYNOPSIS
    ./execute_pwa.py [step]

AUTHOR
    Maoqiang JING <jingmq@ihep.ac.cn>

DATE
    May 2020
\n''')

if __name__ == '__main__':
    args = sys.argv[1:]
    if len(args)<1:
        usage()
        sys.exit()
    step = args[0]
    if step == 'step1': 
        PWAStep1_MCIntegration()
    if step == 'step2': 
        PWAStep2_MaxLikelihoodFit()
    if step == 'step3': 
        PWAStep3_OutputMCWeight() 
    if step == 'step4': 
        PWA_NumericalResult()
