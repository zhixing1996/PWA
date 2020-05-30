#!/usr/bin/env python
import os
import sys 
import math
import random
import linecache
from ROOT import * 

def PWAStep1_MCIntegration():
    '''To do MC Integraton.
    istep=1 ==> calculate fun(i,j) using MC sample. '''
    f=open('temp', 'w')
    f.write('1' + '\n')
    f.close()
    os.system('mkdir -p ../start')
    os.system('mkdir -p ../output/out')
    os.system('mv temp ../start/step.in')
    os.system('cp ../input/init/para.in ../input/res/para.in')
    os.system('cp ../input/init/res.in ../input/res/res.in')
    os.system('cd .. && ./bin/exe pwa')
    os.system('cd .. && ./pwa.fexe')

def PWAStep2_MaxLikelihoodFit():
    '''To do maximum likelihood fit(fumili) and get the S(-log(likelihood)) value.
     istep=2 ==> optimize fractions of resonances. '''
    PWA_DynamicOptimization(nopt = 500, deltas = 0.5, ngood = 20)

def PWAStep3_OutputMCWeight():
    '''To do PWA.
     istep=3 ==> output weight of MC sample. '''
    PWA_ParaFix2Best()
    # MC 
    f=open('temp', 'w')
    f.write('0' + '\n')
    f.close()
    os.system('mv temp ../start/step.in')
    os.system('cd .. && ./bin/exe pwa')
    os.system('cd .. && ./pwa.fexe')
    os.system('mv ../output/out/res.out ../output/out/res_mc.out')
    f = open('temp', 'w')
    f.write('3' + '\n')
    f.close()
    os.system('mv temp ../start/step.in')
    os.system('cd .. && ./pwa.fexe')
    os.system('mv ../output/out/reweight.wt ../output/out/mc.wt')

def PWA_DynamicOptimization(nopt = 1000, deltas = 0.5, ngood = 5): 
    '''To do optimization DYNAMICLY.
    To get the minimum S value(maximum likelihood) by changing para.in randomly. '''
    f = open('temp', 'w')
    f.write('2' + '\n')
    f.close()
    os.system('mv temp ../start/step.in')
    os.system('cd .. && ./bin/exe pwa')
    tempfile = '../output/out/Result_opt.txt'
    if (os.path.exists(tempfile)):
        os.system('rm -f %s ' % tempfile )   
        file_result = open(tempfile, 'w')
    else: 
        file_result = open(tempfile, 'w')
    file_result.write('fit result of dynamic optimization \n')
    file_result.close()
    mini_i = 0 
    mini_like = 10000 
    igood = 0 
    for i in range(nopt): 
        if (igood >= ngood): 
            break 
        PWA_ParaUpdate()
        os.system('cd .. && ./pwa.fexe')
        like = PWA_FindLikelihood('../output/out/like.out')
        if (like < mini_like): 
            os.system('cp -f ../output/out/like.out ../output/out/like.best')
            os.system('cp -f ../output/out/para.out ../output/out/para.best')
            os.system('cp -f ../output/out/res.out ../output/out/res.best')
            if (abs(like - mini_like) <= deltas):
                igood += 1  
            else: 
                igood = 0 
            mini_i = i 
            mini_like = like 
        else:
            if (abs(like - mini_like) <= deltas):
              igood += 1  
        file_result_line = str('%d'%i).rjust(11) + str(' : ') + str('%f'%like).rjust(24) + str(' :   igood = %d'%igood).rjust(4) + str(' :   best_like = %f'%mini_like).rjust(13) + '\n'
        file_result = open(tempfile, 'a+')
        file_result.write(file_result_line)
        file_result.close()
    file_result_line = str('mini_i = %d'%mini_i).rjust(11) + str(' : ') + str('mini_like = %f'%mini_like).rjust(24) + str(' :   igood = %d'%igood).rjust(4) + '\n'
    file_result = open(tempfile, 'a+')
    file_result.write(file_result_line)
    file_result.close()

def PWA_ParaUpdate(): 
    '''To renew para.inp. '''
    file_para = open('../input/res/para.in')
    file_para_new = open('../input/res/para.in.new', 'w')
    for para_line in file_para.readlines():
        temp = para_line.split()
        if (float(temp[2]) <= 0.0): 
            para_line_new = para_line   
        if (float(temp[2]) > 0.0):   
            temp1 = para_line.split(temp[2])
            pos3 = para_line.find(' r ')
            if (pos3 >= 0): 
                para_value = random.uniform(0, 10)
            pos4 = para_line.find(' i ')
            if (pos4 >= 0): 
                para_value = random.uniform(-4.0, 4.0) 
            para_line_new = temp[0].rjust(4, ' ') + str('%4.2f'%para_value).rjust(10, ' ') + '  ' + '  0.200' + temp1[1]
        file_para_new.write(para_line_new)
    file_para.close()
    file_para_new.close()
    os.system('mv ../input/res/para.in.new  ../input/res/para.in') 

def PWA_ParaFix2Best(): 
    '''To fix input/res/para.in to the BEST para file. '''
    file = open('../output/out/para.best')
    file_new = open('../output/out/para.new', 'w')
    for line in file.readlines():
        temp = line.split()
        if (float(temp[2]) <= 0.0): 
            line_new = line   
        if (float(temp[2]) > 0.0):   
            temp1 = line.split(" "+temp[2]+" ")
            temp2 = list(temp1[0])
            line_new = ''.join(temp2[:14]) + '0.000'.rjust(9, ' ') + ' ' + temp1[1]
        file_new.write(line_new)
    file.close()   
    file_new.close()   
    os.system('cp ../output/out/para.new ../input/res/para.in')

def PWA_FindLikelihood(filename = '../output/out/like.out'): 
    '''To find likelihood in like.out. '''
    file_like = open(filename, 'r')
    for line in file_like.readlines():
        pos = line.find('likelihood') 
        if (pos >= 0):
            temp = line.split('*')
            temp1 = temp[1].split(':')
    result = float(temp1[1])
    return result

def PWA_NumericalResult(): 
    '''To tablize the numberical results. '''
    temp_file = open('../output/out/res.best')
    temp_line = temp_file.readlines()
    res_no = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    res_name = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    nres = 0
    for iline in range(len(temp_line)):
        if (':  1 :' in temp_line[iline]):
            temp1 = (temp_line[iline]).split(':')
            temp2 = (temp1[0]).split()
            res_no[nres] = int((temp1[0])[0:5])
            res_name[nres] = str((temp1[0])[6:17])
            nres += 1
    temp_file.close()
    f_out = open('../output/out/coef.out', 'w')
    row = ''
    for ires in range(nres):
        if not ires == nres - 1:
            row += str(res_name[ires]) + '& '
        if ires == nres - 1:
            row += str(res_name[ires])
    row = row.replace('\n', '') + '\\\\ \\hline\n'
    f_out.write(row)
    ntot = nres*nres/2 + nres
    frac = [0]*ntot
    PWA_ProcessFile('../output/out/like.best')
    file_BEST = open('../output/out/like.BEST')
    i = 0
    for line in file_BEST.readlines(): 
        temp1=line.split()
        frac[i]=temp1[2] 
        i += 1
    i = 0
    for ires in range(nres): 
        row = str(res_name[ires])
        for jres in range(nres):
            if (jres < ires):
                row += '& 0.000 '
            if (jres >= ires):
                row += '& ' + str(frac[i]) + ' '
                i=i+1
        row=row.replace('\n', '') + '\\\\ \\hline\n'
        f_out.write(row)
    f_out.close()

def PWA_ProcessFile(filename = '../output/out/like.best'): 
    '''Process file: ../output/out/like.best --> ../output/out/like.BEST. '''
    file_best = open(filename, 'r')
    file_BEST = open('../output/out/like.BEST', 'w')
    for line in file_best.readlines():
        if line.find('*') > 0:
            if line.find('*(') < 0:
                continue
        if line.find('events') < 0 and line.find('=') < 0 and line.find(':') < 0 and line.find('sum') < 0 and not len(line) == 1 and not len(line) == 2:
            line_new = line
            file_BEST.write(line_new)
    file_BEST.close()
