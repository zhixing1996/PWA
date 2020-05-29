#!/usr/bin/env python
import os
import sys 
import math
import random
import linecache
from ROOT import * 

def PWAStep1_MCIntegration():
    '''To do MC Integraton.
    istep=1 ==> calculate fun(i,j) using MC sample '''
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
     istep=2 ==> optimize fractions of resonances '''
    My_DynamicOptimization(nopt = 500, deltas = 0.5, ngood = 20 )

def PWAStep3_OutputMCWeight():
    '''To do PWA.
     istep=3 ==> output weight of MC sample    '''
    My_ParaFix2Best()
    # MC 
    f=open("temp","w")
    f.write("0"+"\n")
    f.close()
    os.system("mv temp  ./start/step.in")
    os.system("cp ./input/dat/num_mc.dat ./input/dat/num.dat")
    os.system("$HOME/bin/exe/for pwa")
    os.system("./pwa.fexe")
    os.system("mv ./output/out/res.out ./output/out/res_mc.out")
    f=open("temp","w")
    f.write("3"+"\n")
    f.close()
    os.system("mv temp  ./start/step.in")
    os.system("./pwa.fexe")
    os.system("mv ./output/out/reweight.wt ./output/out/mc.wt")

def My_DynamicOptimization(nopt = 1000,deltas = 0.5, ngood = 5 ): 
    '''To do optimization DYNAMICLY.
    To get the minimum S value(maximum likelihood) by changing para.in randomly.'''
    f=open("temp","w")
    f.write("2"+"\n")
    f.close()
    os.system("mv temp  ./start/step.in")
    os.system("cp ./input/dat/num_mc.dat ./input/dat/num.dat")
    os.system("$HOME/bin/exe/for pwa")
    tempfile = "./output/out/Result_opt.txt"
    if (os.path.exists(tempfile)):
       os.system("rm -f %s " % tempfile )   
       file_result=open(tempfile,"w")
    else: 
       file_result=open(tempfile,"w")
    file_result.write("fit result of dynamic optimization \n")
    file_result.close()
    mini_i    = 0 
    mini_like = 10000 
    igood = 0 
    for i in range(nopt): 
        if (igood >= ngood ) : 
           break 
        My_ParaUpdate()
        os.system("./pwa.fexe")
        like = My_FindLikelihood("./output/out/like.out")
        if (like < mini_like) : 
           os.system("cp -f ./output/out/like.out ./output/out/like.best")
           os.system("cp -f ./output/out/para.out ./output/out/para.best")
           os.system("cp -f ./output/out/res.out ./output/out/res.best")
           if (abs(like-mini_like) <= deltas) : 
              igood = igood + 1  
           else: 
              igood = 0 
           mini_i = i 
           mini_like = like 
        else:
           if (abs(like-mini_like) <= deltas) : 
              igood = igood + 1  
        file_result_line = str('%d'%i).rjust(11)+str(" : ") + str('%f'%like).rjust(24) + str(" :   igood = %d"%igood).rjust(4) + str(" :   best_like = %f"%mini_like).rjust(13) + "\n"
        file_result=open(tempfile,"a+")
        file_result.write(file_result_line)
        file_result.close()
    file_result_line = str('mini_i = %d'%mini_i).rjust(11)+str(" : ") + str('mini_like = %f'%mini_like).rjust(24) + str(" :   igood = %d"%igood).rjust(4) + "\n"
    file_result=open(tempfile,"a+")
    file_result.write(file_result_line)
    file_result.close()

def My_ParaUpdate(): 
    '''To renew para.inp.'''

    file_para = open("./input/res/para.in")
    file_para_new = open("./input/res/para.in.new","w")
    for para_line in file_para.readlines():
	temp = para_line.split()
	if ( float(temp[2]) <= 0.0 ) : 
	   para_line_new = para_line   
	if ( float(temp[2]) > 0.0 ) :   
           temp1 = para_line.split(temp[2])
	   pos3 = para_line.find(" r ")
           if (pos3 >= 0) : 
              para_value = random.uniform(0,10)
           pos4 = para_line.find(" i ")
           if (pos4 >= 0) : 
              para_value = random.uniform(-4.0, 4.0) 
           para_line_new = temp[0].rjust(4," ")+str('%4.2f'%para_value).rjust(10," ")+"  "+"  0.200"+temp1[1]
           #print "para_line    =",para_line 
	   #print "para_line_new=",para_line_new
        file_para_new.write(para_line_new)

    file_para.close()
    file_para_new.close()

    os.system("mv input/res/para.in.new  input/res/para.in") 


##############################################################
def My_ParaFix2Best(): 
    '''To fix input/res/para.in to the BEST para file. '''

    file = open("./output/out/para.best")
    file_new = open("./output/out/para.new","w")
    for line in file.readlines():
	temp = line.split()
	if ( float(temp[2]) <= 0.0 ) : 
	   line_new = line   
	if ( float(temp[2]) > 0.0 ) :   
           temp1 = line.split(" "+temp[2]+" ")
	   temp2 = list(temp1[0])
	   line_new = "".join(temp2[:14])+"0.000".rjust(9," ")+" "+temp1[1]
           #print "para_line    =",line 
	   #print "para_line_new=",line_new
        file_new.write(line_new)
    file.close()   
    file_new.close()   

    os.system("cp output/out/para.new input/res/para.in")
##################################################
def My_FindLikelihood(filename="./output/out/like.out"): 
    '''To find likelihood in like.out.'''

    file_like=open(filename,"r")
    for line in file_like.readlines():
        pos = line.find("likelihood") 
        if (pos >= 0 ):
           temp = line.split("*")
           temp1 = temp[1].split(":")
    result = float(temp1[1])
    return result

##############################################################
def My_Optimization(nopt=100): 
    '''To do optimization.
    To get the minimum S value(maximum likelihood) by changing para.in randomly.'''

    f=open("temp","w")
    f.write("2"+"\n")
    f.close()
    os.system("mv temp  ./start/step.in")
    os.system("cp ./input/dat/num_mc.dat ./input/dat/num.dat")
    os.system("$HOME/bin/exe/for pwa")
 
    file_result = open("./output/out/Result_opt.txt","w")
    
    mini_like = 1000 
    igood = 0 
    for i in range(nopt): 
        #os.system("./pwa.fexe")
        like = My_FindLikelihood()
        if (like < mini_like) : 
           mini_i = i 
           mini_like = like 
        #file_result_line = str('mini_i = %d'%mini_i).rjust(11)+str(" : ") + str('mini_like = %f'%mini_like).rjust(24) + str(" :   igood = %d"%igood).rjust(4) + "\n"
        file_result_line = str('%d'%i).rjust(11)+str(" : ") + str('%f'%like).rjust(24) + str(" :   igood = %d"%igood).rjust(4) + str(" :   best_like = %f"%mini_like).rjust(13) + "\n"
        file_result.write(file_result_line)

    file_result_line = str('mini_i = %d'%mini_i).rjust(11)+str(" : ") + str('mini_like = %f'%mini_like).rjust(24) + str(" :   igood = %d"%igood).rjust(4) + "\n"
    file_result.write(file_result_line)
    file_result.close()
##########################################################################################
def My_Result_Eff(): 
    '''To get the efficiency. '''

    file_mc  = open("./output/out/res_mc.out")
    file_mct = open("./output/out/res_mct.out")
    file_like = open("./output/out/like.out")
    file_Eff = open("./output/out/Result_Eff.out","w")

    line_Eff = str("name").rjust(21," ") + "&" + str("frac(\%)").rjust(9," ") + "&" + str("number").rjust(13," ") + "&" + str("eff.(\%)").rjust(10," ") + " \\\\ \\hline \n"
    file_Eff.write(line_Eff)

    line_lik_list = file_like.readlines()  
    line_mct_list = file_mct.readlines()  
 
    nline = len(line_mct_list)    
    iline = 0  
    for line_mc in file_mc.readlines()[:nline-1]: 
        line_mct = line_mct_list[iline]
        temp_mc = line_mc.split(":")
        temp_mct = line_mct.split(":")
        iline = iline + 1 
        if (int(temp_mc[1])==1):
           #print "temp_mc = ",temp_mc[0]
           #print "temp_mct = ",temp_mct[0]
           eff = float(temp_mc[2])/float(temp_mct[2])/100*10000 
           #print "eff = ",eff
           for line_like in line_lik_list : 
               if temp_mc[0] in line_like: 
                  temp_like = line_like.split(":")
                  name  = temp_like[0] 
                  frac  = temp_like[1] 
                  numb  = temp_like[2] 
 
           line_Eff = name + "&" + frac + "&" + numb.replace("\n","  ") + "&" + str('%2.2f'%eff).rjust(10," ") + " \\\\ \\hline \n"
           file_Eff.write(line_Eff)
    file_Eff.close()
    file_like.close()
    file_mct.close()
    file_mc.close()

#######################################################################
def My_ResUpdate(name,mass,width): 
    '''To renew res.inp.'''

    file_res = open("./input/res/res.in")
    file_res_new = open("./input/res/res.in.new","w")
    for res_line in file_res.readlines():
	res_line_new = res_line  
	pos = res_line.find(name)
	if ( pos >= 0 ) : 
	   temp = res_line.split("|")
           temp1 = temp[0].split()
           res_line_new = temp1[0].rjust(5," ")+str('%6.5f'%mass).rjust(13," ")+str('%6.5f'%width).rjust(13," ")+" " + "|"+temp[1] + "|"+temp[2]
           #print "res_line    =",res_line 
	   #print "res_line_new=",res_line_new
        file_res_new.write(res_line_new)

    file_res.close()
    file_res_new.close()

    os.system("mv input/res/res.in.new  input/res/res.in") 
#################################################################
def My_GetResultOpt(filename="output/out/Result_opt.txt"): 
    '''To get the mini_i and mini_like from the DynamicFit.'''

    file = open("./output/out/Result_opt.txt")
    for line in file.readlines():
        pos0 = line.find("mini")
        if pos0>=0 : 
           temp = line.split(":")
           temp1 = temp[0].split("=") 
           temp2 = temp[1].split("=") 
           print "No. = ",temp1[1]," like=",temp2[1]
    file.close()
    result = [temp1[1],temp2[1]]    
    return result

#################################################################

def My_GetMWBest(filename,npoint=7):
    c1=TCanvas("c1","",0,0,500,500);

    epsfile=filename+".eps"  
    pdffile=filename+".pdf"  
    filename=filename+".txt"
    file=open(filename);
    sn_graph=TGraph(npoint)
  
    ipoint = 0 
    for file_line in file.readlines(): 
	temp = file_line.split(",")
	posM = filename.find("M")
	posW = filename.find("W")
	if (posM >= 0) :
           temp1 = temp[1].split("=")
	if (posW >= 0):
           temp1 = temp[2].split("=")
  
	x = temp1[1]
        temp2 = temp[3].split("=")
	y = temp2[1]
	if ipoint==0 : 
	   xmin = x
	if ipoint==npoint-1 : 
	   xmax = x

        sn_graph.SetPoint(ipoint, float(x), float(y))

        ipoint = ipoint + 1 
  
    sn_graph.Draw("A*")  
    c1.Update(); 
  
    f1 = TF1("f1","[0]*x*x*x+[1]*x*x+[2]*x+[3]", float(xmin), float(xmax)) 
    sn_graph.Fit("f1","R")
    f1.Draw("samel")
    c1.Update(); 
    #c1.Print(epsfile)
    #c1.Print(pdffile,"pdf")
    #option = raw_input('Enter an "q" to exit, "enter" to continue : ')
  
    fitted = sn_graph.GetFunction("f1")
    #fitted.Draw()
    #c1.Update()
  
    xbest= fitted.GetMinimumX(float(xmin), float(xmax))
    ybest= fitted.GetMinimum(float(xmin), float(xmax))
    result=[xbest,ybest]

    print "X =",xbest
    print "Y =",ybest
    print "result",result

    xerror1= fitted.GetX(ybest+0.5,float(xbest), float(xmax))-xbest
    xerror2= xbest-fitted.GetX(ybest+0.5,float(xmin),float(xbest))
    xerror = (xerror1+xerror2)/2.0
    print "X + error =",xbest,"+",xerror1
    print "X - error =",xbest,"-",xerror2
    print "X +- error=",xbest,"+-",xerror,"(averaged)"

    pt = TPaveText(0.3,0.8,0.7,0.9,"brNDC");
    pt.SetFillColor(18);
    pt.SetTextAlign(12);
    pt.AddText("%.5f +- %.5f GeV"%(xbest,xerror));
    pt.Draw();
    c1.Update();
    c1.Print(epsfile,"eps")
    c1.Print(pdffile,"pdf")
    #option = raw_input('Input any key to continue: ')

    return result

#################################################################  
#################################################################  
##################################################
def My_ParaFineTune(name="K*( 892)c"): 
    '''To fine tune para.in.'''

    file_para = open("./input/res/para.in.best")
    temp = file_para.readlines()
    nline = len(temp)
    file_para.close() 
 
    file_para = open("./input/res/para.in")
    file_para_new = open("./input/res/para.in.new","w")
    for para_line in file_para.readlines():
	temp = para_line.split()
	if ( float(temp[2]) <= 0.0 ) : 
	   para_line_new = para_line   
	if ( float(temp[2]) > 0.0 ) :   
           temp1 = para_line.split(temp[2])
	   pos3 = para_line.find(" r ")
           if (pos3 >= 0) : 
              para_value = random.uniform(0,10)
           pos4 = para_line.find(" i ")
           if (pos4 >= 0) : 
              para_value = random.uniform(-4.0, 4.0) 
           para_line_new = temp[0].rjust(4," ")+str('%4.2f'%para_value).rjust(10," ")+"  "+"  0.200"+temp1[1]
           #print "para_line    =",para_line 
	   #print "para_line_new=",para_line_new
        file_para_new.write(para_line_new)

    file_para.close()
    file_para_new.close()

    os.system("mv input/res/para.in.new  input/res/para.in") 


##################################################################  
###################################################
def My_NumericalResult(): 
    '''To tablize the numberical results'''

    temp_file = open("./output/out/res.best")
    temp_line = temp_file.readlines()

    res_no=[0,0,0,0,0,0,0,0,0,0,0,0]
    res_name=[0,0,0,0,0,0,0,0,0,0,0,0]
#    res_name=()
    nres=0
    for iline in range(len(temp_line)):
	 if (":  1 :" in temp_line[iline]):
	    #print temp_line[iline]
	    temp1= (temp_line[iline]).split(":")
	    temp2= (temp1[0]).split()
	    res_no[nres]  =int((temp1[0])[0:5])
	    res_name[nres]=str((temp1[0])[6:17])
#	    print res_no[nres], res_name[nres]
	    nres=nres+1
    temp_file.close()
#    print "nres=",nres
  
    row = ""
    for ires in range(nres):
	row = row + "&" + str(res_name[ires]) 
    row=row.replace("\n","")+"\\\\ \\hline"
    print row     


    ntot=nres*nres/2+nres
    frac=[0]*ntot
    file_best = open("./output/out/like.BEST")
    i=0
    for line in file_best.readlines() : 
	 temp1=line.split()
	 frac[i]=temp1[2] 
	 i=i+1

    i=0
    for ires in range(nres) : 
	 row = str(res_name[ires])
	 for jres in range(nres) :
	     if (jres < ires) : 
	         row = row + "& 0 "
             if (jres >= ires) :
	         row= row + "&" + str(frac[i]) 
                 i=i+1
         row=row.replace("\n","")+"\\\\ \\hline"
         print row     
		 
