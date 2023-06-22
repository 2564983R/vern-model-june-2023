#Written by: Dr Rea Laila Antoniou Kourounioti

import numpy as np
import parameters as p
import model as m
from scipy.integrate import odeint
import math

def run_model(param,tiTi,TT):

	if not(testBounds(param)):
		raise ValueError("Parameter values unacceptable - outside bounds")
		
	
	p.s1=param[0]#positive
	p.s2=param[1]#positive
	p.s3=param[2]#positive
	p.sel=param[3]#0to1
	p.r1=param[4]#positive and >p.s1
	p.r2=param[5]#positive
	p.k=param[6]#positive
	p.n0i=param[7]#0to1
	p.Tq1=param[8]#-10to30
	p.Tq2=param[9]#-10to30 and >Tq1
	p.nT1=param[10]#-10to30
	p.nT2=param[11]#-10to30 and >nT1
	p.kFLC=param[12]#0to1
	p.f=param[13]#rate of FLC prod.

	p.init=[p.wi, 0, 0, (1-p.n0i)*p.r1/(p.s1+p.r1), (1-p.n0i)*p.s1/(p.s1+p.r1), 0, 0, p.n0i, 0, 0, 0, 0, 0, 1, (((1-p.n0i)*p.r1/(p.s1+p.r1)*p.f)/p.l)] 
	#print(p.init) #predicted values for the first time set
	#TT=np.arange(tiTi[0][0], tiTi[0][-1], 0.25) #from 0, time [0] till last number [-1], by steps of 0.25
	y=odeint(m.model,p.init,TT,args=(tiTi,),hmax=0.1,mxstep=500000)
	if np.amin(y)<-0.00001:
		print('Problem with: ')
		print('Params: ')
		print(p.s1,p.s2,p.s3,p.sel,p.r1,p.r2,p.k,p.n0i,p.Tq1,p.Tq2,p.nT1,p.nT2,p.kFLC)
		print('Negative values in y')
		raise ValueError('Negative values in y')
	if p.tEnd<(TT[len(TT)-1]-1):
		print('Problem with: ')
		print('Params: ')
		print(p.s1,p.s2,p.s3,p.sel,p.r1,p.r2,p.k,p.n0i,p.Tq1,p.Tq2,p.nT1,p.nT2,p.kFLC)
		print('not reached end. Stopped at t='+str(p.tEnd)+' of '+str(TT[len(TT)-1]))
		return math.inf

	(VIN3,FLC)=m.output(y)
	
	return VIN3, FLC, TT, y
#table of variables at times




def testBounds(param):
	s1=param[0]#positive
	s2=param[1]#positive
	s3=param[2]#positive
	sel=param[3]#0to1
	r1=param[4]#positive and >p.s1
	r2=param[5]#positive
	k=param[6]#positive
	n0i=param[7]#0to1
	Tq1=param[8]#-10to30
	Tq2=param[9]#-10to30 and >Tq1
	nT1=param[10]#-10to30
	nT2=param[11]#-10to30 and >nT1
	kFLC=param[12]#0to1
	
	if (s1<0)|(s2<0)|(s3<0)|(sel<0)|(sel>1)|(r1<0)|(r1<s1)|(r2<0)|(k<0)|(n0i<0)|(n0i>1)|(Tq1<-10)|(Tq2<Tq1)|(Tq2>30)|(nT1<-10)|(nT2<nT1)|(nT2>30)|(kFLC<0)|(kFLC>1):
		#reject
		return bool(0)
	else:
		#accept
		return bool(1)
