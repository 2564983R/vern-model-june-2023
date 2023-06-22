#Written by: Dr Rea Laila Antoniou Kourounioti

import math
import numpy as np
import csv

# parameters 
s1=			0.016 #ps1 #Rate of H to I shutdown constant (VIN3 independent)
s2=			0.0111#0.0111 #rate of I to N shutdown constant (VIN3 dependent)
s3=			0.75 #rate of S to P shutdown constant (VIN3 dependent in dividing cells only)
sel=		0.1 #ps (lowercase p)
r1=			0.05 #pr1 #reversal rate constant VIN3 independent (I to H)
r2=			0 #pr2 #Reversal rate constant VIN3 dependent (P to I (in dividing cells only))
k=			0.18 #pk #constant for decrease in current growth rate in permissive flowering conditions
n0i=		0.2 #proportion of silenced cells initially
Tq1=		-1 #temperature bound for VIN3 dependent transitions #T1
Tq2=		18 #temperature bound for VIN3 dependent transitions #T2
nT1=		11.5 #temperature bound for VIN3 independent transitions
nT2=		15 #temperature bound for VIN3 independent transitions
kFLC=       0.04 


# constant parameters
dn=			32 #constant number of non-dividing cells (based on each non-stem cell dividing 5 times)
g1=         0.4 #growth rate at warm temp (22 degrees C)
growthe=	0.22 #constant used in exponential growth function
f=  		2.77 #pf #degradation rate of the FLC mRNA amount
l=          2.77 #pfl

# VIN3 parameters
vV=4 #ratio of spliced to unspliced VIN3
warm=15 #temperature bound for short term temperature memory
dV=18 #degradation rate of spliced transcript
dv=0
A1=0.75 #S1 #constant for short term memory
BT=17 #temperature bound for long-term temperature memory
dB=0.009 #dLL #parameter for long term temperature memory
CT1=8 #temperature bound for current temperature
CT2=15.4 #temperature bound for current temperature
C1=0.0315 #pC1 #constant for current 
C2=0.03 #pC2 #constant for current 
D1=2.05 #pD #constant for diurnal
sv=vV*dV #splicing rate

# related to light and time of day
lightON=10 #time at dawn
lightOFF=18 #time at dusk
moonSt=16 #resetting time for short term clock
moonEnd=18 
sunrise=np.array([[0,lightON,lightOFF]])
sunrise=sunrise.astype(float)

# initial conditions
tEnd=0 
wi=1;# vi=[0, 0];
COLFRInv=(1-n0i)*r1/(s1+r1) #initial FLC value
n1i=(1-n0i)*s1/(s1+r1)
init=[];#to be defined later

# example parameters and temperature conditions
tiTi=np.array([[-5,-0.01,0,28,28.01,100],[22,22,5,5,22,22]]) #[time][temp]
param=[s1,s2,s3,sel,r1,r2,k,n0i,Tq1,Tq2,nT1,nT2,kFLC, f]

print("%%%%%%%%%%%%%%%%%parameters imported successfully%%%%%%%%%%%%%%")
