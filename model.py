#Written by: Dr Rea Laila Antoniou Kourounioti

import numpy as np
import parameters as p
import math

def model(y,time,tiTi):

# where y = p.init
	# variables
	b=y[0] #is wi
	v=y[1]
	V=y[2]
	Hs=y[3]
	Is=y[4]
	Ns=y[5]
	Ss=y[6]
	Ps=y[7]
	#dividing cells
	H=y[8]
	I=y[9]
	N=y[10]
	S=y[11]
	P=y[12]
	#nondividing cells
	n=Hs+Is+Ns+Ss+Ps+H+I+N+S+P;
	#total cells
	r=y[13]
	#reactivation 
	FLC=y[14]
	#FLC not divided by intial condition



	# light
	if len(p.sunrise)>1: #if array contains more than one value
		if time>p.lightON/24: #if sunrise time is after dawn of the first day
			if p.sunrise[math.floor(time-p.lightON/24),0]==math.floor(time-p.lightON/24):
				day=p.sunrise[math.floor(time-p.lightON/24),] #calculates the day corresponding to the sunrise time
			else:
				print("sunrise position and day don't match")
				print("day "+str(p.sunrise[(p.sunrise[:,0]==math.floor(time-p.lightON/24))])+", time "+str(time))
				day=p.sunrise[(p.sunrise[:,0]==math.floor(time-p.lightON/24))][0]
		else: #if time is before dawn of the first day
			if p.sunrise[math.floor(time),0]==math.floor(time):
				day=p.sunrise[math.floor(time),] #calculates the day corresponding to the sunrise time for times during the first day
			else:
				print("sunrise position and day don't match")
				print("day "+str(p.sunrise[(p.sunrise[:,0]==math.floor(time))])+", time "+str(time))
				day=p.sunrise[(p.sunrise[:,0]==math.floor(time))][0]
	else: #if array only contains one value, then that is the value we use
		day=p.sunrise[0]
		
	if day.size==0: #produces an error if no sunrise values can be found
		print("problem extracting sunrise at time:")
		print(time)
		return -1
	
	p.lightON=day[1] #sets light values to the values needed for the corresponding day (determined by the sunrise time)
	p.lightOFF=day[2]
	p.moonSt=day[2] #as the end of the day and start of the night are the same

	
	# temperature
	T=np.interp(time,tiTi[0],tiTi[1]) #evaluate time as a linear interpolation between tiTi[0] and tiTi[1]
	nT=nightTime12h(time,p.lightON,p.lightOFF,tiTi) #calculates the night time temperature

	# parameters
	a=(1+(-1+p.A1)*(MAXbefore(tiTi,time,p.moonSt)>p.warm))
	#short term temperature memory #if it is cold, a=1, if it is warm,a=p.A1
	c=p.C1-((T-p.CT1)/(p.CT2-p.CT1)*(T<p.CT2)*(T>p.CT1)+(T>=p.CT2))*p.C2;
	#current temperature sensing
	d=(p.D1+math.sin(2*math.pi*(time-(p.lightON-1)/24)))**2;
	#diurnal regulation
	pv=a*b*c*d;
	#VIN3 production
	days1=math.floor(time)

	###### CHANGED S1 HERE!!!! Delete Hash in line below to get new s1 function s1=p.s1+(T>=p.nT2 and days1>=tiTi[0][3])*(1-(Hs+H)/n)*((time**(p.s1))/10)
	s1=p.s1#+(T>=p.nT2 and days1>=tiTi[0][3])*(1-(Hs+H)/n)*((time**(p.s1))/10)
	######
	
	s2=p.s2*(p.Tq2-T)*(T-p.Tq1)*(T>p.Tq1)*(T<p.Tq2)*V #VIN3 dependent transition rate (I to N) #only returns a value if temp is in the correct range, if not, returns 0
	g=r*p.g1*math.exp(p.growthe*(T-22)) 
	dn=p.dn #constant number of non-dividing cells (based on each non-stem cell dividing 5 times)
	if (FLC<p.kFLC)&(T>p.Tq2): #determines value of ratio of initial to current growth rate constant k
		k=p.k
	else:
		k=0
	
	r1=p.r1*((nT-p.nT1)/(p.nT2-p.nT1)*(nT<p.nT2)*(nT>p.nT1)+(nT>=p.nT2)) #determines r1 (VIN3 independent transition rate for H to I)
	s3=p.s3*g #rate of S to P shutdown (VIN3 dependent in dividing cells only)?
	r2=p.r2*g #Reversal rate constant VIN3 dependent (P to I (in dividing cells only))

	## ODEs
	out=[(T<p.BT)-p.dB*b,#dB #dL/dT #ODE for long term temperature memory
		pv-p.sv*v-p.dv*v,#v
		p.sv*v-p.dV*V,#V
		-(s1+s2*p.sel)*Hs+r1*Is,#Hs
		s1*Hs-(r1+s2)*Is+r2*Ps,#Is
		s2*p.sel*Hs+s2*Is-g*Ns,#Ns
		g*Ns-s3*Ss,#Ss
		s3*Ss-r2*Ps,#Ps
		dn*g*Hs-(s1+s2*p.sel)*H+r1*I,#H
		dn*r2*Ps+dn*g*Is+s1*H-(r1+s2)*I,#I
		s2*p.sel*H+s2*I,#N
		dn*g*Ns+dn*(g-s3)*Ss,#S
		dn*s3*Ss+dn*(g-r2)*Ps,#P
		-k*r,#r
		p.f*(Hs+H)/n - p.l*FLC]#FLC rate of change, deg before FLC 
		
	p.tEnd=time
	return out
	


def output(y):
	FLCexp=y[:,14]/p.COLFRInv #divides experimental FLC values by initial FLC value
	VIN3=y[:,2]
	return VIN3,FLCexp
#adjusted to intial condition of FLC in COLFRI


def nightTime12h(time,lightON,lightOFF,tiTi):
	midday=np.mean([lightON,lightOFF])
	#midday is defined as half way between the light going on and it turning back off again
	morningTime=0.25+(midday-12)/24 #morning is defined as being 6hrs before midday
	night=math.floor(time+1-morningTime)-1-0.5+morningTime #night is defined as being 12hrs before morning on the day that we are calculating for if the current time is not during the night (if the current time is during the night we calculate for the current night (day before))
	if night<tiTi[0][0]: #if night occurs before start of the cold period, then night time temp is 22
		nT=22
	else: #if it occurs after the start of the cold period
		morning=math.floor(time+1-morningTime)-1+morningTime #calculate the morning time for the night we are investigating (the morning before if it is currently night, or the next morning if it is not)
		Tnight=np.interp(night,tiTi[0],tiTi[1]) #night time temp is calculated as a linear interpolation between temp values in the tiTi array 
		Tmorning=np.interp(morning,tiTi[0],tiTi[1]) #morning temp is calculated as a linear interpolation between temp values in the tiTi array
		nightIn=np.searchsorted(tiTi[0],night,side='right')#where(tiTi[0]>night)[0][0] #returns the position where the value of night would be inserted into the array
		morningIn=np.searchsorted(tiTi[0],morning)-1#where(tiTi[0]<morning)[0][len(np.where(tiTi[0]<morning)[0])-1] #returns the position where the value of morning would be inserted into the array
		if nightIn<morningIn: #if the position of night is smaller than the morning position (if the calculated night occurs before the calculated morning)
			Tsum=np.array(tiTi[1][nightIn:(morningIn+1)]).cumsum() #creates a running sum of the temperatures occuring between the time of the night and the time of the morning
			Tsum[2:]=Tsum[2:]-Tsum[:-2] 
			nT=((Tnight+tiTi[1][nightIn])*(tiTi[0][nightIn]-night)/2+np.dot(Tsum[1:],np.diff(tiTi[0][nightIn:(morningIn+1)]))/2+(Tmorning+tiTi[1][morningIn])*(morning-tiTi[0][morningIn])/2)*2
		else:
			if morningIn==nightIn:
				nT=((Tnight+tiTi[1][nightIn])*(tiTi[0][nightIn]-night)/2+(Tmorning+tiTi[1][morningIn])*(morning-tiTi[0][morningIn])/2)*2
			else:
				nT=(Tmorning+Tnight)/2
	return nT


	
def MAXbefore(tiTi,t,moonEnd):
	day=math.floor(t)
	if day<tiTi[0][0]:
		maxBefore=22
		#if the day is before the start of the cold treatment, the temp is 22 degrees C
	else:
		if t%1<=(moonEnd/24+10**(-10)): #if it is before dusk
			if day>0: #and after the start of the cold treatment
				thatDay=tiTi[1][(tiTi[0]>(day-1+moonEnd/24))&(tiTi[0]<=t)] # if start of the night we're investigating happens before the earliest value of tiTi but the time we're investigating for is after the start of tiTi thatDay is tiTi[1][0]
				#thatDay =tiTi[1][1] (if night is before tiTi but day is after) or tiTi[1][0]
				if len(thatDay)>0: #if more than one value
					maxBefore=max(np.amax(np.interp([day-1+moonEnd/24,t],tiTi[0],tiTi[1])),np.amax(thatDay)) #interpolate to find the exact time and take the maximum temp found
				else:
					maxBefore=np.amax(np.interp([day-1+moonEnd/24,t],tiTi[0],tiTi[1])) #find the max temp
			else: #if it is on the first day of the cold treatment (day=0)
				thatDay=tiTi[1][(np.array([math.floor(x) for x in tiTi[0]])==(day))&(tiTi[0]<=t)] #creates an array of values for thatDay where the investigated day matches up with the position in the tiTi array and the tiTi position is before the current investigated time
				if len(thatDay)>0: #if there is more than one value, calculate the max of all the max interpolation values
					maxBefore=max(np.amax(np.interp([day,t],tiTi[0],tiTi[1])),np.amax(thatDay))
				else: #otherwise, caluclate the max ineterpolation value
					maxBefore=np.amax(np.interp([day,t],tiTi[0],tiTi[1]))
		else: #if it is after dusk
			thatDay=tiTi[1][(tiTi[0]>(day+moonEnd/24))&(tiTi[0]<=t)] #set thatDay to either tiTi[1][1] (if the start of the night is before the start point of the tiTi array but the array starts before the timepoint) or tiTi[1][0]
			if len(thatDay)>0: 
				maxBefore=max(np.amax(np.interp([day+moonEnd/24,t],tiTi[0],tiTi[1])),np.amax(thatDay)) #if there is more than one value, calculate the max of all the max interpolation values
			else:
				maxBefore=np.amax(np.interp([day+moonEnd/24,t],tiTi[0],tiTi[1])) #otherwise, caluclate the max ineterpolation value
	return maxBefore

