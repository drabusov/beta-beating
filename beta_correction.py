#!/usr/bin/env python


import numpy as np
import random
import sys
from bunch import Bunch
from orbit.teapot import TEAPOT_MATRIX_Lattice


# ATTENTION !!! The python libraries numpy and scipy are required
import numpy as np
from scipy.optimize import minimize as som

class betaCorrection(object):
	""" 
	This routine corrected the closed orbit. Simple method. 
	"""
	
	def __init__(self, lattice, bunch):
		
		self.lattice = lattice
		self.bunch = bunch
		self.quadDict = {}

		self.bpmDict = {}
		self.bpmIndex = [] # positions of bpms in TwissData-like output

		self._solution = None
		
	@property
	def solution(self):
		return self._solution

	@solution.setter
	def solution(self, arr):
		self._solution = arr

	# returns the s-array from TwissData-like output
	def getPositions(self):
		matrix_lattice = TEAPOT_MATRIX_Lattice(self.lattice,self.bunch)
		TwissDataX = matrix_lattice.getRingTwissDataX()
		sArr = [x[0] for x in TwissDataX[0]]
		return sArr


	def hyperSurface(self,theta,qx0,qy0,alpha):

		dk=self.normalize(theta,alpha,self.quadDict)
		print("hyper-surface computation")

		self.setQuads(dk,self.quadDict)

		matrix_lattice = TEAPOT_MATRIX_Lattice(self.lattice,self.bunch)
		(mux,bTmp,aTmp) = matrix_lattice.getRingTwissDataX()
		(muy,bTmp,aTmp) = matrix_lattice.getRingTwissDataY()
		qx,qy = mux[-1][1],muy[-1][1]
		print("current tunes are {} {}".format(qx,qy))
		out = np.sqrt((qx-qx0)**2+(qy-qy0)**2)

		print("the distance from set tunes hor:{} ver:{}".format(qx-qx0,qy-qy0))
		print("out is {}\n".format(out))

		return -out 


	# step size, weak focusing approximation
	@staticmethod
	def normalize(theta, alpha, groupDict):

		print("theta is {}".format(theta))

		tmp = zip(theta,groupDict.items())
		dk =[x*alpha/float(len(d)) for x,(key,d) in tmp]
		
		print("dk is {}".format(dk))
		return np.array(dk)
		
	def findElements(self, patternDict, quadsDict):

		quadsDictUpd = {}
		for key,namesList in quadsDict.items():
			quadsDictUpd[key.lower()] = [name.lower() for name in namesList]

		patternDictUpd = {}
		for pattern,key in patternDict.items():
			patternDictUpd[pattern.lower()] = key.lower()

		positionsDict = self.lattice.getNodePositionsDict()
		nodes = self.lattice.getNodes()
		for node in nodes:
			name = node.getName()
			type=node.getType()

			for key,namesList in quadsDictUpd.items():
				d = {} # tmp dict
				if name in namesList:
					if type == "quad teapot":
						d[name] = node.getParam("kq")
					if type == "multipole teapot":
						d[name] = node.getParam("kls")[1]
					if type not in ["quad teapot", "multipole teapot"]:
						print("problem with the element {}. It's not a quad/multipole".format(name))
				if d:
					if "corr" not in key:
						if self.quadDict.has_key(key):
							self.quadDict[key][name] = d[name]	# the structure is: quadDict = { key:{name1:k_i,name2:k_i,...}, ..., corr1:{corr1:kc} }
						else: 
							self.quadDict[key] = d
					else:
						self.quadDict["corr_{}".format(name)] = d 					 

				 

			for pattern,key in patternDictUpd.items():
				d = {} # tmp dict
				if pattern in name.lower():
					if type == "quad teapot":
						d[name] = node.getParam("kq")
					if type == "multipole teapot":
						d[name] = node.getParam("kls")[1]
					if type not in ["quad teapot", "multipole teapot"]:
						print("problem with the element {}. It's not a quad/multipole".format(name))
				if d:
					if "corr" not in key:
						if self.quadDict.has_key(key):
							self.quadDict[key][name] = d[name]	# the structure is: quadDict = { key:{name1:k_i,name2:k_i,...}, ..., corr1:{corr1:kc} }
						else: 
							self.quadDict[key] = d
					else:
						self.quadDict["corr_{}".format(name)] = d 					 

			if type == "monitor teapot":
				self.bpmDict[name] = np.mean(positionsDict[node])

		if not self.quadDict:
			print("Problem with quadrupoles. Incorect user settings")
			print("You can either provide a quadsDict={'group_i':[name0,name1...], ... , 'correctors':list(names)} or patternDict = {group_name:'pattern'}")

		# test for periodicity is required here 
		sArr = self.getPositions()
		I=0
		for name,position in self.bpmDict.items():
			for idx,s in enumerate(sArr[I:]):
				if np.abs(round(s,6)-round(position,6)) < 10**(-4):
					self.bpmIndex.append(idx)
					I = idx
					#print(I)
					break
		n1,n2 = len(self.bpmDict),len(self.bpmIndex)

#		if n1 != n2:
#			print("Problem with BPMs. There are {} in the lattice, {} found".format(n1,n2))
#			#sys.exit(1)


	def setQuads(self,x,quadDict):

		for xi,(group_name,d) in zip(x,quadDict.items()):

			for name,val in d.items():
				elem =self.lattice.getNodeForName(name)
				type = elem.getType()

				if type=="quad teapot":
					elem.setParam("kq", val+xi)
				if type=="multipole teapot":
					elem.setParam("kls", [0,val+xi])
				if type not in ["quad teapot","multipole teapot"]:
					print("Problem with element {} occured. It's {}, but has to be quad/multipole".format(name,type))


	def getBetaBpm(self,TwissDataX,TwissDataY):

		betaX = [x[1] for i,x in enumerate(TwissDataX[-1]) if i in self.bpmIndex]
		betaY = [y[1] for i,y in enumerate(TwissDataY[-1]) if i in self.bpmIndex]

		return np.array(betaX),np.array(betaY)

	def getBeta(self,TwissDataX,TwissDataY):
		return np.transpose(TwissDataX[-1])[-1],np.transpose(TwissDataY[-1])[-1]


	def correction(self,qx0=0.1234,qy0=0.1234,patternDict=dict(),quadsDict = dict(),betaX0=1.0,betaY0=1.0):

		#method = "COBYLA"
		method = "SLSQP"

		L = self.lattice.getLength()
		sVar = np.linspace(0,L,len(betaX0))

		if not self.quadDict:
			self.findElements(patternDict, quadsDict)
		print(self.quadDict)
		
		theta = np.zeros(len(self.quadDict))

		epsilon,ftol = 1e-1,1e-3
		# a normalization conts
		alpha = (qx0**2+qy0**2)*(2*np.pi/L)**2
		print(alpha)
		cons = [{"type": "eq", "fun":self.hyperSurface, "args":(qx0,qy0,alpha)}]		
		optionsDict = {'disp': True, "eps": epsilon, "ftol": ftol} 
		arg=(sVar,betaX0,betaY0,alpha)

		if method == "COBYLA":
			for con in cons:
				con["type"]="ineq"
			optionsDict = {'rhobeg':epsilon, 'catol':1e-6, "tol":ftol,'disp': True}

		fIni = self.observable(theta,*arg)
		fftX,fftY = self.calculateBeat(self.lattice,self.bunch,sVar,betaX0,betaY0)
		aIni = self.calculateAperiodic(fftX,fftY)

		vec = som(self.observable, theta, method=method, constraints=cons, options=optionsDict,args=arg)

		self.solution = vec.x
		# careful! you set dk here! not an absolute value!
		fCorr = self.observable(vec.x,*arg)
		fftX,fftY = self.calculateBeat(self.lattice,self.bunch,sVar,betaX0,betaY0)
		aCorr = self.calculateAperiodic(fftX,fftY)
	
		self.solution = vec.x

		return fIni,fCorr,aIni,aCorr



	def observable(self,theta,sVar,betaX0,betaY0,alpha):

		print("observable/metric computation")

		dk=self.normalize(theta,alpha,self.quadDict)
		
		self.setQuads(dk,self.quadDict)

		fftBeatX,fftBeatY=self.calculateBeat(self.lattice,self.bunch,sVar,betaX0,betaY0)

		metric=self.calculateMax(fftBeatX,fftBeatY)
#		metric=self.calculateMax(fftBeatX,fftBeatY)
		print("the observable/metric is {}\n".format(metric))

		return metric


	@staticmethod
	def calculateBeat(lattice,bunch,sVar,betaX0,betaY0):

		matrix_lattice = TEAPOT_MATRIX_Lattice(lattice,bunch)
		TwissDataX,TwissDataY = matrix_lattice.getRingTwissDataX(),matrix_lattice.getRingTwissDataY()
		betaX,betaY = np.transpose(TwissDataX[-1]),np.transpose(TwissDataY[-1])
		qx,qy = np.transpose(TwissDataX[0])[-1][-1],np.transpose(TwissDataY[0])[-1][-1]
		print("current tunes are {} {}".format(qx,qy))

		betX = np.interp(x=sVar, xp=betaX[0],fp=betaX[1])
		betY = np.interp(x=sVar, xp=betaY[0],fp=betaY[1])

		beatX = (betX-betaX0)
		beatY = (betY-betaY0)
		n = len(beatX)

		fftBeatX =np.abs(np.fft.rfft(beatX))/n
		fftBeatY =np.abs(np.fft.rfft(beatY))/n

		return fftBeatX,fftBeatY


	@staticmethod
	def calculateMax(fftX,fftY):
		return np.sqrt(np.max(fftX)**2 + np.max(fftY)**2)

	@staticmethod
	def calculateAperiodic(fftX,fftY,S=6):

		AperiodicHarmX = [x for i,x in enumerate(fftX) if i%S]
		AperiodicHarmY = [x for i,x in enumerate(fftY) if i%S]

		return np.sqrt(np.sum(AperiodicHarmX)**2 + np.sum(AperiodicHarmY)**2)


	def matchTunes(self,rhobeg=5e-5,A=0.001,qx0=0.1234,qy0=0.1234,quadsDict = dict(),patternDict= dict(),warm=False):

		if not self.quadDict:
			self.findElements(patternDict,quadsDict)

		groupsDict = {key:val for key,val in self.quadDict.items() if key.lstrip("corr_")==key}
		print(groupsDict)

		if len(groupsDict)>2:
			print("warning! number of variables > number of equations")
		
		theta = np.zeros(len(groupsDict))
		optionsDict = {'rhobeg':rhobeg,'maxiter':2000, 'disp': True}
		arg=(qx0,qy0,groupsDict)

		if warm:
			vec = som(self.getWarmTunes, theta, method="COBYLA", options=optionsDict,args=arg)
		else:
			vec = som(self.getTunesDeviation, theta, method="COBYLA",options=optionsDict,args=arg)

		self.setQuads(vec.x, groupsDict)
		self.solution = vec.x
		print


	def getTunesDeviation(self,x,qx0,qy0,groupsDict):

		kw1=self.lattice.getNodeForName("s52qd11").getParam("kq")
		kw2=self.lattice.getNodeForName("s52qd12").getParam("kq")
		print("warm quads are {} {}".format(kw1,kw2))

		self.setQuads(x,groupsDict)

		kw1=self.lattice.getNodeForName("s52qd11").getParam("kq")
		kw2=self.lattice.getNodeForName("s52qd12").getParam("kq")
		print("warm quads are {} {}".format(kw1,kw2))

		matrix_lattice = TEAPOT_MATRIX_Lattice(self.lattice,self.bunch)
		
		(mux,b,a) = matrix_lattice.getRingTwissDataX()
		(muy,b,a) = matrix_lattice.getRingTwissDataY()
		qx1,qy1 = mux[-1][1],muy[-1][1]
		#print("{} {}".format(qx1,qy1))
		#print("{} {}\n".format(qx0,qy0))
		dqx=1-qx1/qx0
		dqy=1-qy1/qy0

		metric = dqx**2+dqy**2

		return metric

	# for sis100 only
	def getWarmTunes(self,x,qx0,qy0,groupsDict):

		self.setQuads(x,groupsDict)

		kqd = self.lattice.getNodeForName("s11qd11").getParam("kq")
		kqf = self.lattice.getNodeForName("s11qd12").getParam("kq")

		kWarm1=1.0139780*kqd*1.3/1.76
		kWarm2=1.0384325*kqf*1.3/1.76

		kqd = self.lattice.getNodeForName("s52qd11").setParam("kq", kWarm1)
		kqf = self.lattice.getNodeForName("s52qd12").setParam("kq", kWarm2)

		matrix_lattice = TEAPOT_MATRIX_Lattice(self.lattice,self.bunch)
		
		(mux,b,a) = matrix_lattice.getRingTwissDataX()
		(muy,b,a) = matrix_lattice.getRingTwissDataY()
		qx1,qy1 = mux[-1][1],muy[-1][1]
		print("{} {}".format(qx1,qy1))
		print("{} {}\n".format(qx0,qy0))
		dqx=1-qx1/qx0
		dqy=1-qy1/qy0

		metric = dqx**2+dqy**2

		return metric

	def periodicFun(self,theta,qx0,qy0):

		self.setQuads(theta)

		matrix_lattice = TEAPOT_MATRIX_Lattice(self.lattice,self.bunch)
		TwissDataX,TwissDataY = matrix_lattice.getRingTwissDataX(),matrix_lattice.getRingTwissDataY()
		betaX,betaY = self.getBetaBpm(TwissDataX,TwissDataY)

		metric_x = 1.0-np.min(betaX)/np.max(betaX)
		metric_y = 1.0-np.min(betaY)/np.max(betaY)

		arr = [metric_x,metric_y]

		if f_max-f_min > f_max:
			metric = 4*f_max
		else:
			metric = np.sum(arr)    
		return metric

	def hyperBetaBeat(self,rhobeg=5e-5,A=0.1,qx0=0.1234,qy0=0.1234,k=tuple(),patternDict=dict(),quadsDict = dict(),betaX0=1.0,betaY0=1.0):

		#if not self.quadDict:
		self.findElements(patternDict,quadsDict)

		nQuadGroups=len(self.quadDict)

		self.lattice.getNodeForName("s52qd11").setParam("kq",k[0])
		self.lattice.getNodeForName("s52qd12").setParam("kq",k[1])

		kw1=self.lattice.getNodeForName("s52qd11").getParam("kq")
		kw2=self.lattice.getNodeForName("s52qd12").getParam("kq")
		print("warm quads are {} {}".format(kw1,kw2))
		k1=self.lattice.getNodeForName("s11qd11").getParam("kq")
		k2=self.lattice.getNodeForName("s11qd12").getParam("kq")
		print("cold quads are {} {}".format(k1,k2))
				
		theta = np.zeros(nQuadGroups)
		optionsDict = {'rhobeg':rhobeg,'maxiter':2000, 'disp': True}
		arg=(qx0,qy0)

		vec = som(self.getTunesDeviation, theta, method="COBYLA", options=optionsDict,args=arg)

		self.setQuads(vec.x)
		kd=self.lattice.getNodeForName("s11qd11").getParam("kq")
		kf=self.lattice.getNodeForName("s11qd12").getParam("kq")
		print("sollution found\n")

		matrix_lattice = TEAPOT_MATRIX_Lattice(self.lattice,self.bunch)
		TwissDataX,TwissDataY = matrix_lattice.getRingTwissDataX(),matrix_lattice.getRingTwissDataY()
		betaX,betaY = np.transpose(TwissDataX[-1]),np.transpose(TwissDataY[-1])

		L = self.lattice.getLength()
		sVar = np.linspace(0,L,len(betaX0))

		betxSmpl = np.interp(x=sVar, xp=betaX[0],fp=betaX[1])
		betySmpl = np.interp(x=sVar, xp=betaY[0],fp=betaY[1])

		beat_x = (betxSmpl-betaX0)
		beat_y = (betySmpl-betaY0)
		n = len(beat_x)

		beta_x_fft =np.abs(np.fft.rfft(beat_x))/n
		beta_y_fft =np.abs(np.fft.rfft(beat_y))/n
	
		periodicHarmX = [x for i,x in enumerate(beta_x_fft) if not i%6]
		periodicHarmY = [x for i,x in enumerate(beta_y_fft) if not i%6]	
		AperiodicHarmX = [x for i,x in enumerate(beta_x_fft) if i%6]
		AperiodicHarmY = [x for i,x in enumerate(beta_y_fft) if i%6]

		maxPosX = np.argmax(beta_x_fft)
		maxPosY = np.argmax(beta_y_fft)
 
		return kd,kf,maxPosX,maxPosY,beta_x_fft[maxPosX],beta_y_fft[maxPosY],np.sum(periodicHarmX)/np.sum(AperiodicHarmX),np.sum(periodicHarmY)/np.sum(AperiodicHarmY) 



