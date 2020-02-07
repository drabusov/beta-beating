#!/usr/bin/env python


import numpy as np
import random
import sys
from bunch import Bunch
from orbit.teapot import TEAPOT_MATRIX_Lattice


# ATTENTION !!! The python packet numpy and scipy are required
import numpy as np
from scipy.optimize import minimize as som

class betaCorrection(object):
	""" 
	This routine corrected the closed orbit. Simple method. 
	"""
	
	def __init__(self, lattice, bunch):
		
		self.lattice = lattice
		self.bunch = bunch
		self.corrDict = {}
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


	def hyperSurface(self,theta,qx0,qy0):

		theta=self.normalize(theta[:])
		print("hyper-surface computation")

		nQuadGroups = len(self.quadDict)
		self.setQuads(theta[:nQuadGroups])
		self.setCorrectors(theta[nQuadGroups:])

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
	def normalize(theta):
		c0 = 0.01
		coeff =[[c0,0.,0.,0.],
			[0.,c0,0.,0.],
			[0.,0.,1.,0.],
			[0.,0.,0.,1.]]
		theta_out = np.dot(coeff,theta)
		print("dk is {}".format(theta_out))
		return theta_out
		
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
							self.quadDict[key][name] = d[name]	# the structure is: quadDict = { key:{name1:k_i,name2:k_i,...}, ... }
						else: 
							self.quadDict[key] = d
					else:
						self.corrDict[name] = d[name] 					 

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
							self.quadDict[key][name] = d[name]	# the structure is: quadDict = { key:{name1:k_i,name2:k_i,...}, ... }
						else: 
							self.quadDict[key] = d
					else:
						self.corrDict[name] = d[name] 					 

			if type == "monitor teapot":
				# by definition any monitor in PyORBIT has zero len
				self.bpmDict[name] = np.mean(positionsDict[node])

		if not self.corrDict and not self.quadDict:
			print("Problem with quadrupoles. Incorect user settings")
			print("You can either provide a quadsDict={'group_i':[name0,name1...], ... , 'correctors':list(names)} or patternList = ['corr','kf','kd'...]")

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

		if n1 != n2:
			print("Problem with BPMs. There are {} in the lattice, {} found".format(n1,n2))
			#sys.exit(1)


	def setQuads(self,x):

		i=0
		for key,d in self.quadDict.items():

			tmp =zip(d.keys(),np.array(d.values())+x[i]*np.ones(len(d.values())))
			tmpDict = dict(tmp)

			for name,val in tmpDict.items():
				elem =self.lattice.getNodeForName(name)
				type = elem.getType()

				if type=="quad teapot":
					elem.setParam("kq", val)
				if type=="multipole teapot":
					elem.setParam("kls", [0,val])
				if type not in ["quad teapot","multipole teapot"]:
					print("Problem with element {} occured. It's {}, but has to be quad/multipole".format(name,type))
			i+=1


	def setCorrectors(self,theta):

		tmp =zip(self.corrDict.keys(),np.array(self.corrDict.values())+theta)
		tmpDict = dict(tmp)
		for name,val in tmpDict.items():
			elem =self.lattice.getNodeForName(name)
			type = elem.getType()
			if type=="quad teapot":
				elem.setParam("kq", val)
			else:
				if type=="multipole teapot":
					elem.setParam("kls", [0,val])
				else:
					print("Problem with element {} occured. It's {}, but has to be quad/multipole".format(name,type))




	def getBetaBpm(self,TwissDataX,TwissDataY):

		betaX = [x[1] for i,x in enumerate(TwissDataX[-1]) if i in self.bpmIndex]
		betaY = [y[1] for i,y in enumerate(TwissDataY[-1]) if i in self.bpmIndex]

		return np.array(betaX),np.array(betaY)

	def getBeta(self,TwissDataX,TwissDataY):
		return np.transpose(TwissDataX[-1])[-1],np.transpose(TwissDataY[-1])[-1]


	def correction(self,qx0=0.1234,qy0=0.1234,patternDict=dict(),quadsDict = dict(),betaX0=1.0,betaY0=1.0):

		method = "COBYLA"

		print(self.quadDict)

		self.findElements(patternDict,quadsDict)
		nCorr=len(self.corrDict)
		nQuadGroups=len(self.quadDict)

		L = self.lattice.getLength()
		sVar = np.linspace(0,L,len(betaX0))

		print("Numbers of monitors found {}".format(len(self.bpmDict)))
		print("Numbers of correctors found {}".format(nCorr))
		print("Numbers of main families found {}".format(nQuadGroups))

		print(self.quadDict)
		
		theta = np.zeros(nCorr+nQuadGroups)

		# reasonable initial changes is on the level of 0.1% of quad strength (jac computations or the first step)  
		epsilon = 1e-3*np.mean(np.abs([np.mean(group.values()) for group in self.quadDict.values()]))
		cons = [{"type": "eq", "fun":self.hyperSurface, "args":(qx0,qy0)}]		
		optionsDict = {'disp': True, "eps": epsilon, "ftol": 1e-5} 
		arg=(sVar,betaX0,betaY0)

		if method == "COBYLA":
			for con in cons:
				con["type"]="ineq"
			optionsDict = {'rhobeg':epsilon, 'catol':1e-6, "tol":1e-3,'disp': True}

		fIni = self.observable(theta,*arg)
		fftX,fftY = self.calculateBeat(self.lattice,self.bunch,sVar,betaX0,betaY0)
		aIni = self.calculateAperiodic(fftX,fftY)

		vec = som(self.observable, theta, method=method, constraints=cons, options=optionsDict,args=arg)

		self.solution = vec.x
		# careful! you set dk here! not an absolute value!
		fCorr = self.observable(vec.x,*arg)

		if fCorr>fIni:
			# recursion? + random coeff*epsilon
			optionsDict["rhobeg"]=5*epsilon
			vec = som(self.observable, theta, method=method, constraints=cons, options=optionsDict,args=arg)
			
		self.solution = vec.x
		# careful! you set dk here! not an absolute value!
		fCorr = self.observable(vec.x,*arg)

		fftX,fftY = self.calculateBeat(self.lattice,self.bunch,sVar,betaX0,betaY0)
		aCorr = self.calculateAperiodic(fftX,fftY)

		return fIni,fCorr,aIni,aCorr



	def observable(self,theta,sVar,betaX0,betaY0):

		print("observable/metric computation")

		theta=self.normalize(theta[:])
		
		nQuadGroups = len(self.quadDict)
		self.setQuads(theta[:nQuadGroups])
		self.setCorrectors(theta[nQuadGroups:])

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

		nQuadGroups=len(self.quadDict)
		
		theta = np.zeros(nQuadGroups)
		optionsDict = {'rhobeg':rhobeg,'maxiter':2000, 'disp': True}
		arg=(qx0,qy0)

		if warm:
			vec = som(self.getWarmTunes, theta, method="COBYLA", options=optionsDict,args=arg)
		else:
			vec = som(self.getTunesDeviation, theta, method="COBYLA",options=optionsDict,args=arg)

		self.setQuads(vec.x)
		self.solution = vec.x
		print


	def getTunesDeviation(self,x,qx0,qy0):

		kw1=self.lattice.getNodeForName("s52qd11").getParam("kq")
		kw2=self.lattice.getNodeForName("s52qd12").getParam("kq")
		print("warm quads are {} {}".format(kw1,kw2))

		self.setQuads(x)

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
	def getWarmTunes(self,x,qx0,qy0):

		self.setQuads(x)

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

		self.setCorrectors(theta[:-2])
		self.setQuads(theta[-2:])

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



