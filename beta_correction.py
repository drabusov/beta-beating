#!/usr/bin/env python


import numpy as np
import random
import sys
from bunch import Bunch


# ATTENTION !!! The python packet numpy and scipy are required
import numpy as np
from scipy.optimize import minimize as som

from orbit.teapot import TEAPOT_MATRIX_Lattice


class betaCorrection:
	""" 
	This routine corrected the closed orbit. Simple method. 
	"""
	
	def __init__(self, lattice, bunch):
		
		self.lattice = lattice
		self.bunch = bunch
		self.corrDict = {}
		self.quadDict = {}
		self.bpmList = []
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


	# A is the constraint for an individual corrector
	# -A < kick_i < A 
	def constraintFunc(self,x,A=0.05):
	    return A*len(x)-np.sum(x**2)

		
	def correction(self,rhobeg=1e-5,A=0.05,pattern="corr"):
	
		self.findElements(pattern)

		print("Numbers of monitors found {}".format(len(self.bpmList)))
		print("Numbers of correctors found {}".format(len(self.corrDict)))
		
		n_corr = len(self.corrDict.keys())
		theta = np.zeros(n_corr)
		cons = [{"type": "ineq", "fun": lambda x: A*len(x)-np.sum(x**2)}]
		optionsDict = {'rhobeg':rhobeg, 'disp': True}

		try:
			vec = som(self.periodicBeta, theta, method="COBYLA", constraints=cons, options=optionsDict)
			self.setCorrectors(vec.x)
			self._solution = vec.x
			self.success = True
		except:
			self.setCorrectors(np.zeros(n_corr))
			self.success = False

	def findElements(self, patternList = ["corr"], quadsDict={}):		

		positionsDict = self.lattice.getNodePositionsDict()
		nodes = self.lattice.getNodes()
		for node in nodes:
			name = node.getName()
			type=node.getType()

			if quadsDict:
				for key,namesList in quadsDict.items():
					d = {} # tmp dict
					if name in namesList:
						if type == "quad teapot":
							d[name] = node.getParam("kq")
						if type == "multipole teapot":
							d[name] = node.getParam("kls")[1]
						if type not in ["quad teapot", "multipole teapot"]:
							print("problem with the element {}. It's not a quad/multipole".format(name))
						
					if "corr" not in key:
						self.quadDict[key] = d	# the structure is: quadDict = { key:{name:k_i}, ... }
					else:
						self.corrDict = d 					 
			else:
				for pattern in patternList:
					d = {}
					if pattern.lower() in name.lower():
						if type == "quad teapot":
							d[name] = node.getParam("kq")
						if type == "multipole teapot":
							d[name] = node.getParam("kls")[1]
						if type not in ["quad teapot", "multipole teapot"]:
							print("problem with the element {}. It's not a quad/multipole".format(name))
					if "corr" not in key:
						self.quadDict[key] = d	# the structure is: quadDict = { key:{name:k_i}, ... }
					else:
						self.corrDict = d 					 

			if type == "monitor teapot":
				s_i,s_f = positionsDict[node]
				self.bpmList.append((name,round(s_i,6),round(s_f,6)))

		if not self.corrDict and not self.quadDict:
			print("Problem with quadrupoles. Incorect user settings")
			print("You can either provide a quadsDict={'group_i':[name0,name1...], ... , 'correctors':list(names)} or patternList = ['corr','kf','kd'...]")
		print(self.corrDict)
		print(self.quadDict)


		# test for periodicity is required here 
		sArr = self.getPositions()
		s,i=0.0,0
		for bpm in self.bpmList:
			#s = sArr[i]
			while round(s,5)<round(bpm[1],5) or round(s,5)>round(bpm[2],5):
				i+=1
				s = sArr[i]
			self.bpmIndex.append(i)

		n1,n2 = len(self.bpmList),len(self.bpmIndex)

		if n1 != n2:
			print("Problem with BPMs. There are {} in the lattice, {} found".format(n1,n2))
			#sys.exit(1)


	def setQuads(self,x):

		i=0
		for key,d in self.quadDict.items():
			print(key)

			tmp =zip(d.keys(),np.array(d.values())+x[i]*np.ones(len(d.values())))
			tmpDict = dict(tmp)

			for name,val in d.items():
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
		for name,val in tmpDict:
			elem =self.lattice.getNodeForName(name)
			type = elem.getType()
			if type=="quad teapot":
				elem.setParam("kq", val)
			else:
				if type=="multipole teapot":
					elem.setParam("kls", [0,val])
				else:
					print("Problem with element {} occured. It's {}, but has to be quad/multipole".format(name,type))


#		for i,name in enumerate(self.corrDict.keys()):
#			elem =self.lattice.getNodeForName(name)
#			type = elem.getType()
#			if type=="quad teapot":
#				elem.setParam("kq", theta[i])
#			else:
#				if type=="multipole teapot":
#					elem.setParam("kls", [0,theta[i]])
#				else:
#					print("Problem with element {} occured. It's {}, but has to be quad/multipole".format(name,type))




	def getBetaBpm(self,TwissDataX,TwissDataY):

		betaX = [x[1] for i,x in enumerate(TwissDataX[-1]) if i in self.bpmIndex]
		betaY = [y[1] for i,y in enumerate(TwissDataY[-1]) if i in self.bpmIndex]

		return np.array(betaX),np.array(betaY)

	def getBeta(self,TwissDataX,TwissDataY):
		return np.transpose(TwissDataX[-1])[-1],np.transpose(TwissDataY[-1])[-1]

	# calculates the value which characterises the periodicity (symmetry) of the lattice
	# BPMs have to be located periodically
	def periodicBeta(self,theta):

		self.setCorrectors(theta)
		print(theta)

		matrix_lattice = TEAPOT_MATRIX_Lattice(self.lattice,self.bunch)
		TwissDataX,TwissDataY = matrix_lattice.getRingTwissDataX(),matrix_lattice.getRingTwissDataY()
		betaX,betaY = self.getBetaBpm(TwissDataX,TwissDataY)

		metric_x = (1.0-np.min(betaX)/np.max(betaX))
		metric_y = (1.0-np.min(betaY)/np.max(betaY))

		f_min,f_max=sorted([metric_x,metric_y])

		if f_max-f_min > f_min:
			metric = 2*f_max
		else:
			metric = metric_x+metric_y    
		return metric


	def getTunesDeviation(self,x,qx0,qy0):
		self.setQuads(x)
		matrix_lattice = TEAPOT_MATRIX_Lattice(self.lattice,self.bunch)
		
		(mux,b,a) = matrix_lattice.getRingTwissDataX()
		(muy,b,a) = matrix_lattice.getRingTwissDataY()
		qx1,qy1 = mux[-1][1],muy[-1][1]
		print("{} {}".format(qx1,qy1))
		print("{} {}".format(qx0,qy0))
		dqx=1-qx1/qx0
		dqy=1-qy1/qy0

		metric = dqx**2+dqy**2

		return metric


	def fixTunes(self,rhobeg=5e-5,A=0.0069,qx0=0.1234,qy0=0.1234,patternQuad= "quad",patternCorr = "corr",betaX0=1.,betaY0=1.):

		self.findElements(patternQuad,patternCorr)
		nCorr=len(self.corrDict)
		nQuad=len(self.quadDict)

		L = self.lattice.getLength()
		sVar = np.linspace(0,L,len(betaX0))

		print("Numbers of monitors found {}".format(len(self.bpmList)))
		print("Numbers of correctors found {}".format(nCorr))
		print("Numbers of main quads found {}".format(nQuad))

		
		theta = np.zeros(nCorr+nQuad)
#		theta = np.random.normal(0,10**(-5),nCorr+nQuad)
		print(theta)
		cons = [{"type": "ineq", "fun": lambda x: A*len(x)-np.sum(x**2)}]
		optionsDict = {'rhobeg':rhobeg,'maxiter':2000, 'disp': True}
		arg=(qx0,qy0,sVar,betaX0,betaY0,True)

		vec = som(self.FFTbeta, theta, method="COBYLA", constraints=cons, options=optionsDict,args=arg)

		self.setQuads(vec.x[:nQuadGroups])
		self.setCorrectors(vec.x[nQuadGroups:])
		self.solution = vec.x


	def periodicFun(self,theta,qx0,qy0):

		self.setCorrectors(theta[:-2])
		self.setQuads(theta[-2:])

		matrix_lattice = TEAPOT_MATRIX_Lattice(self.lattice,self.bunch)
		TwissDataX,TwissDataY = matrix_lattice.getRingTwissDataX(),matrix_lattice.getRingTwissDataY()
		betaX,betaY = self.getBetaBpm(TwissDataX,TwissDataY)

		(mux,b,a) = TwissDataX
		(muy,b,a) = TwissDataY
		qx1,qy1 = mux[-1][1],muy[-1][1]
		print("{} {}".format(qx1,qy1))
		print("{} {}".format(qx0,qy0))
		f_x=np.abs(1-qx1/qx0)
		f_y=np.abs(1-qy1/qy0)

		metric_x = 1.0-np.min(betaX)/np.max(betaX)
		metric_y = 1.0-np.min(betaY)/np.max(betaY)

		arr = [metric_x,metric_y,f_x,f_y]

		f_min,f_max=np.min(arr),np.max(arr)

		if f_max-f_min > f_max:
			metric = 4*f_max
		else:
			metric = np.sum(arr)    
		return metric


	def FFTbeta(self,theta,qx0,qy0,sVar,betaX0,betaY0,suppressAll):
		
		nQuadGroups = len(self.quadDict)

		self.setQuads(theta[:nQuadGroups])
		self.setCorrectors(theta[nQuadGroups:])

		matrix_lattice = TEAPOT_MATRIX_Lattice(self.lattice,self.bunch)
		TwissDataX,TwissDataY = matrix_lattice.getRingTwissDataX(),matrix_lattice.getRingTwissDataY()

		betaX,betaY = np.transpose(TwissDataX[-1]),np.transpose(TwissDataY[-1])

		betxSmpl = np.interp(x=sVar, xp=betaX[0],fp=betaX[1])
		betySmpl = np.interp(x=sVar, xp=betaY[0],fp=betaY[1])


		(mux,a,b_x) = TwissDataX
		(muy,a,b_y) = TwissDataY
		qx1,qy1 = mux[-1][1],muy[-1][1]
		#print("{} {}".format(qx1,qy1))
		#print("{} {}".format(qx0,qy0))
		f_x=np.abs(1-qx1/qx0)
		f_y=np.abs(1-qy1/qy0)


		beat_x = (betxSmpl/betaX0-1)*100
		beat_y = (betySmpl/betaY0-1)*100
		n = len(beat_x)

		beta_x_fft =np.abs(np.fft.rfft(beat_x))/n
		beta_y_fft =np.abs(np.fft.rfft(beat_y))/n


		periodicHarmX = [x for i,x in enumerate(beta_x_fft) if not i%6]
		periodicHarmY = [x for i,x in enumerate(beta_y_fft) if not i%6]
	
		AperiodicHarmX = [x for i,x in enumerate(beta_x_fft) if i%6]
		AperiodicHarmY = [x for i,x in enumerate(beta_y_fft) if i%6]

		if suppressAll:
			metric_x = np.sum(AperiodicHarmX)
			metric_y = np.sum(AperiodicHarmY) 
		else:
			metric_x = np.sum(AperiodicHarmX)-np.sum(periodicHarmX)
			metric_y = np.sum(AperiodicHarmY)-np.sum(periodicHarmY) 

		arr = [metric_x,metric_y] # only periodicity 
		metric = np.sum(arr)    
#		arr = [metric_x,metric_y,f_x,f_y]
		
		f_min,f_max=np.min(arr),np.max(arr) # you can rewrite it using sorted~!

#		if f_max-f_min > f_min:
#			metric = 2*f_max
#		else:
#			metric = np.sum(arr)    

		print(metric)

		return metric


