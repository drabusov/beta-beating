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
		self.solution = 0.0
		self.success = False
		
	def getSolution(self):
		print(self.solution)
		return self.solution

	def getLattice(self):
		return self.lattice

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
			self.solution = vec.x
			self.success = True
		except:
			self.setCorrectors(np.zeros(n_corr))
			self.success = False

	def findElements(self,pattern):		

		nodes = self.lattice.getNodes()
		for node in nodes:
			name = node.getName()
			type=node.getType()
			if type == "quad teapot" or type == "multipole teapot":
				if pattern in name.lower():
					position = self.lattice.getNodePositionsDict()[node]
					self.corrDict[name] = position
				else:
					self.quadDict[name] = node.getParam("kq")

			if type == "monitor teapot":
				s_i,s_f = self.lattice.getNodePositionsDict()[node]
				self.bpmList.append((name,round(s_i,6),round(s_f,6)))


		# test for periodicity is required here 
		sArr = self.getPositions()
		s,i=0.0,0
		for bpm in self.bpmList:
			#s = sArr[i]
			while s<bpm[1] or s>bpm[2]:
				i+=1
				s = sArr[i]
			self.bpmIndex.append(i)

		n1,n2 = len(self.bpmList),len(self.bpmIndex)

		if n1 != n2:
			print("Problem with BPMs. There are {} in the lattice, {} found".format(n1,n2))
			#sys.exit(1)


	def setCorrectors(self,theta):

		for i,name in enumerate(self.corrDict.keys()):
			elem =self.lattice.getNodeForName(name)
			type = elem.getType()
			if type=="quad teapot":
				elem.setParam("kq", theta[i])
			else:
				if type=="multipole teapot":
					elem.setParam("kls", [0,theta[i]])
				else:
					print("Problem with element {} occured. It's {}, but has to be quad/multipole".format(name,type))

	def getBetaBpm(self,TwissDataX,TwissDataY):

		betaX = [x[1] for i,x in enumerate(TwissDataX[-1]) if i in self.bpmIndex]
		betaY = [y[1] for i,y in enumerate(TwissDataY[-1]) if i in self.bpmIndex]

		return np.array(betaX),np.array(betaY)


	# calculates the value which characterises the periodicity (symmetry) of the lattice
	# BPMs have to be located periodically
	def periodicBeta(self,theta):

		self.setCorrectors(theta)

		matrix_lattice = TEAPOT_MATRIX_Lattice(self.lattice,self.bunch)
		TwissDataX,TwissDataY = matrix_lattice.getRingTwissDataX(),matrix_lattice.getRingTwissDataY()
		betaX,betaY = self.getBetaBpm(TwissDataX,TwissDataY)

#		metric_x = (np.max(betaX)-np.min(betaX))*np.std(betaX)
#		metric_y = (np.max(betaY)-np.min(betaY))*np.std(betaY)
		metric_x = (1.0-np.min(betaX)/np.max(betaX))
		metric_y = (1.0-np.min(betaY)/np.max(betaY))

		f_min,f_max=sorted([metric_x,metric_y])

		if f_max-f_min > f_min:
			metric = 2*f_max
		else:
			metric = metric_x+metric_y    
		return metric


	def setQuads(self,x):
		kd,kf = x[0],x[-1]
		for name in self.quadDict.keys():
			elem =self.lattice.getNodeForName(name)
			type = elem.getType()
			if type=="quad teapot" and "d" in name:
				elem.setParam("kq", kd)
			if type=="quad teapot" and "f" in name:
				elem.setParam("kq", kf)

	def getTunesDeviation(self,x,qx0,qy0):
		self.setQuads(x)
		matrix_lattice = TEAPOT_MATRIX_Lattice(self.lattice,self.bunch)

		dqx=matrix_lattice.getRingTwissDataX()[0][-1][1]-qx0
		dqy=matrix_lattice.getRingTwissDataY()[0][-1][1]-qy0

		metric = dqx**2+dqy**2

		return metric


	def fixTunes(self,rhobeg=5e-5,A=0.0069,qx0=0.1234,qy0=0.1234):
	
		print("Numbers of main quads found {}".format(len(self.quadDict)))
		
		for key,val in self.quadDict.items():
			if "d" in key:
				kd=val
			if "f" in key:
				kf=val
		x = [kd,kf]
		cons = [{"type": "ineq", "fun": lambda x: A*len(x)-np.sum(x**2)}]
		optionsDict = {'rhobeg':rhobeg, 'disp': True}
		arg=(qx0,qy0)

		try:
			vec = som(self.getTunesDeviation, x, method="COBYLA", constraints=cons, options=optionsDict,args=arg)
			self.setQuads(vec.x)
			self.success = True
		except:
			self.setQuads(kd,kf)
			self.success = False





