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
		self.bpmList = []
		self.bpmIndex = [] # positions of bpms in TwissData-like output
		self.solution = 0.0
		
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
		

		theta = np.zeros(12)
		cons = [{"type": "ineq", "fun": lambda x: A*len(x)-np.sum(x**2)}]
		optionsDict = {'rhobeg':rhobeg, 'disp': True}

		vec = som(self.periodicBeta, theta, method="COBYLA", constraints=cons, options=optionsDict)

		self.setCorrectors(vec.x)
		self.solution = vec.x


	def findElements(self,pattern):		

		nodes = self.lattice.getNodes()
		for node in nodes:
			name = node.getName()
			type=node.getType()
			if pattern in name.lower():
				if type == "quad teapot" or type == "multipole teapot":
					position = self.lattice.getNodePositionsDict()[node]
					self.corrDict[name] = position
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

		metric_x = (np.max(betaX)-np.min(betaX))*np.std(betaX)
		metric_y = (np.max(betaY)-np.min(betaY))*np.std(betaY)
	    
		if np.abs(metric_x - metric_y) > np.max([metric_x,metric_y]):
			if metric_x < metric_y:
				metric = 2*metric_y
			else:
				metric = 2*metric_x
		else:
			metric = metric_x+metric_y    
		return metric









