import os
import string
import sys
from numpy import *
from scipy.optimize import fsolve
from scipy.optimize import root
from scipy.integrate import odeint
from scipy.constants import c
from matplotlib.pyplot import * # why do you need pyplot here?

#from orbit.orbit_correction import orbit

from copy import copy, deepcopy

from orbit.teapot import TEAPOT_MATRIX_Lattice



class Twiss:
	# Create a simple MAD-like twiss object:  
	def __init__(self):
		self.data = { 'keyword': '',
						's': 0.0,
						'L': 0.0,
						'alfx': 0.0,
						'alfy': 0.0,
						'betx': 0.0,
						'bety': 0.0,
						'mux' : 0.0,
						'muy' : 0.0,
						'Dx': 0.0,
						'Dpx': 0.0,
						'angle': 0.0, 
						'k1': 0.0 }


class Optics:
# An container class for twiss objects:
	def __init__(self):
		self.line = []

	def __len__(self):
		return len(self.line)

	def __getitem__(self,j):
		return self.line[j]

	def __setitem__(self,j,x):
		self.line[j]=x

	def add(self, x):
		self.line.append(x)

	def print_line(self):
		for j in xrange(0,len(self.line)):
			print j, self.line[j].data['keyword'], "s:", self.line[j].data['s'], "L:", self.line[j].data['L'], 360.0*self.line[j].data['mux'],self.line[j].data['bety'],self.line[j].data['alfy'], self.line[j].data['k1'] 

	def get_element(self, s):
		Nb=len(self.line)
		if self.line[0].data['s'] >= s and s >= 0.0: 
			return 0
		for j in xrange(1,Nb):
			if self.line[j-1].data['s'] < s and self.line[j].data['s'] >=s :
				return j
			if self.line[Nb-1].data['s'] < s : 
				return 0
		if s < 0.0 : 
			return Nb-1
		else:
			print "error: s not in range"
			print "STOP."
			sys.exit(1)

	def get_length(self):
		Nb=len(self.line)
		return self.line[Nb-1].data['s']


	# calculates the space charge perveance using given bunch and (Ex,Ey). Is it better to read (Ex,Ey) from the bunch? 
	def getPerveance(self,bunch,circum,emitx,emity):

		# global const 
		clight=2.988e8
		eps0=8.8542e-12
		mu0=1.0/(eps0*clight**2)
		qe=1.6022e-19
		mp=1.6605e-27
		me=9.1094e-31
		# global const 

		syncPart=bunch.getSyncParticle()
		e_kin=syncPart.kinEnergy()
		N = bunch.macroSize()

		gamma0=1.0+(e_kin*1e6*qe)/(mp*clight*clight)
		beta0=np.sqrt((gamma0*gamma0-1.0)/(gamma0*gamma0))
		print("gamma={} beta={}".format(gamma0,beta0))

		current=beta0*clight*N*qe/circum # dc current

		R=circum/(2.0*np.pi)
		Ksc=2.0*current/((4.0*np.pi*mp*eps0*clight**3)/(qe)*(beta0*gamma0)**3)
		print("perveance is {}".format(Ksc))
		return Ksc

		
	def readtwiss_teapot(self,lattice, bunch):
		
		beamline=Optics()
		matrix_lattice = TEAPOT_MATRIX_Lattice(lattice,bunch)

		(arrmuX, arrPosAlphaX, arrPosBetaX) = matrix_lattice.getRingTwissDataX()
		(arrmuY, arrPosAlphaY, arrPosBetaY) = matrix_lattice.getRingTwissDataY()

		(DispersionX, DispersionXP) = matrix_lattice.getRingDispersionDataX()
		(DispersionY, DispersionYP) = matrix_lattice.getRingDispersionDataY()

		nodes = lattice.getNodes()

		matrixNodes = matrix_lattice.getNodes()	
		idx = 0
		Md = np.identity(6) # Matrix of downstream element

		for node in nodes:
			for j in range(len(arrPosBetaX)):
				if (round(lattice.getNodePositionsDict()[node][1],4)==round(arrPosBetaX[j][0],4)):
					muX = arrmuX[j][1]
					betaX = arrPosBetaX[j][1]
					alphaX =  arrPosAlphaX[j][1]
					dx = DispersionX[j][1]
					dmux = DispersionXP[j][1]
					muY = arrmuY[j][1]
					betaY = arrPosBetaY[j][1]
					alphaY = arrPosAlphaY[j][1]
					dmuy = DispersionYP[j][1]

			if node.getType() == "quad teapot":
				k1l = node.getParam("kq")*node.getLength()
			else:


		    		if node.getType()=="multipole teapot":
			    		k1l = node.getParam("kls")[1]
				else:
					k1l = 0.0
			if node.getType() == "bend teapot":
				angle = node.getParam("theta")
			else:
				angle = 0.0
			beamline.add(1)
			j=len(beamline)-1
			beamline[j]=Twiss()
			name = node.getName()
			beamline[j].data['keyword']=name
			beamline[j].data['marker']=node.getType()
			beamline[j].data['s']=round(lattice.getNodePositionsDict()[node][1],4)
			beamline[j].data['L']=node.getLength()
			beamline[j].data['alfx']=alphaX
			beamline[j].data['alfy']=alphaY
			beamline[j].data['betx']=betaX
			beamline[j].data['bety']=betaY
			beamline[j].data['Dx']=dx
			beamline[j].data['Dpx']=dmux
			beamline[j].data['mux']=muX
			beamline[j].data['muy']=muY
			beamline[j].data['angle']=angle
			beamline[j].data['k1']=k1l

			matName = matrixNodes[idx].getName()

			M = np.identity(6)

			while name in matName and idx < len(matrixNodes)-1:
				matNode = matrixNodes[idx]
				matName = matNode.getName()
				mt = matNode.getMatrix()
				idx+=1
				for index in range(36):
					Md[index//6,index%6] = mt.get(index//6,index%6)

				M = np.dot(M[:],Md)

			beamline[j].data['map'] = M

		return beamline
#------------------------------------------------------
# Read MADX TFS file
#-------------------------------------------------------



   
#------------------------------------------------------
# Envelope solver: 
# x0, xs0, y0, ys0: initial values
# emitx/y: rms emittance
# Ksc: space charge perveance
#-------------------------------------------------------
class EnvelopeSolver:
	
	def __init__(self,beamline):
		self.beamline = beamline


	def func_odeint(self,y,s,emitx,emity,sigma_p,Ksc):
		jb=self.beamline.get_element(s)
		k1=self.beamline[jb].data['k1']
		lj=self.beamline[jb].data['L']
		anglej=self.beamline[jb].data['angle']

		f0=y[1]
		f1=-(k1/lj+(anglej/lj)**2)*y[0]+emitx**2/y[0]**3+0.5*Ksc/(y[0]+y[2])+y[4]*sigma_p**2*anglej/(y[0]*lj) 
		f2=y[3]
		f3=(k1/lj)*y[2]+emity**2/y[2]**3+0.5*Ksc/(y[0]+y[2]) # -
		f4=y[5]
		f5=-(k1/lj+(anglej/lj)**2)*y[4]+0.5*Ksc/(y[0]*(y[0]+y[2]))*y[4]+anglej/lj

		return [f0,f1,f2,f3,f4,f5]
	
	def Dfunc_odeint(self,y,s,emitx,emity,sigma_p,Ksc):
		jb=self.beamline.get_element(s)
		k1=self.beamline[jb].data['k1']
		lj=self.beamline[jb].data['L']
		anglej=self.beamline[jb].data['angle']
		a0=-(k1/lj+(anglej/lj)**2)*y[0]+emitx**2/y[0]**3+0.5*Ksc/(y[0]+y[2])+y[4]*sigma_p**2*anglej/(y[0]*lj)
		a1=-(k1/lj+(anglej/lj)**2)*y[1]-3.0*y[1]*emitx**2/y[0]**4-0.5*Ksc*(y[1]+y[3])/(y[0]+y[2])**2+y[5]*sigma_p**2*anglej/(y[0]*lj)-y[4]*y[1]*sigma_p**2*anglej/(y[0]**2*lj) 
		a2=(k1/lj)*y[2]+emity**2/y[2]**3+0.5*Ksc/(y[0]+y[2]) # -
		a3=(k1/lj)*y[3]-3.0*y[3]*emity**2/y[2]**4-0.5*Ksc*(y[1]+y[3])/(y[0]+y[2])**2 # -
		a4=-(k1/lj+(anglej/lj)**2)*y[4]+0.5*Ksc/(y[0]*(y[0]+y[2]))*y[4]+anglej/lj
		a5=-(k1/lj+(anglej/lj)**2)*y[5]+0.5*Ksc/(y[0]*(y[0]+y[2]))*y[5]-0.5*Ksc/(y[0]*(y[0]+y[2]))**2*y[4]*(y[1]*(y[0]+y[2])+y[0]*(y[1]+y[3]) )

		return [a0,a1,a2,a3,a4,a5]

	
	def envelope_odeint(self, emitx, emity, sigma_p, Ksc, x0, xs0, y0, ys0, Dx0, Dxs0):  
		Np=1000	
		Nb=len(self.beamline)
		Lb=self.beamline[Nb-1].data['s']
		s=linspace(0.0,Lb,num=Np)

#		sol=odeint(self.func_odeint,[x0,xs0,y0,ys0,Dx0,Dxs0],s,args=(emitx,emity,sigma_p,Ksc),Dfun=self.Dfunc_odeint,rtol=1.0e-12,atol=1.0e-12)
		sol=odeint(self.func_odeint,[x0,xs0,y0,ys0,Dx0,Dxs0],s,args=(emitx,emity,sigma_p,Ksc),rtol=1.0e-12,atol=1.0e-12)
		envx=sol[:,0]
		envxs=sol[:,1]
		envy=sol[:,2]
		envys=sol[:,3]
		Dx=sol[:,4]
		Dxs=sol[:,5]
		return envx,envxs,envy,envys,Dx,Dxs,s



	#------------------------------------------------------
	# Match: Periodic solution starting from MADX result
	#-------------------------------------------------------
 
	# this is the function for the root searching routine (fsolve)
 
	def func_fsolve(self,x,emitx,emity,sigma_p,Ksc):
		envx,envxs,envy,envys,Dx,Dxs,s = self.envelope_odeint(emitx,emity,sigma_p,Ksc,x[0],x[1],x[2],x[3],x[4],x[5])
		Nb=len(envx)
		return [envx[Nb-1]-x[0],envxs[Nb-1]-x[1],envy[Nb-1]-x[2],envys[Nb-1]-x[3],Dx[Nb-1]-x[4],Dxs[Nb-1]-x[5]] 

	# root searching using fsolve and initial values from MADX
	# returns matched envelopes

	def match_root(self, emitx, emity, sigma_p, Ksc):
		Nb=len(self.beamline)
		# start values
		x0=sqrt(self.beamline[Nb-1].data['betx']*emitx)
		gamx=(1.0+(self.beamline[Nb-1].data['alfx'])**2)/self.beamline[Nb-1].data['betx']
		xs0=-copysign(sqrt(gamx*emitx),self.beamline[Nb-1].data['alfx'])
		y0=sqrt(self.beamline[Nb-1].data['bety']*emity)
		gamy=(1.0+(self.beamline[Nb-1].data['alfy'])**2)/self.beamline[Nb-1].data['bety']
		ys0=-copysign(sqrt(gamy*emity),self.beamline[Nb-1].data['alfy'])
		Dx0=self.beamline[Nb-1].data['Dx']
		Dxs0=self.beamline[Nb-1].data['Dpx']
		# solver
		sol = root(self.func_fsolve,[x0,xs0,y0,ys0,Dx0,Dxs0], args=(emitx,emity,sigma_p,Ksc),method='hybr')
		x0=sol.x[0]
		xs0=sol.x[1]
		y0=sol.x[2]
		ys0=sol.x[3]
		Dx0=sol.x[4]
		Dxs0=sol.x[5]
		envx,envxs,envy,envys,Dx,Dxs,s = self.envelope_odeint(emitx,emity,sigma_p,Ksc,x0,xs0,y0,ys0,Dx0,Dxs0)
		return envx, envxs, envy, envys, Dx, Dxs, s


	# returns the matchted twiss parameter at cell entrance 

	def match_twiss(self, emitx, emity, sigma_p, Ksc):
		Nb=len(self.beamline)
		# start values
		x0=sqrt(self.beamline[Nb-1].data['betx']*emitx)
		gamx=(1.0+(self.beamline[Nb-1].data['alfx'])**2)/self.beamline[Nb-1].data['betx']
		xs0=-copysign(sqrt(gamx*emitx),self.beamline[Nb-1].data['alfx'])
		y0=sqrt(self.beamline[Nb-1].data['bety']*emity)
		gamy=(1.0+(self.beamline[Nb-1].data['alfy'])**2)/self.beamline[Nb-1].data['bety']
		ys0=-copysign(sqrt(gamy*emity),self.beamline[Nb-1].data['alfy'])
		Dx0=self.beamline[Nb-1].data['Dx']
		Dxs0=self.beamline[Nb-1].data['Dpx']
		# solver
		sol = root(self.func_fsolve, [x0,xs0,y0,ys0,Dx0,Dxs0], args=(self.beamline,emitx,emity,sigma_p,Ksc),method='krylov')
		x0=sol.x[0]
		xs0=sol.x[1]
		y0=sol.x[2]
		ys0=sol.x[3]
		Dx0=sol.x[4]
		Dxs0=sol.x[5]
		return x0**2/emitx,y0**2/emity,-copysign(sqrt(x0**2*xs0**2/emitx**2),xs0),-copysign(sqrt(y0**2*ys0**2/emity**2),ys0), Dx0, Dxs0


	#------------------------------------------------------
	# Smooth focusing 
	#-------------------------------------------------------


	def func_smooth(self,x,phase0x,phase0y,length,emitx,emity,Ksc):
		kx=(phase0x/length)**2
		ky=(phase0y/length)**2
		return[emitx**2/x[0]**3-kx*x[0]+0.5*Ksc/(x[0]+x[1]),emity**2/x[1]**3-ky*x[1]+0.5*Ksc/(x[0]+x[1])]


	def match_smooth(self,phase0x,phase0y,length,emitx,emity,Ksc):
		kx=(phase0x/length)**2
		ky=(phase0y/length)**2

		x0=(emitx**2/kx)**(1.0/4.0) 
		y0=(emity**2/ky)**(1.0/4.0) 

		sol = root(self.func_smooth,[x0,y0],args=(phase0x,phase0y,length,emitx,emity,Ksc),method='hybr')
		
		return sol.x[0]**2/emitx,sol.x[1]**2/emity     # beta functions


	#------------------------------------------------------
	# Calculate phase advance for given envelopes
	#-------------------------------------------------------

	def phase_advance(self,envx,envy,Dx,emitx,emity,sigma_p,s):
		Np=len(s)
		phasex=0.0
		phasey=0.0
		ds=s[1]-s[0]
		for j in xrange(0,Np):
			phasex+=ds*emitx/(envx[j]**2-(Dx[j]*sigma_p)**2)
			phasey+=ds*emity/envy[j]**2
		return phasex, phasey	

	# analytic phase advance depression
	# lc: length of the cell

	def phase_analytic(self,emitx,emity,Ksc,lc):
		return 0.5*Ksc*lc/(4.0*emitx), 0.5*Ksc*lc/(4.0*emity) 


	#------------------------------------------------------
	# Entropy growth rate: pre-factor
	#-------------------------------------------------------

	def entropy_rate(self,envx,envy,emitx,emity,s,beta0):
		Np=len(s)
		ratet=0.0
		ds=s[1]-s[0]
		for j in xrange(0,Np):
			Tx=envx[j]**2/emitx**2
			Ty=envy[j]**2/emity**2
			ratet+=ds/(beta0*c)*0.5*(Tx-Ty)**2/(Tx*Ty)
		return ratet

	#------------------------------------------------------
	# Kick-by-kick for thin lattices
	#------------------------------------------------------

	def twiss_transport(self, M,tw):
		gamx0=(1.0+tw.data['alfx']**2)/tw.data['betx']
		betx=M[0,0]**2*tw.data['betx']-2.0*M[0,0]*M[0,1]*tw.data['alfx']+M[0,1]**2*gamx0
		alfx=-M[0,0]*M[1,0]*tw.data['betx']+(M[0,0]*M[1,1]+M[0,1]*M[1,0])*tw.data['alfx']-M[1,1]*M[0,1]*gamx0
		gamx=M[1,0]**2*tw.data['betx']-2.0*M[1,1]*M[1,0]*tw.data['alfx']+M[1,1]**2*gamx0
		gamy0=(1.0+tw.data['alfy']**2)/tw.data['bety']
		bety=M[2,2]**2*tw.data['bety']-2.0*M[2,2]*M[2,3]*tw.data['alfy']+M[2,3]**2*gamy0
		alfy=-M[2,2]*M[3,2]*tw.data['bety']+(M[2,2]*M[3,3]+M[2,3]*M[3,2])*tw.data['alfy']-M[3,3]*M[2,3]*gamy0
		gamy=M[3,2]**2*tw.data['bety']-2.0*M[3,3]*M[3,2]*tw.data['alfy']+M[3,3]**2*gamy0
		Dx=M[0,0]*tw.data['Dx']+M[0,1]*tw.data['Dpx']+M[0,5]
		Dpx=M[1,0]*tw.data['Dx']+M[1,1]*tw.data['Dpx']+M[1,5]

		tw.data['betx'] = betx 
		tw.data['alfx'] = alfx
		tw.data['bety'] = bety
		tw.data['alfy'] = alfy
		tw.data['Dx'] = Dx
		tw.data['Dpx'] = Dpx 

		return tw



	def twiss_evolution(self,tw,Ksc,emitx,emity,sigma_p):
		Nelements=len(self.beamline)
		tw0=deepcopy(tw)
		twiss_vec=zeros((Nelements,6))


		Msc=np.identity(6)

		ax_lens=0.0
		ay_lens=0.0
		for j, b in enumerate(self.beamline):
			ds=b.data['L']

			tw0 =self.twiss_transport(b.data['map'],tw0) # transport without space charge

			ax=sqrt(tw0.data['betx']*emitx + (tw0.data['Dx']*sigma_p)**2)				
			ay=sqrt(tw0.data['bety']*emity)

			kick_gradient_x=0.5*Ksc/(ax*(ax+ay))*ds
			kick_gradient_y=0.5*Ksc/(ay*(ax+ay))*ds

			Msc[0,0]=1.0
			Msc[0,1]=0.0
			Msc[1,0]=kick_gradient_x
			Msc[1,1]=1.0
			Msc[2,2]=1.0
			Msc[2,3]=0.0
			Msc[3,2]=kick_gradient_y
			Msc[3,3]=1.0

	
			# update tw0
			tw0=self.twiss_transport(Msc,tw0)     # space charge kick
			twiss_vec[j,0]=tw0.data['betx']
			twiss_vec[j,1]=tw0.data['alfx']
			twiss_vec[j,2]=tw0.data['bety']
			twiss_vec[j,3]=tw0.data['alfy']
			twiss_vec[j,4]=tw0.data['Dx']
			twiss_vec[j,5]=tw0.data['Dpx']

		return twiss_vec

	def match_twiss_matrix(self,emitx,emity,sigma_p,Ksc):
		# init tw
		tw=deepcopy(self.beamline[-1])
		bline = self.beamline[-1]
	   	# start value
		x0=sqrt(bline.data['betx']*emitx)
		gamx=(1.0+(bline.data['alfx'])**2)/bline.data['betx']
		xs0=-copysign(sqrt(gamx*emitx),bline.data['alfx'])
		y0=sqrt(bline.data['bety']*emity)
		gamy=(1.0+(bline.data['alfy'])**2)/bline.data['bety']
		ys0=-copysign(sqrt(gamy*emity),bline.data['alfy'])
		Dx0=bline.data['Dx']
		Dxs0=bline.data['Dpx']
		# solver
		print("[x0,xs0,y0,ys0,Dx0,Dxs0] = {}".format([x0,xs0,y0,ys0,Dx0,Dxs0]))
		sol=root(self.func_fsolve_matrix, [x0,xs0,y0,ys0,Dx0,Dxs0], args=(Ksc,emitx,emity,sigma_p),method='lm')

		if not sol.success:
			print("Sollution not found try again")
	
		# upd twiss
		tw.data['betx']=sol.x[0]
		tw.data['alfx']=sol.x[1]
		tw.data['bety']=sol.x[2]
		tw.data['alfy']=sol.x[3]
		tw.data['Dx']=sol.x[4]
		tw.data['Dpx']=sol.x[5]

		twiss_vec=self.twiss_evolution(tw,Ksc,emitx,emity,sigma_p)
		return twiss_vec


	def func_fsolve_matrix(self,x,Ksc,emitx,emity,sigma_p):
		tw=deepcopy(self.beamline[-1])
		tw.data['betx']=x[0]
		tw.data['alfx']=x[1]
		tw.data['bety']=x[2]
		tw.data['alfy']=x[3]
		tw.data['Dx']=x[4]
		tw.data['Dpx']=x[5]

		twiss_vec=self.twiss_evolution(tw,Ksc,emitx,emity,sigma_p)
		return [twiss_vec[-1,0]-x[0],twiss_vec[-1,1]-x[1],twiss_vec[-1,2]-x[2],twiss_vec[-1,3]-x[3],twiss_vec[-1,4]-x[4],twiss_vec[-1,5]-x[5]]












