#------------------------------------------------
#pyORBIT error and correction example
#------------------------------------------------

import sys
import math
import matplotlib.pyplot as plt
import numpy as np

import orbit_mpi

from orbit.teapot import teapot
from orbit.teapot import TEAPOT_Lattice

from orbit.lattice import AccLattice, AccNode, AccActionsContainer
from bunch import Bunch

from orbit.matching import Optics, EnvelopeSolver

from orbit.teapot import TEAPOT_MATRIX_Lattice

from orbit.utils.orbit_mpi_utils import bunch_pyorbit_to_orbit

# ATTENTION !!! The python packages numpy and scipy are required
from orbit.orbit_correction import orbit, correction   
from orbit.beta_correction import betaCorrection   

from spacecharge import Boundary2D
from orbit.space_charge.sc2p5d import scAccNodes, scLatticeModifications
#from spacecharge import SpaceChargeCalc2p5D
from spacecharge import SpaceChargeForceCalc2p5D

from orbit.bunch_generators import TwissContainer, TwissAnalysis
from orbit.bunch_generators import GaussDist3D, KVDist3D

from orbit.errors import AddErrorNode,AddErrorSet

from bunch import BunchTwissAnalysis

from orbit.diagnostics import TeapotTuneAnalysisNode,addTeapotDiagnosticsNode

#---------------------------------------------Bunch init---------------------------------------------
print "Start."

emitx = 1e-9
emity = 1e-9
sigma_p = 1e-5

b = Bunch()
total_macroSize=1e+11
b.macroSize(total_macroSize)

lostbunch = Bunch()
trackDict = {}
trackDict["lostbunch"]=lostbunch
trackDict["bunch"]= b
lostbunch.addPartAttr("LostParticleAttributes")

energy = 1
syncPart=b.getSyncParticle()
syncPart.kinEnergy(energy)
EE=syncPart.kinEnergy()
print(EE)
#---------------------------------------------Make a Teapot Lattice----------------------------------
print "Generate Lattice."

lattice = TEAPOT_Lattice("lattice")

lattice.readMADX("FODO_x6.seq","fodo")
lattice.setUseRealCharge(useCharge = 1)


#------------------------------------------------------------------------------------------------------

f = open("tunes.dat")
tunes=[[float(x) for x in line.split(",")] for line in f]

def setQuads(kd,kf):
	for elem in lattice.getNodes():
		name,type = elem.getName(),elem.getType()
		if type=="quad teapot" and "d" in name:
			elem.setParam("kq", kd)
		if type=="quad teapot" and "f" in name:
			elem.setParam("kq", kf)




#--------------------------------------------------------------------------------
matrix_lattice = TEAPOT_MATRIX_Lattice(lattice,b)

TwissDataX0,TwissDataY0 = matrix_lattice.getRingTwissDataX(),matrix_lattice.getRingTwissDataY()

beta_x0 = np.transpose(TwissDataX0[-1])
beta_y0 = np.transpose(TwissDataY0[-1])

#---------------------------------------------------- ERRORS-------------------------------------------

setDict = {}
paramsDict = {}
positioni = 0.0
positionf = lattice.getLength()

setDict["elementtype"] = "quad"
setDict["ringline"] = "ring"

paramsDict["errtype"]  = "FieldError"
paramsDict["subtype"]  = "QuadField"

paramsDict["sample"]      = "Gaussian"
paramsDict["fracerr"]      = 0.012
paramsDict["mean"]        = 0.0
paramsDict["sigma"]       = 1.0

ESet  = AddErrorSet(lattice, positioni, positionf, setDict, paramsDict, seed_value=40)
#ESet  = AddErrorSet(lattice, positioni, positionf, setDict, paramsDict) # Random
print "FIELD ERRORS IN THE QUADRUPOLES/MULTIPOLES ARE APPLIED"

matrix_lattice = TEAPOT_MATRIX_Lattice(lattice,b)
TwissDataX1,TwissDataY1 = matrix_lattice.getRingTwissDataX(),matrix_lattice.getRingTwissDataY()

#--------------------------------------------------------------------------------

n_particles = 5000
twissZ = TwissContainer(alpha = 0., beta = 100000000., emittance = sigma_p)
b.macroSize(total_macroSize/n_particles)


#---------------------------------------------SPLIT LONG ELEMENTS--------------------------------------

print("split long elements")
for node in lattice.getNodes():
	if node.getLength() > 1.0:
		node.setnParts(int(node.getLength()//1+1))


#----------------------------Add Space Charge nodes----------------------------------------------------

sc_path_length_min = 0.000001
sizeX = 32   #number of grid points in horizontal direction
sizeY = 32   #number of grid points in vertical direction
sizeZ = 16   #number of longitudinal slices

calc2p5d = SpaceChargeForceCalc2p5D(sizeX,sizeY,sizeZ)
scLatticeModifications.setSC2p5DAccNodes(lattice, sc_path_length_min, calc2p5d)

print "SC nodes appied to the lattice"

#--------------------------------Tracking-----------------------------------------

bunchtwissanalysis = BunchTwissAnalysis()
n_turns = 250
strIn = "_{}_{}_{}_{:.0e}_{}".format("gauss",n_particles,emitx,total_macroSize,n_turns)

w = open("tune_scan_corr_tmp1{}.dat".format(strIn),"w")
#w = open("tune_scan_corrected{}.dat".format(strIn),"w")
w.write("qx,qy,emitx,emity\n")

for k,wp in enumerate(tunes[10:]):
	setQuads(wp[0],wp[1])

	bc = betaCorrection(lattice,b)
#	if wp[2]==1.95 or wp[3]== 1.95:
#		print("warn")
#		bc.correction(rhobeg=1e-3, A=0.001)
#	else:
	bc.correction(rhobeg=5e-3, A=0.01)

	bc.getSolution()
	bc.fixTunes(rhobeg=5e-5,qx0=wp[2],qy0=wp[3])

	if bc.success == True:
		matrix_lattice = TEAPOT_MATRIX_Lattice(lattice,b)
		TwissDataX,TwissDataY = matrix_lattice.getRingTwissDataX(),matrix_lattice.getRingTwissDataY()
		print("{} {}".format(TwissDataX[0][-1][1],TwissDataY[0][-1][1]))

		twissX = TwissContainer(alpha = TwissDataX[1][0][1], beta = TwissDataX[2][0][1], emittance = emitx)
		twissY = TwissContainer(alpha = TwissDataY[1][0][1], beta = TwissDataY[2][0][1], emittance = emity)
		dist = GaussDist3D(twissX,twissY,twissZ)

		for i in range(n_particles):
			b.addParticle(*dist.getCoordinates())


		print("tracking")
		for i in range(n_turns):
			lattice.trackBunch(b,trackDict)
		bunchtwissanalysis.analyzeBunch(b)
		Ex = bunchtwissanalysis.getEmittance(0)
		Ey = bunchtwissanalysis.getEmittance(1)
		w.write("{},{},{},{}\n".format(TwissDataX[0][-1][1],TwissDataY[0][-1][1],Ex,Ey))

		b.deleteAllParticles()

		print("we are at the {}-th position".format(k))
	else:
		print("the problem occured: {}-th position".format(k))


w.close()

