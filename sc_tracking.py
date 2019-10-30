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
from spacecharge import SpaceChargeCalc2p5D

from orbit.bunch_generators import TwissContainer, TwissAnalysis
from orbit.bunch_generators import GaussDist3D, KVDist3D

from orbit.errors import AddErrorNode,AddErrorSet

#---------------------------------------------Bunch init---------------------------------------------
print "Start."

emitx = 1e-7
emity = 1e-7
sigma_p = 1e-3

b = Bunch()
total_macroSize=1e+11
b.macroSize(total_macroSize)

energy = 1
syncPart=b.getSyncParticle()
syncPart.kinEnergy(energy)
EE=syncPart.kinEnergy()
print(EE)
#---------------------------------------------Make a Teapot Lattice----------------------------------
print "Generate Lattice."

lattice = TEAPOT_Lattice("lattice")
#lattice.readMADX("fodo_thin.seq","fodo")
lattice.readMADX("FODO_x6.seq","fodo")
#lattice.readMADX("cryring.madx","cryring")
#lattice.readMADX("sis100_full_thin_fix.seq","sis100ring")
lattice.setUseRealCharge(useCharge = 1)

#---------------------------------------------SPLIT LONG ELEMENTS--------------------------------------

for node in lattice.getNodes():
	if node.getLength() > 1.0:
		node.setnParts(int(node.getLength()//1+1))

#------------------------------------------------------------------------------------------------------

matrix_lattice = TEAPOT_MATRIX_Lattice(lattice,b)

TwissDataX0,TwissDataY0 = matrix_lattice.getRingTwissDataX(),matrix_lattice.getRingTwissDataY()

beta_x0 = np.transpose(TwissDataX0[-1])
beta_y0 = np.transpose(TwissDataY0[-1])

#----------------------------Add Space Charge nodes----------------------------------------------------

#---------------------------------------------------- ERRORS-------------------------------------------
# WE INTRODUCE MISALIGNEMENT IN THE QUADRUPOLES; dx, dy = HOR AND VER DISPLACEMENT OF QUADRUPOLES

setDict = {}
paramsDict = {}
positioni = 0.0
positionf = lattice.getLength()

setDict["elementtype"] = "quad"
setDict["ringline"] = "ring"

paramsDict["errtype"]  = "FieldError"
paramsDict["subtype"]  = "QuadField"

paramsDict["sample"]      = "Gaussian"
paramsDict["fracerr"]      = 0.01
paramsDict["mean"]        = 0.0
paramsDict["sigma"]       = 1.0

#ESet  = AddErrorSet(lattice, positioni, positionf, setDict, paramsDict, seed_value=50)
ESet  = AddErrorSet(lattice, positioni, positionf, setDict, paramsDict) # Random
print "FIELD ERRORS IN THE QUADRUPOLES/MULTIPOLES ARE APPLIED"
#------------------------------------------------------------------------------------------------------

matrix_lattice = TEAPOT_MATRIX_Lattice(lattice,b)

TwissDataX1,TwissDataY1 = matrix_lattice.getRingTwissDataX(),matrix_lattice.getRingTwissDataY()

beta_x1 = np.transpose(TwissDataX1[-1])
beta_y1 = np.transpose(TwissDataY1[-1])

beat_x0=100*(beta_x1[1]/beta_x0[1] - 1)	
beat_y0=100*(beta_y1[1]/beta_y0[1] - 1)	


bc = betaCorrection(lattice,b)
bc.correction()
bc.getSolution()


matrix_lattice = TEAPOT_MATRIX_Lattice(lattice,b)

TwissDataX2,TwissDataY2 = matrix_lattice.getRingTwissDataX(),matrix_lattice.getRingTwissDataY()

betaX0,betaY0 = bc.getBetaBpm(TwissDataX0,TwissDataY0)
betaX1,betaY1 = bc.getBetaBpm(TwissDataX1,TwissDataY1)
betaX2,betaY2 = bc.getBetaBpm(TwissDataX2,TwissDataY2)

beta_x2 = np.transpose(TwissDataX2[-1])
beta_y2 = np.transpose(TwissDataY2[-1])

beat_x1=100*(beta_x2[1]/beta_x0[1] - 1)	
beat_y1=100*(beta_y2[1]/beta_y0[1] - 1)	



#----------------------------Add Space Charge nodes----------------------------------------------------
'''
sc_path_length_min = 0.000001
sizeX = 32   #number of grid points in horizontal direction
sizeY = 32   #number of grid points in vertical direction
sizeZ = 16   #number of longitudinal slices

calc2p5d = SpaceChargeCalc2p5D(sizeX,sizeY,sizeZ)
scLatticeModifications.setSC2p5DAccNodes(lattice, sc_path_length_min, calc2p5d)

print "SC nodes appied to the lattice"
'''
#--------------------------------------------------------------------------------
'''
n_particles = 1000
twissX = TwissContainer(alpha = arrPosAlphaX0[0][1], beta = arrPosBetaX0[0][1], emittance = emitx)
twissY = TwissContainer(alpha = arrPosAlphaY0[0][1], beta = arrPosBetaY0[0][1], emittance = emity)
twissZ = TwissContainer(alpha = 0., beta = 100000000., emittance = sigma_p)
dist = GaussDist3D(twissX,twissY,twissZ)
#dist = KVDist3D(twissX,twissY,twissZ)


for i in range(n_particles):
	(x,xp,y,yp,z,zp) = dist.getCoordinates()
	b.addParticle(x, xp, y, yp, z, zp)

print "Bunch Generated"
#--------------------------------Tracking-----------------------------------------

x_space0 =[ (b.x(i), b.px(i)) for i in range(n_particles)]

n_turns = 100
for i in range(2,n_turns):
	lattice.trackBunch(b)

x_space =[ (b.x(i), b.px(i)) for i in range(n_particles)]
'''
#--------------------------------------------------------------------------------

plt.figure()
plt.subplot(211)
#plt.plot(beta_x1[0],beta_x1[1], label="horizontal beta-beat")
plt.plot(beta_x0[0],beat_x0, label="horizontal beta-beat before corr")
plt.plot(beta_x0[0],beat_x1, label="horizontal beta-beat corrected")

[plt.axhline(x,color = "blue", ls=":") for x in (betaX1/betaX0-1)*100]
[plt.axhline(x,color = "red") for x in (betaX2/betaX0-1)*100]

plt.xlabel(r'$s$ [m]')
plt.ylabel(r'$\beta_x$ [m]')

plt.subplot(212)
#plt.plot(beta_y1[0],beta_y1[1], label="vertical beta-beat")
plt.plot(beta_y0[0],beat_y0, label="vertical beta-beat")
plt.plot(beta_y0[0],beat_y1, label="vertical beta-beat corrected")

[plt.axhline(y,color = "blue", ls=":") for y in (betaY1/betaY0-1)*100]
[plt.axhline(y,color = "red") for y in (betaY2/betaY0-1)*100]

plt.xlabel(r'$s$ [m]')
plt.ylabel(r'$\beta_y$ [m]')
plt.subplots_adjust(bottom=0.15,left=0.15,hspace=0.5)

#--------------------------------------------------------------------------------
'''
plt.figure()
plt.scatter(*np.transpose(x_space0))
plt.scatter(*np.transpose(x_space))
plt.xlabel('x, m')
plt.ylabel('$p_x$, rad')
plt.grid(True, ls = '--')
plt.title('Phase Space in (x, $p_x$) plane, N = {} particles, turn = 0'.format(n_particles))
'''
#--------------------------------------------------------------------------------
plt.show()
#--------------------------------------------------------------------------------


