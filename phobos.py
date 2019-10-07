#------------------------------------------------
#pyORBIT error and correction example
#------------------------------------------------

import sys
import math
import matplotlib.pyplot as plt
import numpy as np
#import pandas as pd

import orbit_mpi

from orbit.teapot import teapot
from orbit.teapot import TEAPOT_Lattice
from orbit.lattice import AccLattice, AccNode, AccActionsContainer
from bunch import Bunch
from orbit.errors import AddErrorNode
from orbit.errors import AddErrorSet

from orbit.matching import Optics, EnvelopeSolver


from orbit.teapot import TEAPOT_MATRIX_Lattice


from orbit.utils.orbit_mpi_utils import bunch_pyorbit_to_orbit

# ATTENTION !!! The python packet numpy and scipy are required
from orbit.orbit_correction import orbit, correction   


from orbit.bunch_generators import TwissContainer, TwissAnalysis
from orbit.bunch_generators import GaussDist3D, KVDist3D

from spacecharge import Boundary2D
from orbit.space_charge.sc2p5d import scAccNodes, scLatticeModifications
from spacecharge import SpaceChargeCalc2p5D



print "Start."
#---------------------------------------------Bunch init---------------------------------------------
b = Bunch()
total_macroSize=4e+10
b.macroSize(total_macroSize)

energy = 2  
syncPart=b.getSyncParticle()
syncPart.kinEnergy(energy)
EE=syncPart.kinEnergy()
print(EE)
#---------------------------------------------Bunch init---------------------------------------------

print "Generate Lattice."
#---------------------------------------------Make a Teapot Lattice----------------------------------
lattice = TEAPOT_Lattice("lattice")
lattice.readMADX("fodo_thin.seq","fodo")
lattice.setUseRealCharge(useCharge = 1)

#-----------------------------------------------
matrix_lattice = TEAPOT_MATRIX_Lattice(lattice,b)
(muX0, arrPosAlphaX0, arrPosBetaX0) = matrix_lattice.getRingTwissDataX()
(muY0, arrPosAlphaY0, arrPosBetaY0) = matrix_lattice.getRingTwissDataY()

beamline0 = Optics().readtwiss_teapot(lattice, b)
solve = EnvelopeSolver(beamline0)
twiss0 = solve.match_twiss_matrix(12.5e-6,12.5e-6,0.0,0.0) 
#-----------------------------------------------

#---------------------------------------------Make a Teapot Lattice----------------------------------

print "INTRODUCE MISALIGNEMENT IN THE QUADRUPOLES"
#---------------------------------------------ORBIT ERRORS-------------------------------------------
# WE INTRODUCE MISALIGNEMENT IN THE QUADRUPOLES; dx, dy = HOR AND VER DISPLACEMENT OF QUADRUPOLES

setDict = {}
paramsDict = {}
positioni = 0.0
positionf = lattice.getLength()

paramsDict["errtype"]  = "StraightError"
paramsDict["subtype"]  = "TransDisp"
paramsDict["sample"]      = "ReadFromFile"
paramsDict["filename"] = "/home/dmitrii/py-orbit/examples/AccLattice_Tests/err.txt"

setDict["elementtype"] = "mult"
setDict["ringline"] = "ring"

print "INTRODUCE FIELD ERRORS IN THE QUADRUPOLES"

#ESet  = AddErrorSet(lattice, positioni, positionf, setDict, paramsDict, seed_value=50)

paramsDict["errtype"]  = "FieldError"
paramsDict["subtype"]  = "MultipoleField"
paramsDict["errtype"]  = "FieldError"

#paramsDict["sample"]      = "Gauss"
#paramsDict["fracerr"]      = 0.001
#paramsDict["maximum"]        = 0.001
#paramsDict["minimum"]       = 0.0

#ESet  = AddErrorSet(lattice, positioni, positionf, setDict, paramsDict, seed_value=50)
#ESet  = AddErrorSet(lattice, positioni, positionf, setDict, paramsDict) # Random
#---------------------------------------------ORBIT ERRORS-------------------------------------------


#tmp2=[node.getParam("kls")[1] for node in lattice.getNodes() if node.getType()=="multipole teapot"]
#print(np.array(tmp2)-np.array(tmp1))


#---------------------------------------------CALCULATE DISTORTED BETA (BEAT) ------------------------------
print "CALCULATE DISTORTED BETA"

matrix_lattice = TEAPOT_MATRIX_Lattice(lattice,b)
(muX, arrPosAlphaX, arrPosBetaX) = matrix_lattice.getRingTwissDataX()
(muY, arrPosAlphaY, arrPosBetaY) = matrix_lattice.getRingTwissDataY()

beta_x0 = np.transpose(arrPosBetaX0)
beta_x = np.transpose(arrPosBetaX)
beta_y0 = np.transpose(arrPosBetaY0)
beta_y = np.transpose(arrPosBetaY)

beat_x=100*(beta_x[1]/beta_x0[1] - 1)	
beat_y=100*(beta_y[1]/beta_y0[1] - 1)	


#===================== m a t c h i n g =====================#

beamline = Optics().readtwiss_teapot(lattice, b)

beamline.print_line()

#N = total_macroSize
emitx = 12.5e-6     # emittance_x
emity = 12.5e-6     # emittance_y
sigma_p = 1.0e-3    # rms momentum spread


#===================== m a t c h i n g =====================#

solve = EnvelopeSolver(beamline)
Ksc = Optics().getPerveance(b,beamline[-1].data['s'],emitx,emity)

twiss = solve.match_twiss_matrix(emitx,emity,0.0,0.0) 
twiss_sc = solve.match_twiss_matrix(emitx,emity,sigma_p,Ksc)

s_tw = [x.data['s'] for x in beamline]

#===================== m a t c h i n g =====================#


n_particles = 1000
twissX = TwissContainer(alpha = arrPosAlphaX[0][1], beta = arrPosBetaX[0][1], emittance = emitx)
#twissX = TwissContainer(alpha = twiss[0,1], beta = twiss[0,0], emittance = emitx)
twissY = TwissContainer(alpha = arrPosAlphaY[0][1], beta = arrPosBetaY[0][1], emittance = emity)
twissZ = TwissContainer(alpha = 0., beta = 100000000., emittance = 0.001)
dist = GaussDist3D(twissX,twissY,twissZ)
dist = KVDist3D(twissX,twissY,twissZ)


for i in range(n_particles):
	(x,xp,y,yp,z,zp) = dist.getCoordinates()
	b.addParticle(x, xp, y, yp, z, zp)


x_space0 =[ (b.x(i), b.px(i)) for i in range(n_particles)]



print "Bunch Generated."

b.dumpBunch("bunch_2p5D_turn0.dat")

#-----------------------------------
# Tracking
#-----------------------------------

paramsDict = {}
paramsDict["bunch"]= b


n_turns = 2
for i in range(n_turns):
	lattice.trackBunch(b, paramsDict)

b.dumpBunch("bunch_turn" + str(n_turns) + ".dat")

#-----------------------------------
# Add Tune Analysis node
#-----------------------------------

#tunes = TeapotTuneAnalysisNode("tune_analysis")
#tunes.assignTwiss(10.207, 0.0469, -0.05, 0.0061, 10.639, 0.056)
#addTeapotDiagnosticsNode(lattice, 0, tunes)

#-----------------------------------
# Add Space Charge nodes
#-----------------------------------

sc_path_length_min = 0.00000001
sizeX = 32   #number of grid points in horizontal direction
sizeY = 32   #number of grid points in vertical direction
sizeZ = 16   #number of longitudinal slices

calc2p5d = SpaceChargeCalc2p5D(sizeX,sizeY,sizeZ)
scLatticeModifications.setSC2p5DAccNodes(lattice, sc_path_length_min, calc2p5d)



n_turns = 2
for i in range(n_turns):
	lattice.trackBunch(b, paramsDict)
b.dumpBunch("bunch_2p5D_turn" + str(n_turns) + ".dat")

n_turns = 100
for i in range(2,n_turns):
	lattice.trackBunch(b, paramsDict)

x_space =[ (b.x(i), b.px(i)) for i in range(n_particles)]
b.dumpBunch("bunch_2p5D_turn" + str(n_turns) + ".dat")

print("STOP.")


plt.figure()

plt.plot(beta_x[0],beat_x, label="horizontal beta-beat")
plt.plot(beta_y[0],beat_y, label="vertical beta-beat")

plt.plot(s_tw,100*(twiss[:,0]/twiss0[:,0]-1),color='k', ls=":",label=r'beatx SC = 0')
plt.plot(s_tw,100*(twiss_sc[:,0]/twiss0[:,0]-1),color='blue', ls=":",label=r'beatx non-zero SC')

plt.plot(s_tw,100*(twiss[:,2]/twiss0[:,2]-1),color='k', ls=":",label=r'beaty SC = 0')
plt.plot(s_tw,100*(twiss_sc[:,2]/twiss0[:,2]-1),color='blue', ls=":",label=r'beaty non-zero SC')

plt.xlabel(r'$s$ [m]')
plt.ylabel(r'$\beta-beat$ [/%]')
plt.legend(loc=0)
plt.grid()

#--------------------------------------------------------------------------------
plt.figure()
plt.subplot(211)
plt.scatter(s_tw,twiss[:,0],color='b',label=r'$\beta_x$')
plt.scatter(s_tw,twiss_sc[:,0],color='r',label=r'$\beta_x$')
plt.plot(beta_x[0],beta_x[1], label="horizontal beta-beat")

plt.xlabel(r'$s$ [m]')
plt.ylabel(r'$\beta_x$ [m]')
#plt.legend()

plt.subplot(212)
plt.scatter(s_tw,twiss[:,2],color='b',label=r'$\beta_x$')
plt.scatter(s_tw,twiss_sc[:,2],color='r',label=r'$\beta_x$', marker="x")
plt.plot(beta_y[0],beta_y[1], label="vertical beta-beat")

plt.xlabel(r'$s$ [m]')
plt.ylabel(r'$\beta_y$ [m]')
#plt.legend()
plt.subplots_adjust(bottom=0.15,left=0.15,hspace=0.5)

#--------------------------------------------------------------------------------

plt.figure()
plt.scatter(*np.transpose(x_space0))
plt.scatter(*np.transpose(x_space))
plt.xlabel('y, m')
plt.ylabel('$p_y$, rad')
plt.grid(True, ls = '--')
plt.title('Phase Space in (y, $p_y$) plane, N = {} particles, turn = 0'.format(n_particles))
plt.show()

