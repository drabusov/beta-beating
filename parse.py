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




print "Start."
#---------------------------------------------Bunch init---------------------------------------------
b = Bunch()
total_macroSize=1e+11
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
#lattice.readMADX("fodo_thin.seq","fodo")
lattice.readMADX("cryring.madx","cryring")
#lattice.readMADX("sis100_full_thin_fix.seq","sis100ring")
lattice.setUseRealCharge(useCharge = 1)
#-----------------------------------------------
print("check the quads")
for node in lattice.getNodes():
#	print(node.getName())
	if node.getType() == "quad teapot":
		print("parameters. kq {} length {}".format(node.getParam("kq"),node.getLength()))

print(lattice.getLength())
#-----------------------------------------------
matrix_lattice = TEAPOT_MATRIX_Lattice(lattice,b)
(muX0, arrPosAlphaX0, arrPosBetaX0) = matrix_lattice.getRingTwissDataX()
(muY0, arrPosAlphaY0, arrPosBetaY0) = matrix_lattice.getRingTwissDataY()

#print("tunes:\nQx={}\nQy".format(muX0[-1]/2/np.pi,muY0,[-1]/2/np.pi))
print("tunes:\nQx={}\nQy={}".format(muX0[-1][1],muY0[-1][1]))

svar,betax =np.transpose(arrPosBetaX0)
s_inter = np.linspace(svar[0], svar[-1], len(svar)*2)
betax_inter = np.interp(s_inter, svar, betax, period = lattice.getLength()/2)

#--------------------------------------------------------------------------------
plt.figure()
plt.subplot(211)
plt.scatter(svar,betax, label="horizontal beta")
plt.plot(s_inter,betax_inter, label="horizontal beta")
plt.xlabel(r'$s$ [m]')
plt.ylabel(r'$\beta_x$ [m]')
#plt.legend()

plt.subplot(212)
plt.plot(*np.transpose(arrPosBetaY0), label="vertical beta")

plt.xlabel(r'$s$ [m]')
plt.ylabel(r'$\beta_y$ [m]')
#plt.legend()
plt.subplots_adjust(bottom=0.15,left=0.15,hspace=0.5)

plt.show()
#--------------------------------------------------------------------------------


