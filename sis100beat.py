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

from orbit.beta_correction import betaCorrection

from orbit.errors import AddErrorNode,AddErrorSet

#---------------------------------------------Bunch init---------------------------------------------
print "Start."

b = Bunch()
energy = 1
syncPart=b.getSyncParticle()
syncPart.kinEnergy(energy)

#---------------------------------------------Make a Teapot COLD Lattice ---------------

print "Generate ideal cold Lattice."
lattice = TEAPOT_Lattice("lattice")
lattice.readMADX("SIS100RING_cold.seq","sis100ring")
lattice.setUseRealCharge(useCharge = 1)
[node.setnParts(10) for node in lattice.getNodes() if node.getType() in ["drift teapot","quad teapot"]]	

#---------------------------------------------Make a Teapot WARM Lattice--------------

print "Generate Lattice with warm quads."
latt = TEAPOT_Lattice("lattice")
latt.readMADX("sis100.madx","sis100ring")
latt.setUseRealCharge(useCharge = 1)
[node.setnParts(10) for node in latt.getNodes() if node.getType() in ["drift teapot","quad teapot"]]	

#---------------------------------------------------- ERRORS-------------------------------------------

setDict = {}
paramsDict = {}
positioni = 0.0
positionf = lattice.getLength()

setDict["elementtype"] = "quad"
setDict["ringline"] = "ring"

paramsDict["errtype"]  = "FieldError"
paramsDict["subtype"]  = "QuadField"
paramsDict["fracerr"]      = 1e-3

paramsDict["sample"]      = "Gaussian"
paramsDict["mean"]        = 0.0
paramsDict["sigma"]       = 1.0

#paramsDict["sample"]      = "Uniform"
#paramsDict["minimum"]     = -1.0
#paramsDict["maximum"]     =  1.0

ESet  = AddErrorSet(latt, positioni, positionf, setDict, paramsDict, seed_value=40)
#ESet  = AddErrorSet(latt, positioni, positionf, setDict, paramsDict) # Random
print "FIELD ERRORS IN THE QUADRUPOLES/MULTIPOLES ARE APPLIED"

#------------------------------------------------------------------------------------------------------

matrix_latt = TEAPOT_MATRIX_Lattice(latt,b)

TwissDataX,TwissDataY = matrix_latt.getRingTwissDataX(),matrix_latt.getRingTwissDataY()

beta_x = np.transpose(TwissDataX[-1])
beta_y = np.transpose(TwissDataY[-1])

#------------------------------------------------------------------------------------------------------

matrix_lattice = TEAPOT_MATRIX_Lattice(lattice,b)

TwissDataX0,TwissDataY0 = matrix_lattice.getRingTwissDataX(),matrix_lattice.getRingTwissDataY()

beta_x0 = np.transpose(TwissDataX0[-1])
beta_y0 = np.transpose(TwissDataY0[-1])

#------------------------------------------------------------------------------------------------------

sNew = np.linspace(0,1083.6,2000)

betx_cold = np.interp(x=sNew, xp=beta_x0[0],fp=beta_x0[1],period=1083.6/6.0)
bety_cold = np.interp(x=sNew, xp=beta_y0[0],fp=beta_y0[1],period=1083.6/6.0)

betx_warm = np.interp(x=sNew, xp=beta_x[0],fp=beta_x[1])
bety_warm = np.interp(x=sNew, xp=beta_y[0],fp=beta_y[1])

beatx_warm = (betx_warm/betx_cold-1)*100
beaty_warm = (bety_warm/bety_cold-1)*100

beatx_fft_warm =np.abs(np.fft.rfft(beatx_warm))/len(beatx_warm)
beaty_fft_warm =np.abs(np.fft.rfft(beaty_warm))/len(beaty_warm)

#------------------------------------------------------------------------------------------------------
'''
def save_dict(data,filename):
	w= open(filename,"w")
	w.write(",".join(data.keys())+"\n")
	name=data.keys()[0]
	n = len(data[name])
	for i in range(n):
		out = []
		for key,arr in data.items():
			out.append("{}".format(arr[i]))
		w.write(",".join(out)+"\n")
	w.close()

d = {"s":beta_x[0],"betx":beta_x[1],"bety":beta_y[1]}
save_dict(d,"quadQold.csv") 
'''
#------------------------------------------------------------------------------------------------------

bc = betaCorrection(latt,b)
bc.fixTunes(rhobeg=5e-4,A=0.5,qx0=18.84,qy0=18.73,betaX0=betx_cold,betaY0=bety_cold,patternQuad="warm",patternCorr="QS1J".lower())

#------------------------------------------------------------------------------------------------------

matrix_latt = TEAPOT_MATRIX_Lattice(latt,b)
TwissDataX,TwissDataY = matrix_latt.getRingTwissDataX(),matrix_latt.getRingTwissDataY()

beta_x = np.transpose(TwissDataX[-1])
beta_y = np.transpose(TwissDataY[-1])

betx_corr = np.interp(x=sNew, xp=beta_x[0],fp=beta_x[1])
bety_corr = np.interp(x=sNew, xp=beta_y[0],fp=beta_y[1])

beatx_corr = (betx_corr/betx_cold-1)*100
beaty_corr = (bety_corr/bety_cold-1)*100

beatx_fft_corr =np.abs(np.fft.rfft(beatx_corr))/len(beatx_corr)
beaty_fft_corr =np.abs(np.fft.rfft(beaty_corr))/len(beaty_corr)

#------------------------------------------------------------------------------------------------------

plt.figure()
plt.plot(sNew,betx_cold)
plt.plot(sNew,bety_cold)

plt.plot(sNew,betx_warm, ls="--")
plt.plot(sNew,bety_warm, ls="--")

#------------------------------------------------------------------------------------------------------

plt.figure()
plt.plot(sNew,beatx_warm)
plt.plot(sNew,beaty_warm)

plt.plot(sNew,beatx_corr, ls="--")
plt.plot(sNew,beaty_corr, ls="--")

#------------------------------------------------------------------------------------------------------

plt.figure()
plt.plot(beatx_fft_warm)
plt.plot(beaty_fft_warm)

plt.plot(beatx_fft_corr, ls="--",label="corrX")
plt.plot(beaty_fft_corr, ls="--",label="corrY")

[plt.axvline(x, ls=":") for x in [i for i in range(18*6) if not i%6]]

plt.legend()

#------------------------------------------------------------------------------------------------------

fig, ax1 = plt.subplots()

ax1.plot(beta_x[0],beta_x[1])
ax1.plot(beta_y[0],beta_y[1])

ax2 = ax1.twiny()

ax2.plot(beta_x[0],beta_x[1], ls="--")
ax2.plot(beta_y[0],beta_y[1], ls="--")

ax1.set_xlim([0,1083.6/6.0])
ax2.set_xlim([1083.6/6.0,1083.6/3.0])

#--------------------------------------------------------------------------------
plt.show()
#--------------------------------------------------------------------------------


