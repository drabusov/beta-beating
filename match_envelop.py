import os
import string
import sys
from numpy import *

from scipy.constants import c
from matplotlib.pyplot import *



from orbit.teapot import teapot
from orbit.teapot import TEAPOT_Lattice
from orbit.teapot import TEAPOT_MATRIX_Lattice

from bunch import Bunch
from bunch import BunchTwissAnalysis

from orbit.matching import Optics, EnvelopeSolver
from orbit.errors import AddErrorNode, AddErrorSet

from orbit.space_charge.sc2dslicebyslice import scAccNodes, scLatticeModifications
from spacecharge import SpaceChargeCalcSliceBySlice2D
from orbit.bunch_generators import TwissContainer, GaussDist3D, KVDist3D

rc('font', size=16)
rc('axes', linewidth=2)
rc('lines', linewidth=3)
rc('legend', fontsize=18)
rc('figure', figsize=(9,8))

#------------------------------------------------------
# Read MADX twiss file, match and plot envelopes
#-------------------------------------------------------

bunch = Bunch()

#energy = 11.4e-3  # beta must be on def. of disperion
energy = 0.001
bunch.getSyncParticle().kinEnergy(energy)

lattice = TEAPOT_Lattice("no_sc_lattice")
lattice.readMAD("fodo.lat","CELLA") # lattice start at injection point
lattice.setUseRealCharge(useCharge = 1)


matrix_lattice = TEAPOT_MATRIX_Lattice(lattice,bunch)
(muX, arrPosAlphaX, arrPosBetaX0) = matrix_lattice.getRingTwissDataX()
#(muY, arrPosAlphaY, arrPosBetaY0) = matrix_lattice.getRingTwissDataY()


#---------------------------------------------FIELD ERRORS-------------------------------------------
setDict = {}
paramsDict = {}
positioni = 0.0
positionf = lattice.getLength()

paramsDict["errtype"]  = "FieldError"
paramsDict["subtype"]  = "QuadField"
paramsDict["sample"]      = "Gauss"
paramsDict["fracerr"]      = 0.2
paramsDict["maximum"]        = 1.5
paramsDict["minimum"]       = -1.5

setDict["elementtype"] = "quad"
setDict["ringline"] = "ring"

#ESet  = AddErrorSet(lattice, positioni, positionf, setDict, paramsDict, seed_value=50)
ESet  = AddErrorSet(lattice, positioni, positionf, setDict, paramsDict)

matrix_lattice = TEAPOT_MATRIX_Lattice(lattice,bunch)
(muX, arrPosAlphaX, arrPosBetaX) = matrix_lattice.getRingTwissDataX()
(muY, arrPosAlphaY, arrPosBetaY) = matrix_lattice.getRingTwissDataY()

aux0 = np.transpose(arrPosBetaX0)
aux1 = np.transpose(arrPosBetaX)

figure()
plot(aux0[0],aux0[1])
plot(aux1[0],aux1[1])
#---------------------------------------------ORBIT ERRORS-------------------------------------------

beamline = Optics().readtwiss_teapot(lattice, bunch)

beamline.print_line()



emitx = 12.5e-6     # emittance_x
emity = 12.5e-6     # emittance_y
sigma_p = 1.0e-3    # rms momentum spread
#Ksc = 1.e-6*5       # space charge perveance
circum = 216.0      # 1080.0 # circumference

#=============================================================
clight=2.988e8
eps0=8.8542e-12
mu0=1.0/(eps0*clight**2)
qe=1.6022e-19
mp=1.6605e-27
me=9.1094e-31


N=5.0e11     # total number of ions in the ring
Bf=1     # Bunching factor (Bf=1 for coasting beam)
harmonic_rf=1     # number rf buckets in the ring
Nbunch =1

A=1        # mass number
Z=1   # charge number
e_kin=energy     # energy (MeV/u)




gamma0=1.0+(e_kin*1e6*qe)/(mp*clight*clight)
beta0=sqrt((gamma0*gamma0-1.0)/(gamma0*gamma0))

current0=beta0*clight*N*Z*qe/circum # dc current
current=(harmonic_rf/Nbunch)/Bf*beta0*clight*N*Z*qe/circum  # peak current
R=circum/(2.0*np.pi)
Ksc=2.0*current/((4.0*np.pi*mp*A*eps0*clight**3)/(qe*Z)*(beta0*gamma0)**3)
print("Perveance: {}".format(Ksc))

lc = circum/2.0/pi
dqx = Ksc*lc/(4.0*emitx*(1.0+sqrt(emity/emitx)))/pi/2.0
dqy = Ksc*lc/(4.0*emity*(1.0+sqrt(emitx/emity)))/pi/2.0

print("dQx: {}".format(dqx))
print("dQy: {}".format(dqy))
#=============================================================



solve = EnvelopeSolver(beamline)


envx0,envxs0,envy0,envys0,Dx0,Dxs0,s = solve.match_root(emitx,emity,0.0,0.0)    # matched envelopes (no space charge)
phasex0, phasey0 = solve.phase_advance(envx0,envy0,Dx0,emitx,emity,0.0,s)       # phase advance (no space charge)
print "zero current phase advance x [deg]:", phasex0*180.0/pi
print "zero current phase advance y [deg]:", phasey0*180.0/pi
envx,envxs,envy,envys,Dx,Dxs,s = solve.match_root(emitx,emity,sigma_p,Ksc)      # matched envelopes (with space charge)
phasex, phasey = solve.phase_advance(envx,envy,Dx,emitx,emity,sigma_p,s)        # phase advance (with space charge)

print "finite current phase advance x [deg]:", phasex*180.0/pi
print "finite current phase advance y [deg]:", phasey*180.0/pi


print "tuneshift x :", phasex/(2.0*pi)-phasex0/(2.0*pi)
print "tuneshift y :", phasex/(2.0*pi)-phasey0/(2.0*pi)

dqx,dqy = [x/2.0/pi for x in solve.phase_analytic(emitx,emity,Ksc,lc)]
print("dQx: {}".format(dqx))
print("dQy: {}".format(dqy))

#==================================================================================================
twissX = TwissContainer(alpha = arrPosAlphaX[0][1], beta = arrPosBetaX[0][1], emittance = emitx)
twissY = TwissContainer(alpha = arrPosAlphaY[0][1], beta = arrPosBetaY[0][1], emittance = emity)
twissZ = TwissContainer(alpha = 0., beta = 100000., emittance = sigma_p)


dist = GaussDist3D(twissX,twissY,twissZ)

n_particles = 10000
for i in range(n_particles):
    x,xp,y,yp,z,zp = dist.getCoordinates()
    bunch.addParticle(x, xp, y, yp, z, zp)

x0,px0 = [],[]
for i in range(n_particles):
    x0.append(bunch.x(i))
    px0.append(bunch.px(i))


#b.readBunch(distribution_file, n_particles)
print "Bunch Generated."

nParticlesGlobal = bunch.getSizeGlobal()
print("nParticlesGlobal = {}".format(nParticlesGlobal))
total_macroSize=1.e+14
bunch.macroSize(total_macroSize/nParticlesGlobal)


turn = 8

for i in range(turn):
    lattice.trackBunch(bunch)


x,px = [],[]
for i in range(n_particles):
    x.append(bunch.x(i))
    px.append(bunch.px(i))

#-------------------------------------------------------------------------------------------------------------
#SPACE CHARGE

#nMacrosMin = 1
sc_path_length_min = 0.00000001
sizeX = 8   #number of grid points in horizontal direction
sizeY = 8   #number of grid points in vertical direction
sizeZ = 2   #number of longitudinal slices

sc_calc = SpaceChargeCalcSliceBySlice2D(sizeX,sizeY,sizeZ)
#sc_calc = SpaceChargeCalc2p5D(sizeX,sizeY,sizeZ) # check the difference between slice by slice and 2p5D
scLatticeModifications.setSC2DSliceBySliceAccNodes(lattice, sc_path_length_min, sc_calc)

#SPACE CHARGE
#-------------------------------------------------------------------------------------------------------------


for i in range(turn):
    lattice.trackBunch(bunch)


x_sc,px_sc = [],[]
for i in range(n_particles):
    x_sc.append(bunch.x(i))
    px_sc.append(bunch.px(i))



#==================================================================================================

figure()

#scatter(x,px)
scatter(x_sc,px_sc)
scatter(x0,px0)

figure()
matplotlib.pyplot.subplot(311)
scatter(aux1[0],aux1[1])
plot(s,envx0**2/emitx,color='k',ls=':',label=r'$K=0$') # plot betax function: without space charge
plot(s,envx**2/emitx,color='k',label=r'$K>0$') # plot betax function: with space charge
xlabel(r'$s$ [m]')
ylabel(r'$\beta_x$ [m]')
legend(loc=0)

subplot(312)
#scatter(aux1[0],aux1[1])
plot(s,envy0**2/emity,color='b',ls=':',label=r'$K=0$')  # plot betay function: without space charge
plot(s,envy**2/emity,color='b',label=r'$K>0$')  # plot betay function: with space charge
xlabel(r'$s$ [m]')
ylabel(r'$\beta_y$ [m]')
legend(loc=0)


subplot(313)
plot(s,Dx0,color='k',ls='-',label=r'$\Delta Q_{x}^{sc}=0$')  # plot Dx function: without space charge
plot(s,Dx,color='r',label=r'$\Delta Q_{x}^{sc}=-0.6$')  # plot Dx function: with space charge
xlabel(r'$s$ [m]')
ylabel(r'$D_x$ [m]')
#legend(loc=0)
xlim(0.0,max(s))


subplots_adjust(bottom=0.15,left=0.15,hspace=0.5)
show()

