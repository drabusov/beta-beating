import time
import numpy as np
from matplotlib import rc
from matplotlib import pyplot as plt

from orbit.teapot import teapot
from orbit.teapot import TEAPOT_MATRIX_Lattice

from bunch import Bunch, BunchTwissAnalysis
from orbit.bunch_generators import TwissContainer, GaussDist2D, KVDist2D, KVDist3D

from orbit.orbit_correction import orbit, correction
from orbit.errors import AddErrorNode, AddErrorSet

from spacecharge import Boundary2D
from orbit.space_charge.sc2dslicebyslice import scAccNodes, scLatticeModifications
from spacecharge import SpaceChargeCalcSliceBySlice2D

#----------------------------------------------------------------------------------------------------

rc('text', usetex=True)
rc('font', family='serif')
rc('font', size=20)
rc('axes', linewidth=1.5)
rc('lines', linewidth=1.5)
rc('legend', fontsize=18)
rc('xtick', labelsize='small')
rc('ytick', labelsize='small')
rc('figure', figsize=(8,6))

#----------------------------------------------------------------------------------------------------


begin = time.time()
#----------------------------------------------------------------------------------------------------

# get lattice from madx
teapot_latt = teapot.TEAPOT_Lattice()
teapot_latt.readMAD("sis100.lat","SIS100")

#---------------------------------------------Bunch init---------------------------------------------
b = Bunch()

b.mass(0.93827231)
b.macroSize(1.0)

energy = 1 #Gev
b.getSyncParticle().kinEnergy(energy)
#----------------------------------------------------------------------------------------------------

#---------------------------------------------ORBIT ERRORS-------------------------------------------
# WE INTRODUCE MISALIGNEMENT IN THE QUADRUPOLES; dx, dy = HOR AND VER DISPLACEMENT OF QUADRUPOLES
setDict = {}
paramsDict = {}
positioni = 0.0
positionf = teapot_latt.getLength()
paramsDict["errtype"]  = "StraightError"
paramsDict["subtype"]  = "TransDisp"
paramsDict["sample"]      = "Uniform"
paramsDict["maximum"]        = 0.5
paramsDict["minimum"]       = 0.0
paramsDict["dx"]       = 0.001*0
paramsDict["dy"]       = 0.001*0

setDict["elementtype"] = "quad"
setDict["ringline"] = "ring"

#ESet  = AddErrorSet(teapot_latt, positioni, positionf, setDict, paramsDict, seed_value=50)

paramsDict["errtype"]  = "FieldError"
paramsDict["subtype"]  = "QuadField"
paramsDict["sample"]      = "Uniform"
paramsDict["fracerr"]      = 0.005

setDict["elementtype"] = "quad"
setDict["ringline"] = "ring"

ESet  = AddErrorSet(teapot_latt, positioni, positionf, setDict, paramsDict, seed_value=2)


#----------------------------------------------------------------------------------------------------
'''
# not matched	 3D
twissX = TwissContainer(alpha = 0.0046902, beta = 10.207, emittance = 3.e-5)
twissY = TwissContainer(alpha = 0.056823, beta = 10.639, emittance = 3.e-5)
twissZ = TwissContainer(alpha = 0., beta = 100000., emittance = 0.008)


dist = KVDist3D(twissX,twissY,twissZ)

n = 1000
for i in range(n):
    x,xp,y,yp,z,zp = dist.getCoordinates()
    b.addParticle(x, xp, y, yp, z, zp)




#b.readBunch(distribution_file, n_particles)
print "Bunch Generated."

nParticlesGlobal = b.getSizeGlobal()
total_macroSize=1.e+14
b.macroSize(total_macroSize/nParticlesGlobal)
# 3D
'''
#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
'''
#SPACE CHARGE
nMacrosMin = 1
sc_path_length_min = 0.00000001
sizeX = 8   #number of grid points in horizontal direction
sizeY = 8   #number of grid points in vertical direction
sizeZ = 2   #number of longitudinal slices

# Make a boundary
bpoints = 128
bmodes = 10
xdim = 0.073
ydim = 0.073
boundary = Boundary2D(bpoints, bmodes, "Circle", xdim, ydim)

sc_calc = SpaceChargeCalcSliceBySlice2D(sizeX,sizeY,sizeZ)
scLatticeModifications.setSC2DSliceBySliceAccNodes(teapot_latt, sc_path_length_min, sc_calc)

#SPACE CHARGE
'''
#-------------------------------------------------------------------------------------------------------------
#---------------------------------------------ORBIT ERRORS-------------------------------------------

# get lattice function from transfer matrices
matrix_lattice = TEAPOT_MATRIX_Lattice(teapot_latt,b)
(muX, arrPosAlphaX, arrPosBetaX) = matrix_lattice.getRingTwissDataX()
(muY, arrPosAlphaY, arrPosBetaY) = matrix_lattice.getRingTwissDataY()
#----------------------------------Bunch-Distribusion------------------------------------------------
# machted beam

emittance_x, emittance_y = 35*10**(-6),15*10**(-6) 
twissX = TwissContainer(arrPosAlphaX[0][1], arrPosBetaX[0][1], emittance_x)
twissY = TwissContainer(arrPosAlphaY[0][1], arrPosBetaY[0][1], emittance_y)

#twissZ = TwissContainer(alpha = 0., beta = 1., emittance = 10**(-6)) # not matched

#---------------------------------------------CALCULATE ORBIT------------------------------
OrbitX, OrbitY = orbit(teapot_latt,b).get_orbit() 

x0, px0 = OrbitX[0][1], OrbitX[0][2]
y0, py0 = OrbitY[0][1], OrbitY[0][2]
print('Closed x-orbit at s_0 ({},{})'.format(x0,px0))
#---------------------------------------------------------------------------------------------------
#x0 = 0.001

n=10**(2) # num of particles
#dist = KVDist3D(twissX,twissY,twissZ)
dist = KVDist2D(twissX,twissY)
for i in range(n):
    xi,pxi,yi,pyi = dist.getCoordinates()
    b.addParticle(x0+xi,px0+pxi,y0+yi,py0+pyi,0,0)
#    xi,pxi,yi,pyi,z,zp = dist.getCoordinates()
#    b.addParticle(x0+xi,px0+pxi,y0+yi,py0+pyi,z,zp)

#----------------------------------------------------------------------------------------------------

nodes = teapot_latt.getNodes()
position = [teapot_latt.getNodePositionsDict()[node][0] for node in nodes if node.getType()=='monitor teapot'] # s coordinates of BPMs


#----------------------------------------------------------------------------------------------------
# the method calculates (<x>,<y>) of the beam distribution for all BPMs in the lattice
# the tracking takes place inside this method


def getAllBPMsInfo(beam, lattice_nodes):
    lst = list()
    for node in lattice_nodes:
        node.trackBunch(beam)
        if node.getType()=='monitor teapot':
            x = node.getParam("xAvg")
            y = node.getParam("yAvg")
            lst.append([x,y])
            
    return lst


# there is no difference in run time
'''
#try to change to generator style
def getAllBPMsInfo(beam, lattice_nodes):
    return [[node.getParam("xAvg"),node.getParam("yAvg")] for node in lattice_nodes if not node.trackBunch(beam) and node.getType()=='monitor teapot']
'''

#----------------------------------------------------------------------------------------------------

# in order to get beam oscillations (beamhistory), one can repeat getALLBPMsInfo method n_turn times



n_turn = 1000
start = time.time()
print('TRACKING. PREPARATION TIME = {}'.format(start - begin))
arr = [getAllBPMsInfo(b,nodes) for i in range(n_turn)]
stop = time.time()
print('TRACKING IS FINISHED. DURATION = {}. NUMBER OF PARTICLES = {}. NUMBER OF TURNS = {}'.format(stop-start, n, n_turn))
        



# write all bpm data to a file
f = open('bpm.dat','w')
for a in arr:
    s =' '.join([str(x[0]) for x in a])
    f.write(s+'\n')

f.close()


lst = np.transpose(arr)
x_full = lst[0]
y_full = lst[1]



x_mean = [np.mean(elem) for elem in x_full]
y_mean = [np.mean(elem) for elem in y_full]

mx = [elem-np.mean(elem) for elem in x_full]/np.sqrt(n_turn)
my = [elem-np.mean(elem) for elem in y_full]/np.sqrt(n_turn)

U, s, V = np.linalg.svd(np.transpose(mx), full_matrices=True)
betaX = pow(s[0],2)*(V[0]**2)+pow(s[1],2)*(V[1]**2)
print('singular values are: {}'.format(s))

U, s, V = np.linalg.svd(np.transpose(my), full_matrices=True)
betaY = pow(s[0],2)*(V[0]**2)+pow(s[1],2)*(V[1]**2)

beta_orbit_x = np.transpose(arrPosBetaX)
beta_orbit_y = np.transpose(arrPosBetaY)




beta_bpm = list()
for c in position:
    i=0
    while np.abs(beta_orbit_x[0][i] - c) > 10**(-6):
        i+=1
    beta_bpm.append(beta_orbit_x[1][i])


kappaX = beta_bpm[0]/betaX[0]
print('scaling factor for X = {} \n scaling facor for Y is {}'.format(kappaX, 'unknown'))



# ORBIT
x = []
y = []
s = []
for i in xrange(len(OrbitX)):
	s.append(OrbitX[i][0])
	x.append(OrbitX[i][1])
	y.append(OrbitY[i][1])




def count_metric(arr):
    arr_period = np.transpose(np.reshape(arr, (6,14)))
    return [np.std(b)/np.mean(b) for b in arr_period] 

metric_mia = count_metric(betaX)
metric_pyorbit = count_metric(beta_bpm) 

plt.figure()
plt.plot(metric_mia, label = 'MIA')
plt.plot(metric_pyorbit, label = 'pyORBIT')
plt.legend(loc="right")

#----------------------------------------------------------------------------------------------------
#Comparison between pyORBIT and multi-turn tracking methods


plt.figure()
#plt.plot(beta_orbit_x[0],beta_orbit_x[1], label="X-beta pyORBIT")

plt.scatter(position,beta_bpm, color = 'black', label=r"$\beta_x$ by pyORBIT")
plt.scatter(position,kappaX*betaX, color = 'red', marker ='x', label=r"$\beta_x$  by MIA, {} turns, {} particles". format(n_turn,n))
plt.plot([], [], ' ', label='Gradient error $\Delta G = ${}'.format(paramsDict["fracerr"]))
plt.axvline(s[-1]/6.)

plt.legend(loc="right")
plt.grid(True, ls = '--')
#plt.xlim([0,s[-1]/6])
plt.title(r'SIS100 full lattice. $\beta (s)$ at BPMs')
plt.xlabel('s [m]')
plt.ylabel(r'$\beta(s)$ [m]')
plt.savefig('betaX.pdf')




'''
plt.figure()
plt.plot(s,x,'r-x', label="x orbit, standart pyORBIT")
plt.scatter(position,x_mean, label="x at BPMs, multi-turn tracking method")
plt.legend(loc="lower right")
plt.grid(True, ls = '--')
#plt.ylim([min(x),max(x)])
plt.title('X-plane. {} particles, {} turns'.format(n,n_turn))
#plt.title('X-plane. Drift length = 4 instead of 3.062')
#plt.savefig('X-orbit_drift.pdf')

plt.figure()
plt.plot(s,y,'g-x', label="y, standart pyORBIT")
plt.scatter(position,y_mean,label="x at BPMs, multi-turn tracking method")
plt.legend(loc="lower right")
plt.grid(True, ls = '--')
#plt.ylim([min(y),max(y)])
plt.title('Y-plane. {} particles, {} turns'.format(n,n_turn))
#plt.title('Y-plane. Drift length = 4 instead of 3.062')
#plt.savefig('Y-orbit_drift.pdf')
'''

plt.show()









