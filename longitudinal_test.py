import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as animation

from bunch import Bunch
from orbit.teapot import teapot
from orbit.teapot import TEAPOT_MATRIX_Lattice
from orbit.bunch_generators import TwissContainer, GaussDist2D

from orbit.orbit_correction import orbit, correction

from orbit.errors import AddErrorNode
from orbit.errors import AddErrorSet

#----------------------------------------------------------------------------------------------------

# get lattice from madx
teapot_latt = teapot.TEAPOT_Lattice()
teapot_latt.readMAD("sis100.lat","SIS100")

#---------------------------------------------Bunch init---------------------------------------------
b = Bunch()

b.mass(0.93827231)
b.macroSize(1.0)

energy = 1.0 #Gev
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
paramsDict["dx"]       = 0.006
paramsDict["dy"]       = 0.006

setDict["elementtype"] = "quad"
setDict["ringline"] = "ring"

ESet  = AddErrorSet(teapot_latt, positioni, positionf, setDict, paramsDict, seed_value=50)
#ESet  = AddErrorSet(lattice, positioni, positionf, setDict, paramsDict) # Random
#---------------------------------------------ORBIT ERRORS-------------------------------------------

# get lattice function from transfer matrices
matrix_lattice = TEAPOT_MATRIX_Lattice(teapot_latt,b)
(muX, arrPosAlphaX, arrPosBetaX) = matrix_lattice.getRingTwissDataX()
(muY, arrPosAlphaY, arrPosBetaY) = matrix_lattice.getRingTwissDataY()

#----------------------------------Bunch-Distribusion------------------------------------------------
# machted beam
emittance_x, emittance_y = 10**(-8),10**(-8) # should be checked and corrected
twissX = TwissContainer(arrPosAlphaX[0][1], arrPosBetaX[0][1], emittance_x)
twissY = TwissContainer(arrPosAlphaY[0][1], arrPosBetaY[0][1], emittance_y)



dE = 0*10**(-13)
n=2048 # num of particles
dist = GaussDist2D(twissX,twissY)
for i in range(n):
    x,px,y,py = dist.getCoordinates()
    b.addParticle(x,px,y,py,0,dE)

b.compress()

# the method calculates (<x>,<y>, emittance, num_particles) of the beam distribution
def getBPMinfo(beam):

    N = beam.getTotalCount()
    if not N:
        print('BEAM LOST')
        return 0,0,0,0,0

    m =[[beam.z(i), beam.dE(i)] for i in range(N)] 
    return np.transpose(m)


a = getBPMinfo(b)

fig, ax = plt.subplots()
line, = ax.plot(a[0],a[1], 'o')
ax.set_xlim(-10, 10)
ax.set_ylim(-10**(-13), 10**(-13))



def update(data):
    line.set_xdata(data[0])
    line.set_ydata(data[1])
    return line,


def data_gen():
    i=0
    while i<10:
        teapot_latt.trackBunch(b)
        print('TRACK')
        i+=1
        yield getBPMinfo(b)

ani = animation.FuncAnimation(fig, update, data_gen, interval=500)
plt.show()



















