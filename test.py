import numpy as np
from matplotlib import pyplot as plt

from orbit.teapot import teapot
from orbit.teapot import TEAPOT_MATRIX_Lattice

from bunch import Bunch, BunchTwissAnalysis
from orbit.bunch_generators import TwissContainer, GaussDist2D

from orbit.orbit_correction import orbit, correction

from orbit.errors import AddErrorNode, AddErrorSet

#----------------------------------------------------------------------------------------------------

# get lattice from madx
teapot_latt = teapot.TEAPOT_Lattice()
teapot_latt.readMAD("sis18.lat","SIS18")

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
paramsDict["dx"]       = 0.007
paramsDict["dy"]       = 0.007

setDict["elementtype"] = "quad"
setDict["ringline"] = "ring"

ESet  = AddErrorSet(teapot_latt, positioni, positionf, setDict, paramsDict, seed_value=50)
#ESet  = AddErrorSet(teapot_latt, positioni, positionf, setDict, paramsDict) # Random
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



dE = 0.000197*0
n=128 # num of particles
dist = GaussDist2D(twissX,twissY)
for i in range(n):
    x,px,y,py = dist.getCoordinates()
    b.addParticle(x,px,y,py,0,dE)

#----------------------------------------------------------------------------------------------------


nodes = teapot_latt.getNodes()
position = [teapot_latt.getNodePositionsDict()[node][0] for node in nodes if node.getType()=='monitor teapot'] # s coordinates of BPMs

bunchtwissanalysis = BunchTwissAnalysis()

#----------------------------------------------------------------------------------------------------
# the method calculates (<x>,<y>, emittance, num_particles) of the beam distribution
def getBPMinfo(beam):

    N = beam.getTotalCount()
    if not N:
        print('BEAM LOST')
        return 0,0,0,0,0

#    m =[[beam.x(i), beam.px(i), beam.y(i), beam.py(i)] for i in range(N)] 
    m =[[beam.x(i), beam.y(i), beam.z(i), beam.dE(i)] for i in range(N)] 

    xx = np.transpose(m)
    x_mean = np.mean(xx[0])
#    y_mean = np.mean(xx[2])
    y_mean = np.mean(xx[1])
    z_mean = np.mean(xx[2])
    E_mean = np.mean(xx[3])

#    e_x = np.sqrt(np.dot(xx[0],xx[0])*np.dot(xx[1],xx[1])-(np.dot(xx[0],xx[1]))**2 )
#    e_y = np.sqrt(np.dot(xx[2],xx[2])*np.dot(xx[3],xx[3])-(np.dot(xx[2],xx[3]))**2 )
        
#    return x_mean, y_mean, e_x, e_y, N
    return x_mean, y_mean, z_mean, E_mean, N



# the method calculates (<x>,<y>, emittance, num_particles) of the beam distribution
# The same method, but for large num of particles. dumbBunch -> write-read of files
def getBPMinfoFile(beam, filename):

    head = 13 
    m = list()
    N = beam.getTotalCount()
    if not N:
        print('BEAM LOST')
        return 0,0,0,0,0

    f = open(filename,'r')    
    for i,line in enumerate(f): 
        if i>head:
            m.append([float(x) for x in line.split()])
    f.close()

    xx = np.transpose(m)
    x_mean = np.mean(xx[0])
    y_mean = np.mean(xx[2])

    e_x = np.sqrt(np.dot(xx[0],xx[0])*np.dot(xx[1],xx[1])-(np.dot(xx[0],xx[1]))**2 )
    e_y = np.sqrt(np.dot(xx[2],xx[2])*np.dot(xx[3],xx[3])-(np.dot(xx[2],xx[3]))**2 )
        
    return x_mean, y_mean, e_x, e_y, N


#----------------------------------------------------------------------------------------------------
# the method calculates (<x>,<y>, emittance, num_particles) of the beam distribution for all BPMs in the lattice
# the tracking takes place inside this method
def getAllBPMsInfo_old(beam, lattice_nodes, filename):
    lst = list()
    for node in lattice_nodes:
        node.trackBunch(beam)
#        beam.dumpBunch(filename) 
        if node.getType()=='monitor teapot':
            lst.append(getBPMinfo(beam))
#            lst.append(getBPMinfoFile(beam,filename))
    return lst


def getAllBPMsInfo(beam, lattice_nodes):
    lst = list()
    for node in lattice_nodes:
        node.trackBunch(beam)
        if node.getType()=='monitor teapot':
            bunchtwissanalysis.analyzeBunch(beam)
            x = bunchtwissanalysis.getAverage(0)
            y = bunchtwissanalysis.getAverage(2)
            lst.append([x,y])
    return lst

#----------------------------------------------------------------------------------------------------

# in order to get beam oscillations, one can repeat getALLBPMsInfo method n_turn times

n_turn = 1024
#arr = [getAllBPMsInfo(b,nodes,'final.dat') for i in range(n_turn)]
arr = [getAllBPMsInfo(b,nodes) for i in range(n_turn)]

        
lst = np.transpose(arr)
x_full = lst[0]
y_full = lst[1]



x_mean = [np.mean(x) for x in x_full]
y_mean = [np.mean(y) for y in y_full]

#---------------------------------------------CALCULATE ORBIT------------------------------

OrbitX, OrbitY = orbit(teapot_latt,b).get_orbit() 

x = []
y = []
s = []
for i in xrange(len(OrbitX)):
	s.append(OrbitX[i][0])
	x.append(OrbitX[i][1])
	y.append(OrbitY[i][1])



#----------------------------------------------------------------------------------------------------


plt.figure()
plt.plot(s,x,'r-', label="x")
plt.scatter(position,x_mean)
plt.legend(loc="lower right")
plt.grid(True, ls = '--')
plt.ylim([min(x_mean),max(x_mean)])
#plt.savefig('X-orbit_quad.pdf')

plt.figure()
plt.plot(s,y,'b-', label="y")
plt.scatter(position,y_mean)
plt.legend(loc="lower right")
plt.grid(True, ls = '--')
plt.ylim([min(y_mean),max(y_mean)])
#plt.savefig('Y-orbit_quad.pdf')



'''
plt.figure()
plt.plot(ex_full[0], label="X_emittance")
plt.plot(ey_full[0], label="Y_emittance")
plt.legend(loc="lower right")
plt.grid(True, ls = '--')
#plt.savefig('emittance.pdf')
'''
plt.show()








