import math
import matplotlib.pyplot as plt
import numpy as np

from bunch import Bunch
from orbit.teapot import TEAPOT_Lattice
from orbit.teapot import TEAPOT_MATRIX_Lattice
from orbit.errors import AddErrorNode,AddErrorSet

from orbit.beta_correction import betaCorrection

#-------------------------------------------------------------------------

def calculate_betas(teapot_lattice, bunch):

	matrix_lattice = TEAPOT_MATRIX_Lattice(teapot_lattice,bunch)
	TwissDataX,TwissDataY = matrix_lattice.getRingTwissDataX(),matrix_lattice.getRingTwissDataY()
	
	return np.transpose(TwissDataX[-1]), np.transpose(TwissDataY[-1])

def calculate_beat(s_inter, betx_ini,bety_ini,betx_fin,bety_fin):

	betx_inter = np.interp(x=s_inter, xp=betx_fin[0],fp=betx_fin[1])
	bety_inter = np.interp(x=s_inter, xp=bety_fin[0],fp=bety_fin[1])

	beatx = betx_inter-betx_ini
	beaty = bety_inter-bety_ini

	#beatx_fft =np.abs(np.fft.rfft(beatx))/len(beatx)
	#beaty_fft =np.abs(np.fft.rfft(beaty))/len(beaty)

	return beatx,beaty

def save_dict(data,filename,qx,qy):
	w= open(filename,"w")
	w.write("Qx,Qy = {},{}\n".format(qx,qy))
	w.write(",".join(data.keys())+"\n")
	name=data.keys()[0]
	n = len(data[name])
	for i in range(n):
		out = []
		for key,arr in data.items():
			out.append("{}".format(arr[i]))
		w.write(",".join(out)+"\n")
	w.close()

#-------------------------------------------------------------------------
#useful globals

tmp =["{}".format(y) for y in range(1,10)]+["a","b","c","d","e"]

focusing = ["s{}{}qd12".format(x,y) for x in range(1,7) for y in tmp]
defocusing = ["s{}{}qd11".format(x,y) for x in range(1,7) for y in tmp]

focusingWarm = [x for x in focusing if x != "s52qd12"]
defocusingWarm = [x for x in defocusing if x != "s52qd11"]

qDict = {"corrector":["S52QD12","S52QD11"],"group1":focusing,"group2":defocusing}
qDictWarm = {"corrector":["S52QD12","S52QD11"],"group1":focusingWarm,"group2":defocusingWarm}

#---------------------------------------------Bunch init---------------------------------------------
print "Start."

b = Bunch()
energy = 1
syncPart=b.getSyncParticle()
syncPart.kinEnergy(energy)

#---------------------------------------------Make a Teapot COLD Lattice ---------------

print "Generate ideal cold Lattice."
lattice = TEAPOT_Lattice("lattice")
lattice.readMADX("sis100cold.seq","sis100ring")
lattice.setUseRealCharge(useCharge = 1)
[node.setnParts(10) for node in lattice.getNodes() if node.getType() in ["drift teapot","quad teapot"]]	

#---------------------------------------------Make a Teapot WARM Lattice--------------

print "Generate Lattice with warm quads."
latt = TEAPOT_Lattice("lattice")
latt.readMADX("sis100.madx","sis100ring")
latt.setUseRealCharge(useCharge = 1)
[node.setnParts(10) for node in latt.getNodes() if node.getType() in ["drift teapot","quad teapot"]]	

#-------------------------------------- MAIN PART -------------------------------------------------

arrX = np.linspace(18.55,18.95,50)
arrY = np.linspace(18.55,18.95,50)

#tunes = [(x,y) for x in arr for y in arr]
tunes = [(x,y) for x in arrX for y in arrY]

w = open("summary.csv","w")
w.write("qx,qy,kw1,kw2,kd,kf,kw1_ini,kw2_ini,kd_ini,kf_ini,max0,max,aper0,aper\n")

for idx,wp in enumerate(tunes):
	Qx,Qy = wp
	print("tunes are {} {}".format(Qx,Qy))

	#---------------- match cold lattice
	bc = betaCorrection(lattice,b)
	bc.matchTunes(qx0=Qx,qy0=Qy,quadsDict=qDict,A=0.01)

	#---------------- match warm lattice (kw_i = a_i*k_i)
	bc = betaCorrection(latt,b)
	bc.matchTunes(qx0=Qx,qy0=Qy,quadsDict=qDictWarm,A=0.01,warm=True)

	# to summary
	kw1_ini = latt.getNodeForName("s52qd11").getParam("kq")
	kw2_ini = latt.getNodeForName("s52qd12").getParam("kq")
	kd_ini = latt.getNodeForName("s11qd11").getParam("kq")
	kf_ini = latt.getNodeForName("s11qd12").getParam("kq")
	#----------------

	beta_x0,beta_y0 = calculate_betas(lattice, b)
	beta_x,beta_y = calculate_betas(latt, b)
	
	sNew = np.linspace(0,1083.6,2000)
	betx_cold = np.interp(x=sNew, xp=beta_x0[0],fp=beta_x0[1],period=1083.6/6.0)
	bety_cold = np.interp(x=sNew, xp=beta_y0[0],fp=beta_y0[1],period=1083.6/6.0)
	
	beatx_warm,beaty_warm =calculate_beat(sNew, betx_cold,bety_cold,beta_x,beta_y)
	
	#---------------- rock'n'roll
	
	bcNew = betaCorrection(latt,b)
	f0,f,a0,a = bcNew.correction(qx0=Qx,qy0=Qy,betaX0=betx_cold,betaY0=bety_cold,quadsDict=qDictWarm)
	
	#------------------------------------------------------------------------------------------------------
	
	# to summary
	kw1 = latt.getNodeForName("s52qd11").getParam("kq")
	kw2 = latt.getNodeForName("s52qd12").getParam("kq")
	kd = latt.getNodeForName("s11qd11").getParam("kq")
	kf = latt.getNodeForName("s11qd12").getParam("kq")


	w.write("{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n".format(Qx,Qy,kw1,kw2,kd,kf,kw1_ini,kw2_ini,kd_ini,kf_ini,f0,f,a0,a))

	# to ouput files
	betx,bety = calculate_betas(latt, b)
	beatx_corr,beaty_corr =calculate_beat(sNew, betx_cold,bety_cold,betx,bety)
	
	#d = {"s":sNew,"beatx_ini":beatx_warm, "beaty_ini":beatx_warm, "beatx_corr":beatx_corr, "beaty_corr":beaty_corr}
	#save_dict(d,"./data/beat_output_{}.csv".format(idx), Qx,Qy) 

w.close()

print("rock'n'roll completed successfully")
#---------------------------------------------------- ERRORS-------------------------------------------
'''
fracerr = 1e-3
positioni = 0.0
positionf = latt.getLength()

setDict = {"elementtype":"quad","ringline":"ring"}
paramsDict = {"errtype":"FieldError","subtype":"QuadField","fracerr":fracerr,"sample":"Gaussian","mean":0.0,"sigma":1.0}


#paramsDict["sample"]      = "Uniform"
#paramsDict["minimum"]     = -1.0
#paramsDict["maximum"]     =  1.0

#ESet  = AddErrorSet(latt, positioni, positionf, setDict, paramsDict, seed_value=40)
#ESet  = AddErrorSet(latt, positioni, positionf, setDict, paramsDict) # Random
#print "FIELD ERRORS IN THE QUADRUPOLES/MULTIPOLES ARE APPLIED"
'''

