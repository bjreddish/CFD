import numpy as np 
import matplotlib.pyplot as plt
import h5py


def plotFlow(pres,temp,uVel,vVel):
	"""
	Plotting script to show flowfiled 
	"""
	# Pressure
	plt.figure(1)
	# plt.subplot(2,2,1)
	plt.contourf(pres.transpose())
	plt.colorbar(label='N/m^2')
	plt.title('Pressure')
	# Temperature
	plt.figure(2)
	# plt.subplot(2,2,2)
	plt.contourf(temp.transpose())
	plt.colorbar(label='K')
	plt.title('Temp')
	# U Vel
	plt.figure(3)
	# plt.subplot(2,2,3)
	plt.contourf(uVel.transpose())
	plt.colorbar(label='m/s')
	plt.title('U Vel')
	# V Vel
	plt.figure(4)
	# plt.subplot(2,2,4)
	plt.contourf(vVel.transpose())
	plt.colorbar(label='m/s')
	plt.title('V Vel')
	plt.show()
	return
filename = 'baseCase.h5'
f = h5py.File(filename, 'r')
# List all groups
presKey = list(f.keys())[0]
tempKey = list(f.keys())[2]
uKey = list(f.keys())[3]
vKey = list(f.keys())[4]
residualKey = list(f.keys())[1]
print(list(f.keys())[:])
# Get the data
presDomain = f[presKey][:,:]
tempDomain = f[tempKey][:,:]
uVelDomain = f[uKey][:,:]
vVelDomain = f[vKey][:,:]
residualKey = f[residualKey][:]
plotFlow(presDomain,tempDomain,uVelDomain,vVelDomain)
plt.plot(residual)
plt.show()
