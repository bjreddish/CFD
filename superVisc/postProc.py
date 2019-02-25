import numpy as np 
import matplotlib.pyplot as plt
import h5py


def plotFlow(x,y,pres,temp,uVel,vVel):
	"""
	Plotting script to show flowfiled 
	"""
	# Create mesh for points
	X, Y = np.meshgrid(x, y)
	# Pressure
	plt.figure(1)
	# plt.subplot(2,2,1)
	plt.contourf(X,Y,pres.transpose())
	plt.colorbar(label='N/m^2')
	plt.xlabel('x(m)')
	plt.ylabel('y(m)')
	plt.grid(True)
	plt.title('Pressure')
	# Temperature
	plt.figure(2)
	# plt.subplot(2,2,2)
	plt.contourf(X,Y,temp.transpose())
	plt.colorbar(label='K')
	plt.xlabel('x(m)')
	plt.ylabel('y(m)')
	plt.grid(True)
	plt.title('Temp')
	# U Vel
	plt.figure(3)
	# plt.subplot(2,2,3)
	plt.contourf(X,Y,uVel.transpose())
	plt.colorbar(label='m/s')
	plt.xlabel('x(m)')
	plt.ylabel('y(m)')
	plt.grid(True)
	plt.title('U Vel')
	# V Vel
	plt.figure(4)
	# plt.subplot(2,2,4)
	plt.contourf(X,Y,vVel.transpose())
	plt.colorbar(label='m/s')
	plt.xlabel('x(m)')
	plt.ylabel('y(m)')
	plt.grid(True)
	plt.title('V Vel')
	plt.show()
	return

f = h5py.File(input('File Name:'), 'r')
# List all groups
# print(list(f.keys())[:])
# Get the data
presDomain = f['pres'][:,:]
tempDomain = f['temp'][:,:]
uVelDomain = f['u'][:,:]
vVelDomain = f['v'][:,:]
residual = f['residual'][:]
lengthOfPlate = f['lengthOfPlate']
deltax = f['deltax']
deltay = f['deltay']
xNodes,yNodes = presDomain.shape
x = np.arange(xNodes)*deltax
y = np.arange(yNodes)*deltay
plotFlow(x,y,presDomain,tempDomain,uVelDomain,vVelDomain)
plt.plot(residual)
plt.grid(True)
plt.yscale('log')
plt.xlabel('Iterations')
plt.ylabel('Residual')
plt.show()
