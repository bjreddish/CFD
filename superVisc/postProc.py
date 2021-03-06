import numpy as np 
import matplotlib.pyplot as plt
import h5py


def plotFlow(x,y,pres,temp,uVel,vVel,vTotal):
	"""
	Plotting script to show flowfiled 
	"""
	# Create mesh for points
	X, Y = np.meshgrid(x, y)
	levs = 50
	# Pressure
	plt.figure(1)
	# plt.subplot(2,2,1)
	plt.contourf(X,Y,pres.transpose(),levs)
	plt.colorbar(label='N/m^2')
	plt.xlabel('x(m)')
	plt.ylabel('y(m)')
	plt.grid(True)
	plt.title('Pressure')
	# Temperature
	plt.figure(2)
	# plt.subplot(2,2,2)
	plt.contourf(X,Y,temp.transpose(),levs)
	plt.colorbar(label='K')
	plt.xlabel('x(m)')
	plt.ylabel('y(m)')
	plt.grid(True)
	plt.title('Temp')
	# U Vel
	plt.figure(3)
	# plt.subplot(2,2,3)
	plt.contourf(X,Y,uVel.transpose(),levs)
	plt.colorbar(label='m/s')
	plt.xlabel('x(m)')
	plt.ylabel('y(m)')
	plt.grid(True)
	plt.title('U Vel')
	# V Vel
	plt.figure(4)
	# plt.subplot(2,2,4)
	plt.contourf(X,Y,vVel.transpose(),levs)
	plt.colorbar(label='m/s')
	plt.xlabel('x(m)')
	plt.ylabel('y(m)')
	plt.grid(True)
	plt.title('V Vel')
	plt.show()

	plt.contourf(X,Y,vTotal.transpose(),levs)
	plt.colorbar(label='m/s')
	plt.xlabel('x(m)')
	plt.ylabel('y(m)')
	plt.grid(True)
	plt.title('Velocity')
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
vTotal = (uVelDomain**2 + vVelDomain**2)**0.5
residual = f['residual'][:]
lengthOfPlate = f['lengthOfPlate']
deltax = f['deltax']
deltay = f['deltay']
xNodes,yNodes = presDomain.shape
x = np.arange(xNodes)*deltax
y = np.arange(yNodes)*deltay
plotFlow(x,y,presDomain,tempDomain,uVelDomain,vVelDomain,vTotal)
plt.plot(residual)
plt.grid(True)
plt.yscale('log')
plt.xlabel('Iterations')
plt.ylabel('Residual')
plt.show()
