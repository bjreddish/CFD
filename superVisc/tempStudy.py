import numpy as np 
import matplotlib.pyplot as plt
import h5py
def header(msg):
	print('-' * 50)
	print('['+msg+']')
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

f1 = h5py.File('adiabaticFull.h5', 'r')
presDomain1 = f1['pres'][:,:]
tempDomain1 = f1['temp'][:,:]
rhoDomain1 = presDomain1/(287*tempDomain1)
uVelDomain1 = f1['u'][:,:]
vVelDomain1 = f1['v'][:,:]
residual1 = f1['residual'][:]
lengthOfPlate1 = f1['lengthOfPlate']
deltax1 = f1['deltax']
deltay1 = f1['deltay']
xNodes1,yNodes1 = presDomain1.shape
x1 = np.arange(xNodes1)*deltax1
y1 = np.arange(yNodes1)*deltay1
print(presDomain1.shape)


f2 = h5py.File('nonad.h5', 'r')
presDomain2 = f2['pres'][:,:]
tempDomain2 = f2['temp'][:,:]
rhoDomain2 = presDomain2/(287*tempDomain2)
uVelDomain2 = f2['u'][:,:]
vVelDomain2 = f2['v'][:,:]
residual2 = f2['residual'][:]
lengthOfPlate2 = f2['lengthOfPlate']
deltax2 = f2['deltax']
deltay2 = f2['deltay']
xNodes2,yNodes2 = presDomain2.shape
x2 = np.arange(xNodes2)*deltax2
y2 = np.arange(yNodes2)*deltay2




# nondimensionalizing y
y1bar = np.linspace(0,(y1[-1]/x1[-1]) * (931.944399)**0.5,yNodes1)
y2bar = np.linspace(0,(y2[-1]/x2[-1]) * (931.944399)**0.5,yNodes2)



# Plotting at outflow (trailing edge)
# pressure
plt.plot(presDomain1[int(xNodes1)-1,:]/presDomain1[0,-1],y1bar,label='Adiabatic')
plt.plot(presDomain2[int(xNodes2)-1,:]/presDomain2[0,-1],y2bar,label='Set Wall Temp')
plt.title('Pressure at outflow')
plt.legend()
plt.grid(True)
plt.xlabel('Pressure')
plt.ylabel('Normalized y')
plt.show()


# u Vel
plt.plot(uVelDomain1[int(xNodes1)-1,:]/uVelDomain1[0,-1],y1bar,label='Adiabatic')
plt.plot(uVelDomain2[int(xNodes2)-1,:]/uVelDomain1[0,-1],y2bar,label='Set Wall Temp')
plt.title('Velocity in X at Outflow')
plt.legend(loc=2)
plt.grid(True)
plt.xlabel('u-vel')
plt.ylabel('Normalized y')
plt.show()

# Temp 
plt.plot(tempDomain1[int(xNodes1)-1,:]/tempDomain1[0,-1],y1bar,label='Adiabatic')
plt.plot(tempDomain2[int(xNodes2)-1,:]/tempDomain1[0,-1],y2bar,label='Set Wall Temp')
plt.title('Temperature at Outflow')
plt.legend()
plt.grid(True)
plt.xlabel('Temperature')
plt.ylabel('Normalized y')
plt.show()

# Plotting at plate
plt.plot(x1,presDomain1[:,0]/presDomain1[0,-1],label='Adiabatic')
plt.plot(x2,presDomain2[:,0]/presDomain2[0,-1],label='Set Wall Temp')
plt.title('Pressure at Surface')
plt.legend()
plt.grid(True)
plt.xlabel('x-location on plate')
plt.ylabel('nondimensionalized pressure')
plt.show()