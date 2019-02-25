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

f1 = h5py.File('baseCase.h5', 'r')
presDomain1 = f1['pres'][:,:]
tempDomain1 = f1['temp'][:,:]
uVelDomain1 = f1['u'][:,:]
vVelDomain1 = f1['v'][:,:]
residual1 = f1['residual'][:]
lengthOfPlate1 = f1['lengthOfPlate']
deltax1 = f1['deltax']
deltay1 = f1['deltay']
xNodes1,yNodes1 = presDomain1.shape
x1 = np.arange(xNodes1)*deltax1
y1 = np.arange(yNodes1)*deltay1


f2 = h5py.File('100base.h5', 'r')
presDomain2 = f2['pres'][:,:]
tempDomain2 = f2['temp'][:,:]
uVelDomain2 = f2['u'][:,:]
vVelDomain2 = f2['v'][:,:]
residual2 = f2['residual'][:]
lengthOfPlate2 = f2['lengthOfPlate']
deltax2 = f2['deltax']
deltay2 = f2['deltay']
xNodes2,yNodes2 = presDomain2.shape
x2 = np.arange(xNodes2)*deltax2
y2 = np.arange(yNodes2)*deltay2


f3 = h5py.File('140base.h5', 'r')
presDomain3 = f3['pres'][:,:]
tempDomain3 = f3['temp'][:,:]
uVelDomain3 = f3['u'][:,:]
vVelDomain3 = f3['v'][:,:]
residual3 = f3['residual'][:]
lengthOfPlate3 = f3['lengthOfPlate']
deltax3 = f3['deltax']
deltay3 = f3['deltay']
xNodes3,yNodes3 = presDomain3.shape
x3 = np.arange(xNodes3)*deltax3
y3 = np.arange(yNodes3)*deltay3

f4 = h5py.File('200base.h5', 'r')
presDomain4 = f4['pres'][:,:]
tempDomain4 = f4['temp'][:,:]
uVelDomain4 = f4['u'][:,:]
vVelDomain4 = f4['v'][:,:]
residual4 = f4['residual'][:]
lengthOfPlate4 = f4['lengthOfPlate']
deltax4 = f4['deltax']
deltay4 = f4['deltay']
xNodes4,yNodes4 = presDomain4.shape
x4 = np.arange(xNodes4)*deltax4
y4 = np.arange(yNodes4)*deltay4

f5 = h5py.File('280base.h5', 'r')
presDomain5 = f5['pres'][:,:]
tempDomain5 = f5['temp'][:,:]
uVelDomain5 = f5['u'][:,:]
vVelDomain5 = f5['v'][:,:]
residual5 = f5['residual'][:]
lengthOfPlate5 = f5['lengthOfPlate']
deltax5 = f5['deltax']
deltay5 = f5['deltay']
xNodes5,yNodes5 = presDomain5.shape
x5 = np.arange(xNodes5)*deltax5
y5 = np.arange(yNodes5)*deltay5


# nondimensionalizing y
y1bar = np.linspace(0,(y1[-1]/x1[-1]) * (931.944399)**0.5,yNodes1)
y2bar = np.linspace(0,(y2[-1]/x2[-1]) * (931.944399)**0.5,yNodes2)
y3bar = np.linspace(0,(y3[-1]/x3[-1]) * (931.944399)**0.5,yNodes3)
y4bar = np.linspace(0,(y4[-1]/x4[-1]) * (931.944399)**0.5,yNodes4)
y5bar = np.linspace(0,(y5[-1]/x5[-1]) * (931.944399)**0.5,yNodes5)


# Plotting at outflow (trailing edge)
plt.plot(presDomain1[int(yNodes1)-1,:]/presDomain1[0,-1],y1bar,label='4,900 Cells')
plt.plot(presDomain2[int(yNodes2)-1,:]/presDomain2[0,-1],y2bar,label='10,000 Cells')
plt.plot(presDomain3[int(yNodes3)-1,:]/presDomain3[0,-1],y3bar,label='19,600 Cells')
plt.plot(presDomain4[int(yNodes4)-1,:]/presDomain4[0,-1],y4bar,label='40,000 Cells')
plt.plot(presDomain5[int(yNodes5)-1,:]/presDomain5[0,-1],y5bar,label='80,000 Cells')
plt.title('Pressure at outflow')
plt.legend()
plt.grid(True)
plt.xlabel('Pressure')
plt.ylabel('Normalized y')
plt.show()
plt.plot(uVelDomain1[int(yNodes1)-1,:]/uVelDomain1[0,-1],y1bar,label='4,900 Cells')
plt.plot(uVelDomain2[int(yNodes2)-1,:]/uVelDomain1[0,-1],y2bar,label='10,000 Cells')
plt.plot(uVelDomain3[int(yNodes3)-1,:]/uVelDomain1[0,-1],y3bar,label='19,600 Cells')
plt.plot(uVelDomain4[int(yNodes4)-1,:]/uVelDomain1[0,-1],y4bar,label='40,000 Cells')
plt.plot(uVelDomain5[int(yNodes5)-1,:]/uVelDomain1[0,-1],y5bar,label='80,000 Cells')
plt.title('u-vel at outflow')
plt.legend()
plt.grid(True)
plt.xlabel('u-vel')
plt.ylabel('Normalized y')
plt.show()
plt.plot(tempDomain1[int(yNodes1)-1,:]/tempDomain1[0,-1],y1bar,label='4,900 Cells')
plt.plot(tempDomain2[int(yNodes2)-1,:]/tempDomain1[0,-1],y2bar,label='10,000 Cells')
plt.plot(tempDomain3[int(yNodes3)-1,:]/tempDomain1[0,-1],y3bar,label='19,600 Cells')
plt.plot(tempDomain4[int(yNodes4)-1,:]/tempDomain1[0,-1],y4bar,label='40,000 Cells')
plt.plot(tempDomain5[int(yNodes5)-1,:]/tempDomain1[0,-1],y5bar,label='80,000 Cells')
plt.title('temp at outflow')
plt.legend()
plt.grid(True)
plt.xlabel('u-vel')
plt.ylabel('Normalized y')
plt.show()

# Plotting at plate
plt.plot(x1,presDomain1[:,0]/presDomain1[0,-1],label='4,900 Cells')
plt.plot(x2,presDomain2[:,0]/presDomain2[0,-1],label='10,000 Cells')
plt.plot(x3,presDomain3[:,0]/presDomain3[0,-1],label='19,600 Cells')
plt.plot(x4,presDomain4[:,0]/presDomain4[0,-1],label='40,000 Cells')
plt.plot(x5,presDomain5[:,0]/presDomain5[0,-1],label='80,000 Cells')
plt.title('Pressure at Surface')
plt.legend()
plt.grid(True)
plt.xlabel('x-location on plate')
plt.ylabel('nondimensionalized pressure')
plt.show()






psum1 = (presDomain1[:,0].sum()*deltax1)
psum2 = (presDomain2[:,0].sum()*deltax2)
psum3 = (presDomain3[:,0].sum()*deltax3)
psum4 = (presDomain4[:,0].sum()*deltax4)
psum5 = (presDomain5[:,0].sum()*deltax5)

# change in pressure on plate N/m^2
header('pressure on plate')
d1 = psum2-psum1
d2 = psum3-psum2
d3 = psum4-psum3
d4 = psum5-psum4
print(abs(d1))
print(abs(d2))
print(abs(d3))
print(abs(d4))

#max pressures
header('max pressure')
print(presDomain1.max())
print(presDomain2.max())
print(presDomain3.max())
print(presDomain4.max())
print(presDomain5.max())


