import numpy as np 
import matplotlib.pyplot as plt

################
#mu Calculation#
################
def calcMu(muRef,temp,tempRef):
	"""
	Calculate local dynamic viscosity 
	using Sutherland's Law
	"""
	mu = muRef* ( (temp/tempRef)**(3/2) ) * ( (tempRef+110)/(temp+110) )
	return mu

	
############################
#Convert to U, F, E Vectors#
############################
def getUfromPrim(presDomain,tempDomain,uVelDomain,vVelDomain,muDomain,R,cv):
	U1 = presDomain/(R*tempDomain)
	U2 = U1 * uVelDomain
	U3 = U1 * vVelDomain
	U5 = U1*(cv*tempDomain + ( (uVelDomain**2 + vVelDomain**2)**0.5) /2)
	return U1,U2,U3,U5
#####
# Tau Calculations
def getTauxyfromPrim(uVelDomain,vVelDomain,muDomain,case,deltay,deltax):
	xMax,yMax = uVelDomain.shape
	tauxy = np.zeros(uVelDomain.shape)
	# Case 1: predictor for E
	if case ==1:
		for i in range(0,xMax):
			for j in range(0,yMax):
				# Calculate dudy
				if j==0:# Wall
					dudy = (uVelDomain[i,j+1]-uVelDomain[i,j])/deltay
				elif j == yMax-1: # top
					dudy = (uVelDomain[i,j]-uVelDomain[i,j-1])/deltay
				else:# j values between top and bottom
					dudy = (uVelDomain[i,j+1]-uVelDomain[i,j-1])/(2*deltay)
				# Calculate dvdx
				if i==0: # left inlet
					dvdx = (vVelDomain[i+1,j]-vVelDomain[i,j])/deltax
				else: # middle and right side
					dvdx = (vVelDomain[i,j]-vVelDomain[i,j-1])/deltax
				tauxy[i,j] = muDomain[i,j]*(dudy + dvdx)
	# Case 2: predictor for F
	elif case == 2:
		for i in range(0,xMax):
			for j in range(0,yMax):
				# Calculate dudy
				if j ==0: # if we are at the wall then we cannot do rearward diff
					dudy = (uVelDomain[i,j+1]-uVelDomain[i,j])/deltay
				else: #rearward diff
					dudy = (uVelDomain[i,j] - uVelDomain[i,j-1])/deltay
				# Calculate dvdx
				if i ==0:
					dvdx = (vVelDomain[i+1,j]-vVelDomain[i,j])/deltax
				elif i == xMax-1:
					dvdx = (vVelDomain[i,j]-vVelDomain[i-1,j])/deltax
				else:
					dvdx = (vVelDomain[i+1,j]-vVelDomain[i-1,j])/(2*deltax)
				tauxy[i,j] = muDomain[i,j]*(dudy + dvdx)
	# Case 3: corrector for E
	elif case == 3:
		for i in range(0,xMax):
			for j in range(0,yMax):
				# Calculate dudy
				if j==0:# Wall
					dudy = (uVelDomain[i,j+1]-uVelDomain[i,j])/deltay
				elif j == yMax-1: # top
					dudy = (uVelDomain[i,j]-uVelDomain[i,j-1])/deltay
				else:# j values between top and bottom
					dudy = (uVelDomain[i,j+1]-uVelDomain[i,j-1])/(2*deltay)
				# Calculate dvdx
				if i==xMax-1: # outlet
					dvdx = (vVelDomain[i,j]-vVelDomain[i-1,j])/deltax
				else: # middle and left side
					dvdx = (vVelDomain[i+1,j]-vVelDomain[i,j])/deltax
				tauxy[i,j] = muDomain[i,j]*(dudy + dvdx)
	# Case 4: corrector for F
	elif case == 4:
		for i in range(0,xMax):
			for j in range(0,yMax):
				# Calculate dudy
				if j ==yMax-1: # if we are at the wall then we cannot do rearward diff
					dudy = (uVelDomain[i,j]-uVelDomain[i,j-1])/deltay
				else: #rearward diff
					dudy = (uVelDomain[i,j+1] - uVelDomain[i,j])/deltay
				# Calculate dvdx
				if i ==0:
					dvdx = (vVelDomain[i+1,j]-vVelDomain[i,j])/deltax
				elif i == xMax-1:
					dvdx = (vVelDomain[i,j]-vVelDomain[i-1,j])/deltax
				else:
					dvdx = (vVelDomain[i+1,j]-vVelDomain[i-1,j])/(2*deltax)
				tauxy[i,j] = muDomain[i,j]*(dudy + dvdx)
	else:
		print('ERROR:No case selected. Must be 1, 2, 3, or 4')
	print(tauxy)
	return tauxy
def getTauxx():
	# Calculate Tauxx for predictor
	
	# Calculate Tauxx for corrector
	return

def getTauyy():
	# Calculate Tauxx for predictor
	
	# Calculate Tauxx for corrector
	return
	
	
	
#####################
#Boundary Conditions#
#####################
def initBoundaryCond(presDomain,tempDomain,uVelDomain,vVelDomain,tempWall):
	"""
	Calculation of boundary conditions
	"""
	# Case 1. [0,0] leading edge
	uVelDomain[0,0] = 0
	vVelDomain[0,0] = 0
	# Case 2. [0,1:] inflow and [:.-1] top 
	vVelDomain[0,1:]  = 0  # not including leading edge
	vVelDomain[:,-1]  = 0
	# Case 3. [:,0] bottom (flat plate wall)
	presDomain[1:,0] = 2*presDomain[1:,1] - presDomain[1:,2]
	tempDomain[1:,0] = tempWall
	uVelDomain[1:,0] = 0
	vVelDomain[1:,0] = 0	
	# Case 4. [-1,1:-1 Outflow
	presDomain[-1,1:-1] = 2*presDomain[-1,1:-1] - presDomain[-1,1:-1]
	tempDomain[-1,1:-1] = 2*tempDomain[-1,1:-1] - tempDomain[-1,1:-1]
	uVelDomain[-1,1:-1] = 2*uVelDomain[-1,1:-1] - uVelDomain[-1,1:-1]
	vVelDomain[-1,1:-1] = 2*vVelDomain[-1,1:-1] - vVelDomain[-1,1:-1]
	return presDomain,tempDomain,uVelDomain,vVelDomain

########################
#Time Step Calculations#
########################
def calcTimeStep(presDomain,tempDomain,uVelDomain,vVelDomain,muDomain,deltax,deltay,gamma,Pr,R,K):
	rhoDomain = presDomain/(R*tempDomain)
	vPrime =( (4/3) * muDomain * (gamma * muDomain/Pr) ) / rhoDomain
	vPrime = vPrime.max()
	a = (gamma*R*tempDomain)**0.5
	deltatCFL =( (abs(uVelDomain)/deltax) + (abs(vVelDomain)/deltay) +
		a * ( (1/deltax**2) + (1/deltay**2) )**0.5 
		+ 2* vPrime*((1/deltax**2)+(1/deltay**2)) ) **-1
	deltat = K*deltatCFL.min()
	return deltat

##################
#Iteration Scheme#
##################
def solveFlow(presDomain,tempDomain,uVelDomain,vVelDomain,muDomain,U1,U2,U3,U5,deltay,deltax):
	"""
	We will be using the MacCormack method to solve the 
	Navier-Stokes equation. 
	"""
	# Calculate tau for predictor both F and E
	tauxyE = getTauxyfromPrim(uVelDomain,vVelDomain,muDomain,1,deltay,deltax)
	tauxyF = getTauxyfromPrim(uVelDomain,vVelDomain,muDomain,2,deltay,deltax)
	#tauxxE = 1
	#tauyyF = 2
	#E1,E2,E3,E5 = getE()
	#F1,F2,F3,F5 = getF()
	
	
	
	# Calculate tau for corrector both F and E
	tauxyE = getTauxyfromPrim(uVelDomain,vVelDomain,muDomain,3,deltay,deltax)
	tauxyF = getTauxyfromPrim(uVelDomain,vVelDomain,muDomain,4,deltay,deltax)
	#tauxxE = 1
	#tauyyF = 2
	#E1,E2,E3,E5 = getE()
	#F1,F2,F3,F5 = getF()
	return presDomain,tempDomain,uVelDomain,vVelDomain

#############
#MAIN SCRIPT#
#############
def superVisc(girdPtsX,girdPtsY,machInf,TwTInf,iters,K):
	"""
	This program uses an iterative method to solve
	the full Navier-Stokes equation for super sonic
	flow over a flat plate. 
	Input Parameters:
	girdPtsX    - Nodes in x
	girdPtsY    - Nodes in y
	machInf     - Inlet mach number
	TwTInf      - Relationship of wall temp over free stream 
	iters       - iterations to complete
	K           - Courant number
	**TO DO: set value for desired residual remove iters**
	"""

	# Set flow/geometry parameters
	lengthOfPlate =0.00001
	aInf         = 340.28        # Speed of sound: m/s
	pressureInf  = 101325.0      # Pressure N/m^2
	tempInf      = 288.16        # Kelvin
	gamma        = 1.4           # Ratio of specific heats
	Pr           = 0.71          # Prandtl number
	muRef        = 1.7894*10**-5  # Dynamic viscosity reference kg/(m*s)
	tempRef      = 288.16        # Reference temperature Kelvin
	R            = 287           # Gas constant Number J/(kg*K)
	cv           = R/(gamma-1)
	cp           = gamma * cv
	rhoInf       = pressureInf/(R*tempInf)
	uVelInf      = machInf*aInf
	vVelInf      = 0 
	muInf        = calcMu(muRef,tempInf,tempRef)
	ReL          = rhoInf*uVelInf*lengthOfPlate/muInf
	eInf         = cv*tempInf
	kInf         = muInf*cp/Pr
	tempWall     = TwTInf*tempInf
	# Calculate delta x
	deltax = lengthOfPlate/(girdPtsX-1)

	# Calculate delta y
	delta = (5*lengthOfPlate)/(ReL**.5)
	heightOfGrid = 5*delta
	deltay = heightOfGrid/(girdPtsY-1)
	print('Delta in x: ',deltax,'\nDelta in y: ',deltay)

	# Initialize free stream
	presDomain = np.ones([girdPtsX,girdPtsY])*pressureInf
	tempDomain = np.ones([girdPtsX,girdPtsY])*tempInf
	uVelDomain = np.ones([girdPtsX,girdPtsY])*uVelInf
	vVelDomain = np.ones([girdPtsX,girdPtsY])*vVelInf
	# Initialize BC's
	presDomain,tempDomain,uVelDomain,vVelDomain = initBoundaryCond(
		presDomain,tempDomain,uVelDomain,vVelDomain,tempWall)

	# Visualize initial conditions
	#plotFlow(presDomain,tempDomain,uVelDomain,vVelDomain)

	# Begin Simulations 
	for i in range(iters):
		# Calculate time step
		muDomain = calcMu(muRef,tempDomain,tempRef)
		deltat = calcTimeStep(presDomain,tempDomain,uVelDomain,
			vVelDomain,muDomain,deltax,deltay,gamma,Pr,R,K)

		
		# Calculate terms for differential equation
		U1,U2,U3,U5 = getUfromPrim(
			presDomain,tempDomain,uVelDomain,vVelDomain,muDomain,R,cv)
		#Solve internal points
		presDomain,tempDomain,uVelDomain,vVelDomain = solveFlow(
			presDomain,tempDomain,uVelDomain,vVelDomain,muDomain
			,U1,U2,U3,U5,deltay,deltax)
 
		# Set BC's
		
		# iterate 

		# check for convergence

	return

####################
#Plotting Functions#
####################
def plotFlow(pres,temp,uVel,vVel):
	"""
	Plotting script to show flowfiled 
	"""
	plt.figure()
	plt.subplot(2,2,1)
	plt.contourf(pres.transpose())
	plt.colorbar()
	plt.subplot(2,2,2)
	plt.contourf(temp.transpose())
	plt.colorbar()
	plt.subplot(2,2,3)
	plt.contourf(uVel.transpose())
	plt.colorbar()
	plt.subplot(2,2,4)
	plt.contourf(vVel.transpose())
	plt.colorbar()
	plt.show()
	return
	
def plotParam(param):
	plt.figure()
	plt.contourf(param.transpose())
	plt.colorbar()
	plt.show()
	return

superVisc(70,70,4,1,1,0.6)