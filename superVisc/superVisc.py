import numpy as np 
import matplotlib.pyplot as plt
import pdb
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
	U5 = U1 * (cv*tempDomain + (uVelDomain**2 + vVelDomain**2) / 2 )
	return U1,U2,U3,U5

def U2Prim(U1Bar,U2Bar,U3Bar,U5Bar,R,cv):
	rho = U1Bar    
	uVelDomain = U2Bar/U1Bar
	vVelDomain = U3Bar/U1Bar
	e =   (U5Bar/U1Bar) - ( (uVelDomain**2) + (vVelDomain**2) )/2 
	tempDomain = e/cv
	presDomain = U1Bar*R*tempDomain
	return presDomain,tempDomain,uVelDomain,vVelDomain

def getE(tauxyE,tauxyF,tauxxE,tauyyF,presDomain,tempDomain,uVelDomain,vVelDomain,qxE,R,cv):
	rho = presDomain/(R*tempDomain)
	E1 = uVelDomain*rho
	E2 = rho * uVelDomain**2 + presDomain - tauxxE
	E3 = rho * uVelDomain * vVelDomain - tauxyE
	Et = rho*(cv*tempDomain + ( (uVelDomain**2 + vVelDomain**2)**0.5) /2)
	E5 = (Et + presDomain) * uVelDomain - uVelDomain*tauxxE - vVelDomain*tauxyE + qxE
	return E1,E2,E3,E5
def getF(tauxyE,tauxyF,tauxxE,tauyyF,presDomain,tempDomain,uVelDomain,vVelDomain,qyF,R,cv):
	rho = presDomain/(R*tempDomain)
	F1 = vVelDomain*rho
	F2 = rho * uVelDomain* vVelDomain -tauxyF
	F3 = rho * vVelDomain**2 + presDomain - tauyyF 
	Et = rho*(cv*tempDomain + ( (uVelDomain**2 + vVelDomain**2)**0.5) /2)
	F5 = (Et + presDomain) * vVelDomain - uVelDomain*tauxyF - vVelDomain*tauyyF + qyF
	return F1,F2,F3,F5

#########################
# Tau and q calculations#
#########################
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
	return tauxy
def getTauxx(uVelDomain,vVelDomain,muDomain,case,deltay,deltax):
	# Calculate Tauxx for predictor
	xMax,yMax = uVelDomain.shape
	tauxx = np.zeros(uVelDomain.shape)
	# Case 1: predictor for E
	if case ==1:
		for i in range(0,xMax):
			for j in range(0,yMax):
				# Calculate dvdy with central diff
				if j==0:# Wall
					dvdy = (vVelDomain[i,j+1]-vVelDomain[i,j])/deltay
				elif j == yMax-1: # top
					dvdy = (vVelDomain[i,j]-vVelDomain[i,j-1])/deltay
				else:# j values between top and bottom
					dvdy = (vVelDomain[i,j+1]-vVelDomain[i,j-1])/(2*deltay)

				# Calculate dudx with rearward diff
				if i==0: # left inlet
					dudx = (uVelDomain[i+1,j]-uVelDomain[i,j])/deltax
				else: # middle and right side
					dudx = (uVelDomain[i,j]-uVelDomain[i-1,j])/deltax
				tauxx[i,j] = -(2/3) * muDomain[i,j] * (dudx + dvdy) + 2 * muDomain[i,j] * dudx
	# Case 2: corrector for E
	elif case == 2:
		for i in range(0,xMax):
			for j in range(0,yMax):
				# Calculate dvdy with central diff
				if j==0:# Wall
					dvdy = (vVelDomain[i,j+1]-vVelDomain[i,j])/deltay
				elif j == yMax-1: # top
					dvdy = (vVelDomain[i,j]-vVelDomain[i,j-1])/deltay
				else:# j values between top and bottom
					dvdy = (vVelDomain[i,j+1]-vVelDomain[i,j-1])/(2*deltay)

				# Calculate dudx with forward diff
				if i==xMax-1: # right outlet
					dudx = (uVelDomain[i,j]-uVelDomain[i-1,j])/deltax
				else: # middle and right side
					dudx = (uVelDomain[i+1,j]-uVelDomain[i,j])/deltax
				tauxx[i,j] = -(2/3) * muDomain[i,j] * (dudx + dvdy) + 2 * muDomain[i,j] * dudx
	# Calculate Tauxx for corrector
	return tauxx
def getTauyy(uVelDomain,vVelDomain,muDomain,case,deltay,deltax):
	xMax,yMax = uVelDomain.shape
	tauyy = np.zeros(uVelDomain.shape)
	# Calculate tauyy for predictor
	if case == 1:
		for i in range(0,xMax):
			for j in range(0,yMax):
				# Calculate dvdy with rearward
				if j ==0: # if we are at the wall then we cannot do rearward diff
					dvdy = (vVelDomain[i,j+1]-vVelDomain[i,j])/deltay
				else: #rearward diff
					dvdy = (vVelDomain[i,j] - vVelDomain[i,j-1])/deltay

				# Calculate dudx with central diff
				if i ==0:
					dudx = (uVelDomain[i+1,j]-uVelDomain[i,j])/deltax
				elif i == xMax-1:
					dudx = (uVelDomain[i,j]-uVelDomain[i-1,j])/deltax
				else:
					dudx = (uVelDomain[i+1,j]-uVelDomain[i-1,j])/(2*deltax)
				tauyy[i,j] = -(2/3) * muDomain[i,j] * (dudx + dvdy) + 2 * muDomain[i,j] * dvdy

	# Calculate tauyy for corrector
	elif case == 2:
		for i in range(0,xMax):
			for j in range(0,yMax):
				# Calculate dvdy with forward diff
				if j ==yMax-1: # outlet
					dvdy = (vVelDomain[i,j]-vVelDomain[i,j-1])/deltay
				else: #forward diff
					dvdy = (vVelDomain[i,j+1] - vVelDomain[i,j])/deltay

				# Calculate dudx with central diff
				if i ==0:
					dudx = (uVelDomain[i+1,j]-uVelDomain[i,j])/deltax
				elif i == xMax-1:
					dudx = (uVelDomain[i,j]-uVelDomain[i-1,j])/deltax
				else:
					dudx = (uVelDomain[i+1,j]-uVelDomain[i-1,j])/(2*deltax)
				tauyy[i,j] = -(2/3) * muDomain[i,j] * (dudx + dvdy) + 2 * muDomain[i,j] * dvdy
	return tauyy
def getqxE(tempDomain,cp,muDomain,deltax,case):
	xMax,yMax = tempDomain.shape
	qxE = np.zeros(tempDomain.shape)
	# Case 1: predictor for E
	if case ==1:
		for i in range(0,xMax):
			for j in range(0,yMax):
				# Calculate dudx with rearward diff
				if i==0: # left inlet
					dTdx = (tempDomain[i+1,j]-tempDomain[i,j])/deltax
				else: # middle and right side
					dTdx = (tempDomain[i,j]-tempDomain[i-1,j])/deltax
				qxE[i,j] = -((muDomain[i,j]*cp)/0.71)*dTdx			
	# Case 2: corrector for E
	elif case == 2:
		for i in range(0,xMax):
			for j in range(0,yMax):
				# Calculate dudx with forward diff
				if i==xMax-1: # right outlet
					dTdx = (tempDomain[i,j]-tempDomain[i-1,j])/deltax
				else: # middle and right side
					dTdx = (tempDomain[i+1,j]-tempDomain[i,j])/deltax	
				qxE[i,j] = -((muDomain[i,j]*cp)/0.71)*dTdx	
	return qxE
def getqyF(tempDomain,cp,muDomain,deltay,case):
	xMax,yMax = tempDomain.shape
	qyF = np.zeros(tempDomain.shape)
	# Predictor step
	if case == 1:
		for i in range(0,xMax):
			for j in range(0,yMax):
				# Calculate dvdy with rearward
				if j ==0: # if we are at the wall then we cannot do rearward diff
					dTdy = (tempDomain[i,j+1]-tempDomain[i,j])/deltay
				else: #rearward diff
					dTdy = (tempDomain[i,j] - tempDomain[i,j-1])/deltay
				qyF[i,j] = -((muDomain[i,j]*cp)/0.71)*dTdy	
	# Corrector step calculation
	elif case == 2:
		for i in range(0,xMax):
			for j in range(0,yMax):
				# Calculate dvdy with forward diff
				if j ==yMax-1: # outlet
					dTdy = (tempDomain[i,j]-tempDomain[i,j-1])/deltay
				else: #forward diff
					dTdy = (tempDomain[i,j+1] - tempDomain[i,j])/deltay
				qyF[i,j] = -((muDomain[i,j]*cp)/0.71)*dTdy	
	return qyF
	
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

###############################
#Predictor and Corrector Steps#
###############################
def predictor(U1,U2,U3,U5,E1,E2,E3,E5,F1,F2,F3,F5,deltat,deltax,deltay):
	xMax,yMax = U1.shape
	U1Bar = np.zeros(U1.shape)
	U2Bar = np.zeros(U1.shape)
	U3Bar = np.zeros(U1.shape)
	U5Bar = np.zeros(U1.shape)
	# Loop through all innner nodes
	for i in range(0,xMax):
		for j in range(0,yMax):
			if i==0 or i==xMax-1 or j==0 or j==yMax-1:
				U1Bar[i,j] = U1[i,j]*1 
				U2Bar[i,j] = U2[i,j]*1
				U3Bar[i,j] = U3[i,j]*1
				U5Bar[i,j] = U5[i,j]*1
			else:
				U1Bar[i,j] = U1[i,j] - (deltat/deltax)*(E1[i+1,j]-E1[i,j]) - (deltat/deltay)*(F1[i,j+1]-F1[i,j])
				U2Bar[i,j] = U2[i,j] - (deltat/deltax)*(E2[i+1,j]-E2[i,j]) - (deltat/deltay)*(F2[i,j+1]-F2[i,j])
				U3Bar[i,j] = U3[i,j] - (deltat/deltax)*(E3[i+1,j]-E3[i,j]) - (deltat/deltay)*(F3[i,j+1]-F3[i,j])
				U5Bar[i,j] = U5[i,j] - (deltat/deltax)*(E5[i+1,j]-E5[i,j]) - (deltat/deltay)*(F5[i,j+1]-F5[i,j])
	return U1Bar,U2Bar,U3Bar,U5Bar
def corrector(U1,U2,U3,U5,U1Bar,U2Bar,U3Bar,U5Bar,E1Bar,E2Bar,E3Bar,E5Bar,F1Bar,F2Bar,F3Bar,F5Bar,deltat,deltax,deltay):
	U1,U2,U3,U5 = [1,2,3,5]
	return U1,U2,U3,U5
##################
#Iteration Scheme#
##################
def solveFlow(presDomain,tempDomain,uVelDomain,vVelDomain,muDomain,deltat,deltay,deltax,cp,cv,R,tempRef,muRef):
	"""
	We will be using the MacCormack method to solve the 
	Navier-Stokes equation. 
	"""
	################
	#Predictor Step#
	################
	# Calculate tau for predictor both F and E

	tauxyE = getTauxyfromPrim(uVelDomain,vVelDomain,muDomain,1,deltay,deltax) # predictor for E values (case1:getTauxyfromPrim)
	tauxyF = getTauxyfromPrim(uVelDomain,vVelDomain,muDomain,2,deltay,deltax) # predictor for F values (case2:getTauxyfromPrim)
	tauxxE = getTauxx(uVelDomain,vVelDomain,muDomain,1,deltay,deltax) # predictor for E values (case1:getTauxx)
	tauyyF = getTauyy(uVelDomain,vVelDomain,muDomain,1,deltay,deltax) # predictor for F values (case1:getTauyy)
	qxE = getqxE(tempDomain,cp,muDomain,deltax,1) # predictor for E values (case1:getqxE)
	qyF = getqyF(tempDomain,cp,muDomain,deltay,1) # predictor for F values (case1:getqyF)
	# Calculate the values for U, E and F arrays
	U1,U2,U3,U5 = getUfromPrim(presDomain,tempDomain,uVelDomain,vVelDomain,muDomain,R,cv)
	E1,E2,E3,E5 = getE(tauxyE,tauxyF,tauxxE,tauyyF,presDomain,tempDomain,uVelDomain,vVelDomain,qxE,R,cv)
	F1,F2,F3,F5 = getF(tauxyE,tauxyF,tauxxE,tauyyF,presDomain,tempDomain,uVelDomain,vVelDomain,qyF,R,cv)
	# Calculate predictor step
	U1Bar,U2Bar,U3Bar,U5Bar = predictor(U1,U2,U3,U5,E1,E2,E3,E5,F1,F2,F3,F5,deltat,deltax,deltay)
	# Extract primitives
	presDomain,tempDomain,uVelDomain,vVelDomain = U2Prim(U1Bar,U2Bar,U3Bar,U5Bar,R,cv)
	# Enforce BC's
	################
	#Corrector Step#
	################
	muDomain = muRef* ( (tempDomain/tempRef)**(3/2) ) * ( (tempRef+110)/(tempDomain+110) )
	# muDomain = calcMu(muRef,tempDomain,tempRef) # Recalculate mu base on predicted temp
	tauxyE = getTauxyfromPrim(uVelDomain,vVelDomain,muDomain,3,deltay,deltax) # corrector for E values (case3:getTauxyfromPrim)
	tauxyF = getTauxyfromPrim(uVelDomain,vVelDomain,muDomain,4,deltay,deltax) # corrector for F values (case4:getTauxyfromPrim)
	tauxxE = getTauxx(uVelDomain,vVelDomain,muDomain,2,deltay,deltax) # corrector for E values (case2:getTauxx)
	tauyyF = getTauyy(uVelDomain,vVelDomain,muDomain,2,deltay,deltax) # corrector for F values (case2:getTauyy)
	qxE = getqxE(tempDomain,cp,muDomain,deltax,2)
	qyF = getqyF(tempDomain,cp,muDomain,deltay,2)
	U1Bar,U2Bar,U3Bar,U5Bar = getUfromPrim(presDomain,tempDomain,uVelDomain,vVelDomain,muDomain,R,cv)
	E1Bar,E2Bar,E3Bar,E5Bar = getE(tauxyE,tauxyF,tauxxE,tauyyF,presDomain,tempDomain,uVelDomain,vVelDomain,qxE,R,cv)
	F1Bar,F2Bar,F3Bar,F5Bar = getF(tauxyE,tauxyF,tauxxE,tauyyF,presDomain,tempDomain,uVelDomain,vVelDomain,qyF,R,cv)
	
	# corrector(U1,U2,U3,U5,U1Bar,U2Bar,U3Bar,U5Bar,E1Bar,E2Bar,E3Bar,E5Bar,F1Bar,F2Bar,F3Bar,F5Bar,deltat,deltax,deltay)


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
	print('GRID DIMENSIONS','\nDelta in x: ',deltax,'\nDelta in y: ',deltay)

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

		#Solve internal points
		print('Iter: ',i,' | Max Pressure: ', presDomain.max(),' | Max Temperature: ',
			tempDomain.max(),' | Max X-Vel: ',uVelDomain.max(),' | Max Y_Vel: ',
			vVelDomain.max())
		presDomain,tempDomain,uVelDomain,vVelDomain = solveFlow(presDomain,tempDomain,uVelDomain,vVelDomain,muDomain,deltat,deltay,deltax,cp,cv,R,tempRef,muRef)
 
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


#superVisc(girdPtsX,girdPtsY,machInf,TwTInf,iters,K)
superVisc(70,70,4,1,1,0.6)