##
# Time evolution of shallow water
# A real-time evolving simulation of water circulation under different environments.
# The "shallow water equations" are a set of hyperbolic partial differential equations
# that describe the flow below a pressure surface in a fluid. 
# The simulation shows the elevations and direction of the flow of water over time, 
# initially starting on a "hill". The flow is altered by rotation (as on a rotating planet), 
# by wind, and by friction. 
##


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tkr
import math

# Grid and Variable Initialization 
ncol = 5
nrow = ncol

nSlices = 1000			# maximum number of frames to show in the plot
ntAnim = 1			# number of time steps for each frame

horizontalWrap = True 
interpolateRotation = False			# or True, either works
textOutput = False
plotOutput = True
arrowScale = 30
rotationScheme = "PlusMinus"			# "WithLatitude", "PlusMinus", "Uniform"

windScheme = ""			# "Curled", "Uniform"
initialPerturbation = "Tower"			# "Tower", "NSGradient", "EWGradient"

dT = 600			# seconds
G = 9.8e-4			# artificially low to allow a long time step
HBackground = 4000			# meters

dX = 1.E4			# meters, small ocean on a small low-G planet.  
dY = dX

dxDegrees = dX / 110.e3
flowConst = G			# 1/s2
dragConst = 1.E-6			# about 10 days decay time
meanLatitude = 30			# degrees

latitude = []
rotConst = []
windU = []
for irow in range(0,nrow):
    if rotationScheme is "WithLatitude":
        latitude.append( meanLatitude + (irow - nrow/2) * dxDegrees )
        rotConst.append( -7.e-5 * math.sin(math.radians(latitude[-1])))
    elif rotationScheme is "PlusMinus":
        rotConst.append( -3.5e-5 * (1. - 0.8 * ( irow - (nrow-1)/2 ) / nrow ))
    elif rotationScheme is "Uniform":
        rotConst.append( -3.5e-5 ) 
    else:
        rotConst.append( 0 )

    if windScheme is "Curled":
        windU.append( 1e-8 * math.sin( (irow+0.5)/nrow * 2 * 3.14 ) ) 
    elif windScheme is "Uniform":
        windU.append( 1.e-8 )
    else:
        windU.append( 0 )

itGlobal = 0

U = np.zeros((nrow, ncol+1))
V = np.zeros((nrow+1, ncol))
H = np.zeros((nrow, ncol+1))

dUdT = np.zeros((nrow, ncol))
dVdT = np.zeros((nrow, ncol))
dHdT = np.zeros((nrow, ncol))

dHdX = np.zeros((nrow, ncol+1))
dHdY = np.zeros((nrow, ncol))
dUdX = np.zeros((nrow, ncol))
dVdY = np.zeros((nrow, ncol))

rotV = np.zeros((nrow,ncol))
rotU = np.zeros((nrow,ncol))

tempU = np.zeros((nrow, ncol))
tempV = np.zeros((nrow, ncol))
    
midCell = int(ncol/2)
if initialPerturbation is "Tower":
    H[midCell,midCell] = 1
elif initialPerturbation is "NSGradient":
    H[0:midCell,:] = 0.1
elif initialPerturbation is "EWGradient":
    H[:,0:midCell] = 0.1

def animStep():    

    global stepDump, itGlobal

    # Boundary Conditions
    V[nrow,:] = 0
    V[0,:] = 0

    if horizontalWrap is True:
        U[:,ncol] = U[:,0]
        H[:,ncol] = H[:,0]

    else:
        U[:,0] = 0
        U[:,ncol] = 0

    # Time Loop
    for it in range(0,ntAnim):
        for icol in range(0, ncol):
            for irow in range(0, nrow):
                # Longitudinal derivatives
                dHdX[irow, icol] = (H[irow, icol] - H[irow, icol-1]) / dX
                dUdX[irow, icol] = (U[irow, icol+1] - U[irow, icol]) / dX

                # Latitudinal derivatives
                dHdY[irow, icol] = (H[irow, icol] - H[irow-1, icol]) / dY
                dVdY[irow, icol] = (V[irow+1, icol] - V[irow, icol]) / dY

                # Rotational terms
                if interpolateRotation == True:
                    # Take average to get velocities at cell centers
                    tempU[irow, icol] = (U[irow, icol] + U[irow, icol+1]) / 2 * rotConst[irow]
                    tempV[irow, icol] = (V[irow, icol] + V[irow+1, icol]) / 2 * rotConst[irow]

                    # Vertical rotational component for V: average rotU and interpolate vertically
                    rotU[irow, icol] = (tempU[irow, icol] + tempU[irow-1, icol]) / 2
                    # Horizontal rotational component for U: average rotV and interpolate horizontally
                    rotV[irow, icol] = (tempV[irow, icol] + tempV[irow, icol-1]) / 2

                else:
                    rotU[irow,icol] = rotConst[irow] * U[irow,icol]
                    rotV[irow,icol] = rotConst[irow] * V[irow,icol]

                # Time derivatives
                dUdT[irow,icol] = rotV[irow,icol] - flowConst * dHdX[irow,icol] - dragConst * U[irow,icol] + windU[irow]
                dVdT[irow,icol] = -rotU[irow,icol] - flowConst * dHdY[irow,icol] - dragConst * V[irow,icol]
                dHdT[irow,icol] = -(dUdX[irow,icol] + dVdY[irow,icol]) * HBackground / dX
                
                # Step forward in time
                U[irow, icol] += dUdT[irow, icol] * dT
                V[irow, icol] += dVdT[irow, icol] * dT
                H[irow, icol] += dHdT[irow, icol] * dT
                
    itGlobal = itGlobal + ntAnim

def firstFrame():
    global fig, ax, hPlot
    fig, ax = plt.subplots()
    ax.set_title("H")   
    hh = H[:,0:ncol]
    loc = tkr.IndexLocator(base=1, offset=1)
    ax.xaxis.set_major_locator(loc)
    ax.yaxis.set_major_locator(loc)
    grid = ax.grid(which='major', axis='both', linestyle='-')
    hPlot = ax.imshow(hh, interpolation='nearest', clim=(-0.5,0.5))   
    plotArrows()
    plt.show(block=False) 

def plotArrows():
    global quiv, quiv2
    xx = []
    yy = []
    uu = []
    vv = []
    for irow in range( 0, nrow ):
        for icol in range( 0, ncol ):
            xx.append(icol - 0.5)
            yy.append(irow )
            uu.append( U[irow,icol] * arrowScale )
            vv.append( 0 )
    quiv = ax.quiver( xx, yy, uu, vv, color='white', scale=1)
    for irow in range( 0, nrow ):
        for icol in range( 0, ncol ):
            xx.append(icol)
            yy.append(irow - 0.5)
            uu.append( 0 )
            vv.append( -V[irow,icol] * arrowScale )
    quiv2 = ax.quiver( xx, yy, uu, vv, color='white', scale=1)

def updateFrame():
    global fig, ax, hPlot, quiv, quiv2
    hh = H[:,0:ncol]
    hPlot.set_array(hh)
    quiv.remove()    
    quiv2.remove()
    plotArrows()
    plt.pause(0.00001)
    fig.canvas.draw()
    print("Time: ", math.floor( itGlobal * dT / 86400.*10)/10, "days")

def textDump():
    print("time step ", itGlobal)    
    print("H", H)
    print("dHdX" )
    print( dHdX)
    print("dHdY" )
    print( dHdY)
    print("U" )
    print( U)
    print("dUdX" )
    print( dUdX)
    print("rotV" )
    print( rotV)
    print("V" )
    print( V)
    print("dVdY" )
    print( dVdY)
    print("rotU" )
    print( rotU)
    print("dHdT" )
    print( dHdT)
    print("dUdT" )
    print( dUdT)
    print("dVdT" )
    print( dVdT)
    print('windU')
    print(windU)

if textOutput is True:
    textDump()
if plotOutput is True:
    firstFrame()
for i_anim_step in range(0,nSlices):
    animStep()
    if textOutput is True:
        textDump()
    if plotOutput is True:
        updateFrame()