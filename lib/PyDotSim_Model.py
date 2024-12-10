"""

Model geomtry definitions for PyDotSim

by Christian Heyn

"""


import numpy as np
import math
import lib.PyDotSim_Constants as cnst

T =         4
xAl =       0.34
statesNum = 3
WFflat_threshold = 0.02

MultiParticleLoops  = 1

# colMap = 'gnuplot2'
# colMap = 'jet'
colMap = 'gray'

QDshapeStr =        ['Ellipsoid', 'Disk','Cone','Cone-shell']
QDshapeParaStr =    ['rx,  ry,  rz', 'r,  h,  --', 'r,  h,  --', 'r,  hH,  hQ']
QDpara =            [1, 5e-9, 5e-9, 5e-9]           # QDpara:   [QD shape, r1, r2, r3]
cellNum1D =         61
cellSize =          (2*QDpara[1] + 6e-9)/cellNum1D
cellVol =           cellSize*cellSize*cellSize
meQ, mhQ =          cnst.AlGaAs_effMass(0)
meB, mhB =          cnst.AlGaAs_effMass(xAl)
smoothMaterialBorders = False
EFields =           [0, 0]                          # EFields [Fz, Fx]
BZfield =           [0, 0, 0]                       # BFields [Bz, le, lh]
PointCharges =      []                              # Point charges

MultiParticleStr =  ['X', 'XX', 'X+', 'X-']
MPpara =            [0, 0]                          # MPpara:   [Particle, self-consistent loops]

inputStr =      [ [ 'Size [nm]', 'r1', 'r2', 'r3' ],
                  [ 'Cells num (1D)', 'N'],
                  [ 'Cell size [nm]', '' ],
                  [ 'E fields [V/m]', 'Fz', 'Fx' ],
                  [ 'Bz field', 'Bz[T]', 'le' , 'lh']    ]

inputVals =     []
seriesStr =     ['r1', 'r2', 'r3', 'Fz', 'Fx', 'Bz']
modelCreated =  False

material = []
Ve = []; Vh = []
masse = []; massh = []


def calcMaterialPara():
    global meQ, mhQ, meB, mhB
    meQ, mhQ = cnst.AlGaAs_effMass(0)
    meB, mhB = cnst.AlGaAs_effMass(xAl)


def setInputVals():
    inputVals[0].set(QDpara[1]*1e9)
    inputVals[1].set(QDpara[2]*1e9)
    inputVals[2].set(QDpara[3]*1e9)
    inputVals[3].set(cellNum1D)
    inputVals[4].set( '{0: >#6.4f}'.format(cellSize*1e9) )
    inputVals[5].set(EFields[0])
    inputVals[6].set(EFields[1])
    inputVals[7].set(BZfield[0])
    inputVals[8].set(BZfield[1])
    inputVals[9].set(BZfield[2])


def readInputVals():
    global QDpara, cellNum1D, cellSize, meQ, mhQ, meB, mhB, EFields
    QDpara[1] =  float(inputVals[0].get())/1e9
    QDpara[2] =  float(inputVals[1].get())/1e9
    QDpara[3] =  float(inputVals[2].get())/1e9
    cellNum1D =    int(inputVals[3].get())
    cellSize  =  float(inputVals[4].get())/1e9
    EFields[0] = float(inputVals[5].get())
    EFields[1] = float(inputVals[6].get())
    BZfield[0] = float(inputVals[7].get())
    BZfield[1] = float(inputVals[8].get())
    BZfield[2] = float(inputVals[9].get())



# returns the xyz position of a mesh cell with index cellIndx
def cellPos(cellIndx):
    x0 = -cellSize*(cellNum1D-1)/2
    y0 = -cellSize*(cellNum1D-1)/2
    z0 = -cellSize*(cellNum1D-1)/2
    cx = x0 + cellIndx[0]*cellSize
    cy = y0 + cellIndx[1]*cellSize
    cz = z0 + cellIndx[2]*cellSize
    return cx, cy, cz


x = []; y = []; z =[]
def createXYZ():
    global x, y, z, cellSize
    #cellSize =  (2*QDpara[1] + 6e-9)/cellNum[0]
    for i in range(cellNum1D):
        cx, cy, cz = cellPos([i, 0, 0])
        x.append(cx)
    for i in range(cellNum1D):
        cx, cy, cz = cellPos([0, i, 0])
        y.append(cy)
    for i in range(cellNum1D):
        cx, cy, cz = cellPos([0, 0, i])
        z.append(cz)


# returns True if the mesh cell with position pos is inside the QD
# QDpara: [QD shape, r1, r2, r3, ...]
def cellQD(QDpara, pos):
    inQD = False
    QDshape = QDpara[0]
    cx = pos[0]; cy = pos[1]; cz = pos[2]
    # Ellipse, so far only sphere
    # idx_QD = ((X-x0)/Rx).^2 + ((Y-y0)/Ry).^2 + ((Z-z0)/Rz).^2 <= 1 ;
    if QDshape == 1:
        rxQD = QDpara[1]; ryQD = QDpara[2]; rzQD = QDpara[3];
        cr = np.sqrt(cx**2 + cy**2 + cz**2);
        if cr <= rxQD: inQD = True
    # Disk, r1: r, r2: h
    if QDshape == 2:
        rQD = QDpara[1]; hQD = QDpara[2]
        h0 = hQD/2
        if abs(cz) <= h0:
            cr = np.sqrt(cx**2+cy**2)
            if cr <= rQD: inQD = True
    # Cone, r1: r, r2: h
    if QDshape == 3:
        rQD = QDpara[1]; hQD = QDpara[2]
        h0 = hQD/2
        if abs(cz) <= h0:
            rmax = (cz+h0)*rQD/hQD
            cr = np.sqrt(cx**2+cy**2)
            if cr <= rmax: inQD = True
     # Cone shell, r1: r, r2: hHole, r3: hQd
    if QDshape == 4:
        rQD = QDpara[1]; hHole = QDpara[2]; hQD = QDpara[3]
        h0 = hHole/2
        if abs(cz) <= h0:
            rmax = (cz+h0)*rQD/hHole
            cr = np.sqrt(cx**2+cy**2)
            if cr <= rmax:
                zmax = hQD-h0 + (hHole-hQD)*cr/rQD
                if cz <= zmax: inQD = True
    return inQD


# returns the material composition beween 0=QD and 1=Barrier
def cellMaterial(cellIndx):
    cx, cy, cz = cellPos(cellIndx)
    cmaterial = 1
    if cellQD(QDpara, [cx, cy, cz]): cmaterial = 0
    return cmaterial


# creates the model with:
# material: 0=QD or 1=Barrier
# potentials Ve, Vh
# effective masses: masse, massh
# add fields
def createModel():
    global cellVol, material, Ve, Vh, masse, massh, modelCreated
    cellVol = cellSize*cellSize*cellSize
    N = cellNum1D; N3 = N*N*N
    createXYZ()

    material = np.zeros(N3).reshape(N,N,N)
    masse = np.zeros(N3).reshape(N,N,N)
    massh = np.zeros(N3).reshape(N,N,N)
    Ve = np.zeros(N3).reshape(N,N,N)
    Vh = np.zeros(N3).reshape(N,N,N)

    dEc, dEv = cnst.bandgapDiscontinuity(T, xAl)
    for ix in range(N):
        for iy in range(N):
            for iz in range(N):
                material[ix,iy,iz] = cellMaterial([ix,iy,iz])
    if smoothMaterialBorders:
        mGradientAll = np.gradient(material)
        mGradient = np.add(np.abs(mGradientAll[0]),np.abs(mGradientAll[1]))
        mGradient = np.add(mGradient,np.abs(mGradientAll[2]))
    for ix in range(N):
        for iy in range(N):
            for iz in range(N):
                if smoothMaterialBorders:
                    if abs(mGradient[ix,iy,iz]) != 0:
                        material[ix,iy,iz] = 0.5
                masse[ix,iy,iz] = material[ix,iy,iz]*(meB-meQ)+meQ
                massh[ix,iy,iz] = material[ix,iy,iz]*(mhB-mhQ)+mhQ
                Ve[ix,iy,iz] = material[ix,iy,iz]*dEc
                Vh[ix,iy,iz] = material[ix,iy,iz]*dEv

    addEfields()
    addPCs()
    addBZfield()
    modelCreated = True


# Fz: Electric field [V/m] along z-direction
# Fx: Electric field [V/m] along x-direction
def addEfields():
    global Ve, Vh
    q = 1           # q for V in [eV}
    Fz = float(EFields[0])
    if abs(Fz) > 0:
        for ix in range(len(Ve[0])):
            for iy in range(len(Ve[1])):
                for iz in range(len(Ve[2])):
                    x,y,z = cellPos([ix,iy,iz])
                    Ve[ix, iy, iz] = Ve[ix, iy, iz] + q*Fz*z
                    Vh[ix, iy, iz] = Vh[ix, iy, iz] - q*Fz*z
    Fx = float(EFields[1])
    if abs(Fx) > 0:
        for ix in range(len(Ve[0])):
            for iy in range(len(Ve[1])):
                for iz in range(len(Ve[2])):
                    x,y,z = cellPos([ix,iy,iz])
                    Ve[ix, iy, iz] = Ve[ix, iy, iz] + q*Fx*x
                    Vh[ix, iy, iz] = Vh[ix, iy, iz] - q*Fx*x


# Point charges
def addPCs():
    if len(PointCharges) > 0:
        q = cnst.e
        N = cellNum1D
        #c = 4*cnst.epsilon_GaAs*math.pi
        for iPC in range(len(PointCharges)):
            for ix in range(N):
                for iy in range(N):
                    for iz in range(N):
                        x,y,z = cellPos([ix,iy,iz])
                        xPC = PointCharges[iPC][0]
                        yPC = PointCharges[iPC][1]
                        zPC = PointCharges[iPC][2]
                        d = np.sqrt( pow((x-xPC),2) + pow((y-yPC),2) + pow((z-zPC),2) )
                        if d == 0: d = model.cellSize/4
                        Ve[ix, iy, iz] = Ve[ix, iy, iz] + q / (4*cnst.epsilon_GaAs*math.pi*d)
                        Vh[ix, iy, iz] = Vh[ix, iy, iz] - q / (4*cnst.epsilon_GaAs*math.pi*d)


# Bz: Magnetic field [T] along z-direction
def addBZfield():
    Bz = BZfield[0]
    Le = BZfield[1]
    Lh = BZfield[2]
    if (abs(Bz)>0) or (abs(Le)>0) or (abs(Lh)>0):
        q = cnst.e
        c0_e = (1/q)* cnst.hbar**2*Le**2 / (2*cnst.me_GaAs)     # in [eV]
        c0_h = (1/q)* cnst.hbar**2*Lh**2 / (2*cnst.mhh_GaAs)    # in [eV]
        c1_e =  (1/q)* q*cnst.hbar*Bz*Le / (2*cnst.me_GaAs)     # in [eV]
        c1_h = (1/q)* -q*cnst.hbar*Bz*Lh / (2*cnst.mhh_GaAs)    # in [eV]
        c2_e =  (1/q)* (q*Bz)**2 / (8*cnst.me_GaAs)             # in [eV]
        c2_h = (1/q)* (-q*Bz)**2 / (8*cnst.mhh_GaAs)            # in [eV]
        N = cellNum1D
        for ix in range(N):
            for iy in range(N):
                for iz in range(N):
                    x,y,z = cellPos([ix,iy,iz])
                    rxy2 =  (x**2 + y**2)
                    if rxy2 > 0:
                        Ve[ix, iy, iz] = Ve[ix, iy, iz] + c0_e/rxy2 + c1_e + c2_e*rxy2
                        Vh[ix, iy, iz] = Vh[ix, iy, iz] + c0_h/rxy2 + c1_h + c2_h*rxy2


# Dirichlet boundary condition in 1D: psi(0) = psi(N) = 0
# psi(0) = 0 means V(0) = large
def createBoundaries(V):
    N = cellNum1D
    for ix in range(len(V[0])):
            for iy in range(len(V[1])):
                for iz in range(len(V[2])):
                    if ix==0 or ix==N-1: V[ix, iy, iz] = 100
                    if iy==0 or iy==N-1: V[ix, iy, iz] = 100
                    if iz==0 or iz==N-1: V[ix, iy, iz] = 100
    return V


# shows planes through the calculated model potential
def ShowPotential(figModel, canvasModel):
    figModel.clf()
    ix=cellNum1D//2; iy=cellNum1D//2; iz=cellNum1D//2;
    figModel.add_subplot(131).imshow(Ve[ix, : , : ], interpolation='none', cmap=colMap)
    figModel.add_subplot(132).imshow(Ve[: , iy, : ], interpolation='none', cmap=colMap)
    figModel.add_subplot(133).imshow(Ve[: , : , iz], interpolation='none', cmap=colMap)
    canvasModel.draw()


# cell size = 2*r1 + 2*borderSize[nm]
def autoCellSize():
    global cellSize
    borderSize = 5e-9
    readInputVals()
    cellSize =  (2*QDpara[1] + 2*borderSize)/cellNum1D
    setInputVals()

