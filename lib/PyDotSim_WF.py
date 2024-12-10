"""

Wave function analysis for PyDotSim

by Christian Heyn

"""

import timeit
import math
import numpy as np

import lib.PyDotSim_OS as osprocs
import lib.PyDotSim_Constants as cnst
import lib.PyDotSim_Model as model


psie =[]; psih = []
deh =       0
overlap =   0
tau =       0
Ceh =       0
Cee =       0
Chh =       0
E_PL =      0       # dot emission energy in eV
rxe =       0
rye =       0
rze =       0
rxh =       0
ryh =       0
rzh =       0



# Normalize WF for integrated WF^2 = 1
def normalizeWF(psi):
    int_WF = 0; sum_WF = 0
    for ix in range(len(psi[0])):
        for iy in range(len(psi[1])):
            for iz in range(len(psi[2])):
                sum_WF += psi[ix,iy,iz]
                int_WF = int_WF + psi[ix,iy,iz]**2*model.cellVol
    WF = psi/np.sqrt(int_WF)
    if sum_WF < 1:
        WF = -WF        #  the solver produces randomly inverted psi
    return WF


def WFcenterOfMass(WF):
    sumx = 0; sumy = 0; sumz = 0; sumWF2 = 0
    for ix in range(len(WF[0])):
        for iy in range(len(WF[1])):
            for iz in range(len(WF[2])):
                WF2 = WF[ix,iy,iz]**2
                x,y,z = model.cellPos([ix,iy,iz])
                sumWF2 += WF2
                sumx += x*WF2; sumy += y*WF2; sumz += z*WF2;
    rx = sumx/sumWF2; ry = sumy/sumWF2; rz = sumz/sumWF2;
    return rx, ry, rz


#  Overlap integral
def overlapWF(WFe, WFh):
    integr = 0
    cellVol = pow(model.cellSize,3)
    for ix in range(len(WFe[0])):
        for iy in range(len(WFe[1])):
            for iz in range(len(WFe[2])):
               integr = integr + WFe[ix,iy,iz]*WFh[ix,iy,iz]*cellVol
    return integr**2


# QD radiative lifetime
# according Dalgarno, PRB 77 (2008)
def WFlifetime(overlapp):
    tau = 0;
    if overlapp > 0:
        Eg = cnst.GaAs_bandgap(model.T)
        tau = cnst.QDlifetime(overlapp,E_PL)
    tau_ns = tau*1e9;
    return tau_ns


# Transforms 3D WF into linear list
# considers only elements larger than threshold
def flattenWF(WF, threshold):
    WFmax = np.max(WF)
    threshold = threshold*WFmax
    WFflat = []
    for ix in range(len(WF[0])):
        for iy in range(len(WF[1])):
            for iz in range(len(WF[2])):
                WFval = WF[ix,iy,iz]
                if WFval > threshold:
                    WFval2 = pow(WFval, 2)
                    x,y,z = model.cellPos([ix,iy,iz])
                    WFflat.append([x,y,z,WFval,WFval2,ix,iy,iz])
    return WFflat


# Coulomb integral
def Coulomb(WFaflat, WFbflat):
    dV = model.cellSize*model.cellSize*model.cellSize
    int_ab = 0
    for ib in range(len(WFbflat)):
        xb = WFbflat[ib][0]; yb = WFbflat[ib][1]; zb = WFbflat[ib][2];
        WFb2 = WFbflat[ib][4]
        for ia in range(len(WFaflat)):
            xa = WFaflat[ia][0]; ya = WFaflat[ia][1]; za = WFaflat[ia][2];
            WFa2 = WFaflat[ia][4]
            dab  = np.sqrt( pow((xa-xb), 2) + pow((ya-yb), 2) + pow((za-zb), 2) )
            if dab == 0:
                dab = model.cellSize/4
            int_ab = int_ab + WFa2*WFb2/dab
    C_ab = pow(cnst.e,2)/(4*math.pi*cnst.epsilon_GaAs)*int_ab*pow(dV,2)*cnst.JeV
    return C_ab


# Calculation of the additional potential at a single mesh cell with position (x1,y1,z1) induced by the WF of another particle
# For this WF as sum of point charges at positions (x2,y2,z2) with respctive potential V = 1/(4 pi epsilon) q / r
# d12 = r = distance to point charge
# q = WFsquare: charge at position (x2,y2,z2)
# V = integrate_(x2,y2,z2) [ (q/d12) / (4 pi epsilon_GaAs) ]
def calcCellWFpot(x1, y1, z1, WFflat):
    int_q = 0
    for i in range(len(WFflat)):
        x2 = WFflat[i][0]; y2 = WFflat[i][1]; z2 = WFflat[i][2];
        d12  = np.sqrt( pow((x2-x1), 2) + pow((y2-y1), 2) + pow((z2-z1), 2) );
        WFsquare = WFflat[i][4]     # gives the local charge
        if d12 == 0:
            d12 = model.cellSize/4
        int_q = int_q+WFsquare/d12
    int_q = int_q*model.cellVol
    V = int_q*cnst.e/(4*math.pi*cnst.epsilon_GaAs)
    return V


# Adds a WF as a potential V for self-consistent many-particle
# X  (1e,1h) :  Ve - 1x V_WFh,              Vh - 1x V_WFe
# XX (2e,2h):   Ve - 2x V_WFh + 1x V_WFe,   Vh - 2x V_WFe + 1x V_WFh
# X+ (1e,2h):   Ve - 2x V_WFh,              Vh - 1x V_WFe + 1x V_WFh
# X- (2e,1h):   Ve - 1x V_WFh + 1x V_WFe,   Vh - 2x V_WFe
def addWF2pot(loop):
    WFe = osprocs.loadWF('e',loop)
    WFh = osprocs.loadWF('h',loop)
    WFeflat = flattenWF(WFe, model.WFflat_threshold)
    WFhflat = flattenWF(WFh, model.WFflat_threshold)
    particle = model.MultiParticleStr[model.MPpara[0]]
    # potential for e
    for i in range(len(WFeflat)):
        x = WFeflat[i][0]; y = WFeflat[i][1]; z = WFeflat[i][2];
        V_WFh = calcCellWFpot(x, y, z, WFhflat)
        if particle in ['XX', 'X-']:
            V_WFe = calcCellWFpot(x, y, z, WFeflat)
        ix = WFeflat[i][5]; iy = WFeflat[i][6]; iz = WFeflat[i][7]
        if particle == 'X':
            model.Ve[ix, iy, iz] = model.Ve[ix, iy, iz] - V_WFh
        if particle == 'XX':
            model.Ve[ix, iy, iz] = model.Ve[ix, iy, iz] - 2*V_WFh + V_WFe
        if particle == 'X+':
            model.Ve[ix, iy, iz] = model.Ve[ix, iy, iz] - 2*V_WFh
        if particle == 'X-':
            model.Ve[ix, iy, iz] = model.Ve[ix, iy, iz] - V_WFh + V_WFe
    # potential for h
    for i in range(len(WFhflat)):
        x = WFhflat[i][0]; y = WFhflat[i][1]; z = WFhflat[i][2];
        V_WFe = calcCellWFpot(x, y, z, WFeflat)
        if particle in ['XX', 'X+']:
            V_WFh = calcCellWFpot(x, y, z, WFhflat)
        if particle == 'X':
            model.Vh[ix, iy, iz] = model.Vh[ix, iy, iz] - V_WFe
        if particle == 'XX':
            model.Vh[ix, iy, iz] = model.Vh[ix, iy, iz] - 2*V_WFe + V_WFh
        if particle == 'X+':
            model.Vh[ix, iy, iz] = model.Vh[ix, iy, iz] - V_WFe + V_WFh
        if particle == 'X-':
            model.Vh[ix, iy, iz] = model.Vh[ix, iy, iz] - 2*V_WFe


def analyzeWF():
    global deh, overlap, tau, Ceh, Cee, Chh, rxe, rye, rze, rxh, ryh, rzh
    withMultiexcitons = False
    particle = model.MultiParticleStr[model.MPpara[0]]
    if particle in ['XX', 'X+', 'X-']: withMultiexcitons = True
    start = timeit.default_timer()
    WFe = normalizeWF(psie[0])
    WFh = normalizeWF(psih[0])
    rxe, rye, rze = WFcenterOfMass(WFe)
    rxh, ryh, rzh = WFcenterOfMass(WFh)
    deh = np.sqrt( np.square(rxe-rxh) + np.square(rye-ryh) + np.square(rze-rzh))
    overlap = overlapWF(WFe, WFh)
    tau = WFlifetime(overlap)
    WFeflat = flattenWF(WFe, model.WFflat_threshold)
    WFhflat = flattenWF(WFh, model.WFflat_threshold)
    Ceh = Coulomb(WFeflat, WFhflat)
    if withMultiexcitons:
        Cee = Coulomb(WFeflat, WFeflat)
        Chh = Coulomb(WFhflat, WFhflat)
    stop = timeit.default_timer()
    runTime = stop - start
    return runTime

