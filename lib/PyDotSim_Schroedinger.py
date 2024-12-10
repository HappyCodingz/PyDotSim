"""

3D Schroedinger solver for PyDotSim

Conversion of my modified Version of Schroedinger 3D (Octave),
https://github.com/LaurentNevou/Q_Schrodinger3D_demo
to Python

by Christian Heyn

"""

import numpy as np
import numpy.matlib
import math
import scipy.sparse
import scipy.sparse.linalg

import lib.PyDotSim_Constants as cnst
import lib.PyDotSim_Model as model


# Initial mass m and potential V matrix have shape (N,N,N)
# Hamiltonian matrix creation using scipy.sparse.lil_matrix
# Transformation of Hamiltonian matrix into csr for computation
def calcSchroedingerFEM(eh):
    a = model.cellSize
    N = model.cellNum1D; N2 = N*N; N3 = N*N*N

    # mass m and potential V matrix with shape (N,N,N)
    if eh == 'e':
        masseh = model.masse
        Veh = np.copy(model.Ve)*cnst.e
    if eh == 'h':
        masseh = model.massh
        Veh = np.copy(model.Vh)*cnst.e
    # transformation of m and V matrix into linear shape (N3)
    mass = masseh.reshape(1, N3)[0]
    # !!! With Boundaries: solver finds no solution
    # Veh = model.createBoundaries(Veh)
    V = Veh.reshape(1, N3)[0]

    # Hamiltonian H = k*Laplace operator + V
    # k = (-cnst.hbar**2/(2*cnst.me0*Mass)) * (1/a**2)
    k0 = (-cnst.hbar**2/(2*cnst.me0)) * (1/a**2)
    H = scipy.sparse.lil_matrix((N3, N3))
    for i in range(N3):
        f = 0
        j = i - N2
        if j >= 0: H[i,j]    = k0/mass[i]; f += 1
        j = i - N
        if j >= 0: H[i,j]    = k0/mass[i]; f += 1
        j = i - 1
        if j >= 0: H[i,j]    = k0/mass[i]; f += 1
        j = i + 1
        if j <= N3-1: H[i,j] = k0/mass[i]; f += 1
        j = i + N
        if j <= N3-1: H[i,j] = k0/mass[i]; f += 1
        j = i + N2
        if j <= N3-1: H[i,j] = k0/mass[i]; f += 1
        j = i
        H[i,j]               = -f*k0/mass[i] + V[i]

    # transformation of Hamiltonian into csr matrix for faster computation
    H = H.tocsr()

    # Compute Schroedinger equation via Hamiltonian matrix
    # constant seed for reproducible results via v0
    np.random.seed(0)
    v0 = np.random.rand(min(H.shape))
    Energy, PSI = scipy.sparse.linalg.eigsh(H, k=model.statesNum, which='SM', v0=v0 )
    E = Energy/cnst.e
    psi = []
    for i in range(model.statesNum):
        psi_temp = PSI[0:, i]
        #psi_temp = PSI[0:,i:i+1]
        psi_temp = psi_temp.reshape(N,N,N)
        psi_temp[0,0,0] = psi_temp[1,1,1]
        psi.append(psi_temp)
    psi = np.array(psi)
    return E, psi


#  correction of randomly inverted psi calculated by the solver
def correctInvertedPSI(psi):
    psi0 = psi[0]
    sum_WF = 0
    for ix in range(len(psi0[0])):
        for iy in range(len(psi0[1])):
            for iz in range(len(psi0[2])):
                sum_WF += psi0[ix,iy,iz]
    if sum_WF < 1:
        psi = -psi
    return psi
