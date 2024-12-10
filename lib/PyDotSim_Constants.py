"""

Constants for SimQDStates

by Christian Heyn

"""

import math


h = 6.62606876e-34;             # Planck's constant [J*s]
hbar = h/(2*math.pi);           # Planck's constant reduced [J*s]
epsilon0 = 8.854187817e-12;     # Vacuum permittivity [F/m] = [A*s/(V*m)]
me0 = 9.10938188e-31;           # Electron mass [kg]
c = 299792458;                  # speed of light in m/s
e = 1.60217662e-19;             # elementary charge in coulombs
JeV = 1/1.60217646e-19;         # Joule in eV

# GaAs QDs related constants
me = 0.067
me_GaAs = me*me0                                        # GaAs effective electron mass [kg]
mh = 0.51
mhh_GaAs = mh*me0                                       # GaAs effective heavy hole mass [kg]
epsilon_GaAs = 13.1*epsilon0;                           # GaAs permittivity
m_Ex = 1/(1/me_GaAs + 1/mhh_GaAs);                      # Exciton effective mass
r_Ex_GaAs = 4*math.pi*epsilon_GaAs*hbar**2/(m_Ex*e**2); # GaAs exciton Bohr radius in [m]
E_Ry_GaAs = m_Ex*e**4/(2*(2*epsilon_GaAs*h)**2);        # Rydberg constant in GaAs bulk in [J]
E_p = 28.8;                                             # Kane energy in eV
n = 3.6;                                                # refractive index


def AlGaAs_effMass(x):
    me = 0.067+0.057*x
    mhh = 0.51+0.25*x
    return round(me,4), round(mhh,4)

def GaAs_bandgap(T): # Varshni
    # Vurgaftman, Thurmond
    EGaAs0 = 1.519; a = 5.405e-4; b = 204;
    EGaAs = EGaAs0 - a*T*T/(T+b);
    return EGaAs

def AlAs_bandgap(T): # Varshni
    # Vurgaftman
    EGaAs0 = 3.099; a = 8.85e-4; b = 530;
    EGaAs = EGaAs0 - a*T*T/(T+b);
    return EGaAs

def AlGaAs_bandgap(T, x): # Varshni
    # Vurgaftman
    C = -0.127 + 1.310*x
    EGaAs = GaAs_bandgap(T)
    EAlAs = AlAs_bandgap(T)
    EAlGaAs = EGaAs + (EAlAs-EGaAs-C)*x+C*x*x
    return EAlGaAs

# GaAs - AlGaAs valence band edge discontinuity
# from Adachi, Properties of semiconductor alloys, 2009: dEc/dEv = 63/37 = f
# dEc+dEv = dE = EAlGaAs-EGaAs
# dEc = dEv*f = (dE-dEc)*f = dE*f-dEc*f
# dEc+dE*f = dEc*(1+f) = dE*f
# dEc = dE*f/(1+f)
def bandgapDiscontinuity(T, x):
    EGaAs = GaAs_bandgap(T)
    EAlGaAs = AlGaAs_bandgap(T, x)
    dE = EAlGaAs-EGaAs
    f = 63/37
    dEc = dE*f/(1+f)
    dEv = dE-dEc
    return dEc, dEv

def QDlifetime(overlapp, E_PL):
    tau = (1/2) * 3 * h**2 * c**3 * epsilon0 * me0 / (math.pi * n * e**2 * E_p/JeV * E_PL/JeV * overlapp)
    return tau
