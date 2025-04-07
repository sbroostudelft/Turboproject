import numpy as np
import matplotlib.pyplot as plt
from numpy.typing import NDArray
from ambiance import Atmosphere
from scipy.optimize import fsolve

def getMeanLineRadius(rTip: float, hubToTip: float): #area average
    rHub: float =  rTip * hubToTip
    return np.sqrt( (rHub**2 + rTip**2)/ 2 )


def computeVelocityTrianglesWithAlpha1Known(psi: float, phi: float, alpha1: float = 0) -> tuple[float, float, float, float]: 
    #assuming the flow enters the rotor and leaves the stator axially (a1 = a3 = 0)
    #I believe this also assumes that the radius remains the same?
    
    beta2  : float = np.atan( (psi + phi * np.tan(alpha1) - 1) / phi)
    beta1  : float = np.atan( np.tan(alpha1) - 1 / phi)
    alpha2 : float = np.atan( np.tan(beta2) + 1 / phi)
    
    R : float = -psi / 2 - phi * np.tan(alpha1) + 1
    
    return alpha2, beta1, beta2, R

def computeVelocityTrianglesWithRKnown(psi: float, phi: float, R: float) -> tuple[float, float, float, float]: 
    alpha1 = np.atan(-1 * ( ((R + psi/2 - 1)/ phi)  ))
    beta2  : float = np.atan( (psi + phi * np.tan(alpha1) - 1) / phi)
    beta1  : float = np.atan( np.tan(alpha1) - 1 / phi)
    alpha2 : float = np.atan( np.tan(beta2) + 1 / phi)
    
    return alpha1, alpha2, beta1, beta2

def drawArrow(vector : NDArray, origin : NDArray, ax,color="black") -> None:
    ax.arrow(origin[0],origin[1],vector[0],vector[1],color=color)

def drawVelocityTriangles(alpha1: float, alpha2: float, beta1: float, beta2: float) -> None:
    C1 = 1 * np.array([np.cos(alpha1),np.sin(alpha1)])
    W1 = C1[0] * np.array([1,np.tan(beta1)])
    U1 = -W1 + C1
    
    
    C2 = C1[0] * np.array([1,np.tan(alpha2)]) # axial the same
    W2 = C1[0] * np.array([1,np.tan(beta2)])
    
    
    fig,ax = plt.subplots(1,1)
    
    ax.set_ylim((-3,3))
    ax.set_xlim((0,3))
    
    drawArrow(vector=C1, origin=np.array([0,0]), ax = ax, color="red")
    drawArrow(vector=U1, origin=W1, ax = ax,color="green")
    drawArrow(vector=W1, origin=np.array([0,0]), ax = ax,color="blue")
    
    drawArrow(vector=C2, origin=W1, ax = ax,color="red")
    drawArrow(vector=U1, origin=W1+W2, ax = ax,color="green")
    drawArrow(vector=W2, origin=W1, ax = ax,color="blue")
    
    drawArrow(vector=C1, origin=W1 + C2, color = "red",ax = ax)
    
    plt.show()


def getStaticProperties(C: float, T0: float, P0: float, k : float = 1.4, R : float = 287.05):
    def solve(T):
        return T0 / T - (1 + (k-1)/2 * (C / np.sqrt(k*R*T))**2 )
    
    T = fsolve(solve,T0)
    P = (T/T0) ** (k / (k-1)) * P0
    rho = P / R / T
    
    return T,P,rho
    
def getPropertiesAfterStage(T0in: float, P0in: float, C1mag: float, omega: float, rMean: float, alpha1: float, alpha2: float, reaction: float ,cp: float, k: float = 1.4, etaIso: float = 1, R: float = 287.05) -> tuple[float,float,float,float]:
    U = omega * rMean #assuming constant radius
    Ca = C1mag * np.cos(alpha1)
    
    C2mag = Ca / np.cos(alpha2)
    
    DeltaTis = U * Ca / cp * (np.tan(alpha2) - np.tan(alpha1)) 
    DeltaT = DeltaTis / etaIso
    
    T03is = T0in + DeltaTis
    T03   = T0in + DeltaT
    
    T0RotorIs = T0in + DeltaTis * reaction
    T0Rotor   = T0in + DeltaT * reaction
    
    P03 = (T03is/T0in) ** (k / (k-1)) * P0in
    rho03 = P03 / (R * T03)
    
    
    P0Rotor = (T0RotorIs/T0in) ** (k / (k-1)) * P0in
    rho0Rotor = P0Rotor / (R * T0Rotor)
    
    
    T3, P3, rho3 = getStaticProperties(C1mag,T03,P03)
    Trotor, Protor, rhoRotor = getStaticProperties(C2mag,T0Rotor,P0Rotor)
    
    return T03, P03, rho03, P03/P0in , T3, P3, rho3, T0Rotor , P0Rotor, rho0Rotor, P0Rotor/P0in, Trotor, Protor, rhoRotor


def getStagnationInletProperties(h: float, M: float, k: float = 1.4, R: float = 287.05) -> tuple[float, float, float]:
    atmoConditions = Atmosphere(h)
    T = atmoConditions.temperature
    P = atmoConditions.pressure
    rho = atmoConditions.density
    
    T0 = (1 + (k-1)/2 * M**2) * T
    P0 = (T0/T) ** (k/(k-1)) * P
    rho0 = (T0/T) ** (1/(k-1)) * rho
    
    return T0, P0, rho0
    
def calculatePower(mdot: float, T0out: float, T0in: float, cp: float):
    return mdot * cp * (T0out - T0in)

def getStageEfficiencies(P0in: float, P03: float, P3: float):
    #returns total-to-total and total-to-static
    return P03/P0in , P3/P0in

if __name__ == "__main__":
    alpha1, alpha2, beta1, beta2 = computeVelocityTrianglesWithRKnown(0.5,0.5,0.5)
    
    T0, P0, rho0 = getStagnationInletProperties(10e3,0.78)
    
    T03, P03, rho03, PRTotal , T3, P3, rho3, T0Rotor , P0Rotor, rho0Rotor, PRRotor, Trotor, Protor, rhoRotor =  getPropertiesAfterStage(
        T0in = T0,
        P0in = P0,
        C1mag = 0.6 * 295, # this is not correct exactly just for reference
        omega = 5000 * 2 * np.pi / 60,
        rMean = 0.7,
        alpha1 = alpha1,
        alpha2 = alpha2,
        reaction= 0.5,
        cp = 1006,
    )
    
    etaTtT, etaTtS = getStageEfficiencies(P0,P03,P3)
    

    print(f"-------------------------------------------------")
    print(f"alpha1      [deg]: {np.degrees(alpha1):>10.2f}")
    print(f"beta1       [deg]: {np.degrees(beta1):>10.2f}")
    print(f"alpha2      [deg]: {np.degrees(alpha2):>10.2f}")
    print(f"beta2       [deg]: {np.degrees(beta2):>10.2f}")
    print(f"T0            [K]: {T0[0]:>10.2f}")
    print(f"P0           [Pa]: {P0[0]:>10.2f}")
    print(f"rho0      [kg/m3]: {rho0[0]:>10.2f}")
    print(f"T0Rotor       [K]: {T0Rotor[0]:>10.2f}")
    print(f"P0Rotor      [Pa]: {P0Rotor[0]:>10.2f}")
    print(f"rho0Rotor [kg/m3]: {rho0Rotor[0]:>10.2f}")
    print(f"TRotor        [K]: {Trotor[0]:>10.2f}")
    print(f"PRotor       [Pa]: {Protor[0]:>10.2f}")
    print(f"rhoRotor  [kg/m3]: {rhoRotor[0]:>10.2f}")
    print(f"PR Rotor      [-]: {PRRotor[0]:>10.2f}")
    print(f"T03           [K]: {T03[0]:>10.2f}")
    print(f"P03          [Pa]: {P03[0]:>10.2f}")
    print(f"rho03     [kg/m3]: {rho03[0]:>10.2f}")
    print(f"T3            [K]: {T3[0]:>10.2f}")
    print(f"P3           [Pa]: {P3[0]:>10.2f}")
    print(f"rho3      [kg/m3]: {rho3[0]:>10.2f}")
    print(f"PR Total      [-]: {PRTotal[0]:>10.2f}")
    print(f"eta TtT       [-]: {etaTtT[0]:>10.2f}")
    print(f"eta TtS       [-]: {etaTtS[0]:>10.2f}")
    print(f"-------------------------------------------------")
    #drawVelocityTriangles(alpha1,alpha2,beta1,beta2)