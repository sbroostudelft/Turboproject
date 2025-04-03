import numpy as np
import matplotlib.pyplot as plt
from numpy.typing import NDArray
from ambiance import Atmosphere

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

def  computeVelocityTrianglesWithRKnown(psi: float, phi: float, R: float) -> tuple[float, float, float, float]: 
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


def getPropertiesAfterStage(Tin: float, Pin: float, rhoIn: float, C1mag: float, omega: float, rMean: float, alpha1: float, alpha2: float, cp: float, k: float = 1.4, etaIso: float = 1) -> tuple[float,float,float,float]:
    U = omega * rMean #assuming constant radius
    Ca = C1mag * np.cos(alpha1)
    
    DeltaTis = U * Ca / cp * (np.tan(alpha2) - np.tan(alpha1)) 
    DeltaT = DeltaTis / etaIso
    
    T03is = Tin + DeltaTis
    T03   = Tin + DeltaT
    
    P03 = (T03is/T03) ** (k / (k-1)) * Pin
    rho03 = (T03is/T03) ** (1 / (k-1)) * rhoIn
    
    
    return T03, P03, rho03, P03/Pin


def getStagnationInletProperties(h: float, M: float, k: float = 1.4) -> tuple[float, float, float]:
    atmoConditions = Atmosphere(h)
    T = atmoConditions.temperature
    P = atmoConditions.pressure
    rho = atmoConditions.density
    
    T0 = (1 + (k-1)/2 * M**2) * T
    P0 = (T0/T) ** (k/(k-1)) * P
    rho0 = (T0/T) ** (1/(k-1)) * rho
    
    return T0, P0, rho0
    

if __name__ == "__main__":
    alpha1, alpha2, beta1, beta2 = computeVelocityTrianglesWithRKnown(0.3954,0.5633,0.5)
    print(np.degrees(alpha1),np.degrees(alpha2),np.degrees(beta1),np.degrees(beta2))
    
    
    drawVelocityTriangles(alpha1,alpha2,beta1,beta2)