import numpy as np
import matplotlib.pyplot as plt
from numpy.typing import NDArray
from ambiance import Atmosphere
from scipy.optimize import fsolve
import math
import os
import subprocess
from textwrap import dedent
import time
from Meangen_control import run_meangen, run_stagen, run_multall
from Blade_Angles import howell_loading_criterion, diffusion_factor, calculate_incidence_deflection


def getMeanLineRadius(rTip: float, hubToTip: float): #area average
    rHub: float =  rTip * hubToTip
    # return np.sqrt( (rHub**2 + rTip**2)/ 2 )
    return 0.5 * (rHub + rTip)


def computeVelocityTrianglesWithAlpha1Known(psi: float, phi: float, alpha1: float = 0) -> tuple[float, float, float, float]: 
    #assuming the flow enters the rotor and leaves the stator axially (a1 = a3 = 0)
    #I believe this also assumes that the radius remains the same?
    
    beta2  : float = np.atan( (psi + phi * np.tan(alpha1) - 1) / phi)
    beta1  : float = np.atan( np.tan(alpha1) - 1 / phi)
    alpha2 : float = np.atan( np.tan(beta2) + 1 / phi)
    
    R : float = -psi / 2 - phi * np.tan(alpha1) + 1
    
    return alpha2, beta1, beta2, R

def computeVelocityTrianglesWithRKnown(psi: float, phi: float, R: float) -> tuple[float, float, float, float]: 
    alpha1 = math.atan(-1 * ( ((R + psi/2 - 1)/ phi)  ))
    beta2  : float = math.atan( (psi + phi * np.tan(alpha1) - 1) / phi)
    beta1  : float = math.atan( np.tan(alpha1) - 1 / phi)
    alpha2 : float = math.atan( np.tan(beta2) + 1 / phi)
    
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

    
def getPropertiesAfterRow(h0: float, deltaH0: float ,
                            S : float, eta: float, p0Inlet: float, h0Inlet: float,
                            v: float, Sinlet: float, entropyFraction: float,
                            cp: float = 1006,
                            k: float = 1.4, 
                            R: float = 287.05):
    
    h0out = h0 + deltaH0
    hOut  = h0out - 0.5 *  v**2
    
    T0out = (h0 + deltaH0) / cp
    Tout  = hOut/cp
    
    Sout  = S + (1-eta) * deltaH0 / T0out * entropyFraction
    
    pOut  = p0Inlet * (hOut/h0Inlet)**(k/(k-1)) * np.exp((Sinlet-Sout)/R)
    p0Out = pOut * (h0out/hOut)**(k/(k-1))
    
    rhoOut = pOut / R / Tout
    rho0Out = p0Out / R / T0out
    
    
    return h0out, hOut, T0out, Tout, Sout, p0Out, pOut, rho0Out, rhoOut
    

def getStagePropertiesAfterStage(omega: float, rMean: float, Caxial: float,
                       alpha1: float, alpha2:float,
                       psi: float, DOR: float, eta: float,
                       T0inlet: float, p0Inlet: float,
                            cp: float = 1006,
                            k: float = 1.4, 
                            R: float = 287.05):
    h0 = T0inlet * cp
    h0Inlet = h0


    U = omega * rMean

    deltaH0Total  = psi * U**2
    deltaH0Rotor  = deltaH0Total * DOR
    deltaH0Stator = deltaH0Total * (1-DOR)

    h0Rotor, hRotor, T0rotor, Trotor, Srotor, p0Rotor, pRotor, rho0Rotor, rhoRotor = getPropertiesAfterRow(
        h0 = h0,
        deltaH0 = deltaH0Total,
        S = 0,
        eta = eta,
        p0Inlet = p0Inlet,
        h0Inlet=h0Inlet,
        v = abs(Caxial/np.cos(alpha2)),
        Sinlet=0,
        entropyFraction=0.5
        
    )
    
    h03, h3, T03, T3, S3, p03, p3, rho03, rho3 = getPropertiesAfterRow(
        h0 = h0,
        deltaH0 = deltaH0Total,
        S = 0,
        eta = eta,
        p0Inlet = p0Inlet,
        h0Inlet=h0Inlet,
        v = abs(Caxial/np.cos(alpha1)),
        Sinlet=0,
        entropyFraction=1
        
    )
    
    deltaH0Isen = deltaH0Total/eta
    
    w = deltaH0Isen/eta
    etaTT = deltaH0Isen/w
    etaTS = (deltaH0Isen - 1/2 * Caxial**2) / w
    


    return T03, p03, rho03, p03/p0Inlet, T3, p3, rho3, T0rotor , p0Rotor, rho0Rotor, p0Rotor/p0Inlet, Trotor, pRotor, rhoRotor, etaTT, etaTS

def getStagnationInletProperties(h: float, M: float, k: float = 1.4, R: float = 287.05) -> tuple[float, float, float]:
    atmoConditions = Atmosphere(h)
    Tinf = atmoConditions.temperature
    Pinf = atmoConditions.pressure
    rhoInf = atmoConditions.density
    
    T0 = (1 + (k-1)/2 * M**2) * Tinf
    P0 = (T0/Tinf) ** (k/(k-1)) * Pinf
    rho0 = (T0/Tinf) ** (1/(k-1)) * rhoInf
    
    return T0, P0, rho0, Tinf, Pinf, rhoInf
    
def calculatePower(mdot: float, T0out: float, T0in: float, cp: float):
    return mdot * cp * (T0out - T0in)

def getStageEfficiencies(DeltaT: float, DeltaTis: float, c1: float, cp : float = 1006 ):
    #returns total-to-total and total-to-static
    tt = DeltaTis/DeltaT
    ts = (cp * DeltaTis) / (cp * DeltaT + 1/2 * c1**2) #!!!!!!!!!!!!!!!!!!!!!!!!!! I am really not sure this is actually correct
    
    return tt, ts


def getMassFlow(hubToTip: float, rTip: float, rho:float, AxialSpeed:float):
    rHub = rTip * hubToTip
    return np.pi * (rTip**2 - rHub**2) * AxialSpeed * rho

def getRelativeMachNumbers(hubToTip: float, rTip: float, AxialSpeed:float, omega: float, a: float):
    rHub = rTip * hubToTip   
    rMean = getMeanLineRadius(rTip,hubToTip)
    
    omegaRad = (omega * 2 * np.pi / 60)
    
    VrelHub = np.sqrt(AxialSpeed**2 + (omegaRad * rHub)**2) 
    VrelMean = np.sqrt(AxialSpeed**2 + (omegaRad * rMean)**2) 
    VrelTip = np.sqrt(AxialSpeed**2 + (omegaRad * rTip)**2) 
    
    
    return VrelHub/a, VrelMean/a, VrelTip/a

def getTipRadiusFromMassFlow(mdot: float, hubToTip: float, AxialSpeed: float, rho: float, blockage=1) -> float:
    A = mdot / AxialSpeed / rho / (1-blockage)
    return float(np.sqrt(A / (np.pi  * ( 1 - hubToTip**2))))

def getPitchOverCord(Vm: float, rho: float, p0: float, p: float, angleIn: float, angleOut: float, Zw: float = 0.5) -> float:
    Vt1 = np.tan(angleIn)  * Vm
    Vt2 = np.tan(angleOut) * Vm
    
    
    return Zw * (p0-p) / (Vm * rho * abs(Vt1 - Vt2))

def getRelativeTotalPressureOfRotor(T0:float,Vabs:float,Vrel:float,p:float,k:float=1.4,cp:float=1006):
    h0 = cp * T0
    hstatic = h0 - 1/2 * Vabs**2
    h0rel = hstatic + 1/2 * Vrel**2
    return p * (h0rel/hstatic)**(k/(k-1))

def getBladeNumberAndAxialCord(Z: float, Rmean: float, pitchOverCord: float, rTip: float, Tstatic: float, psi: float, phi: float, k: float = 1.4, R: float = 287.05 ):
    # Mtip = (rTip * omega * 2 * np.pi / 60) / np.sqrt(k * R * Tstatic)

    # def system(X):
    #     Cx = X[0]
    #     Z  = X[1]

    #     return [
    #         float((2 * np.pi * Rmean / Z) / Cx - pitchOverCord),
    #         float(Z * ( ( (Z * 1 / 2 / np.pi * np.sqrt((k-1)/k) * Cx / rTip * Mtip )**2 + 1  )**(k/(k-1)) - 1 ) - 2 * np.pi * k * rTip / Cx * psi * phi * Mtip**2)
            
    #     ]
        
    # X = fsolve(system,[10e-2,20])
    # Z = round(X[1])
    # Z = (X[1])
    Cx = (2 * np.pi * Rmean / Z ) / pitchOverCord 
    return Z, Cx

def getAeroForces(A: float, Z: float, p1: float, p2:float, mdot:float, Caxial:float, angle1: float, angle2: float):
    Fm =  A/Z * (p2-p1)
    Ft =  mdot/Z * Caxial * (np.sin(angle1) - np.sin(angle2))
    return np.sqrt(Fm**2 + Ft**2)
