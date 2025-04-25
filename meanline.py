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
    
def getPropertiesAfterStage(T0in: float, P0in: float, C1mag: float, omega: float, rMean: float, alpha1: float, alpha2: float, reaction: float ,cp: float = 1006, k: float = 1.4, etaIso: float = 1, R: float = 287.05) -> tuple[float,float,float,float]:
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
    
    return T03, P03, rho03, P03/P0in , T3, P3, rho3, T0Rotor , P0Rotor, rho0Rotor, P0Rotor/P0in, Trotor, Protor, rhoRotor, DeltaT, DeltaTis


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

def getTipRadiusFromMassFlow(mdot: float, hubToTip: float, AxialSpeed: float, rho: float) -> float:
    return float(np.sqrt(mdot / (np.pi  * ( 1 - hubToTip**2) * AxialSpeed * rho)))

def getPitchOverCord(Vm: float, rho: float, p0: float, p: float, angleIn: float, angleOut: float, Zw: float = 0.5) -> float:
    Vt1 = np.tan(angleIn)  * Vm
    Vt2 = np.tan(angleOut) * Vm
    
    
    return Zw * (p0-p) / (Vm * rho * abs(Vt1 - Vt2))

def getRelativeTotalPressureOfRotor(T0:float,Vabs:float,Vrel:float,p:float,k:float=1.4,cp:float=1006):
    h0 = cp * T0
    hstatic = h0 - 1/2 * Vabs**2
    h0rel = hstatic + 1/2 * Vrel**2
    return p * (h0rel/hstatic)**(k/(k-1))

def getBladeNumberAndAxialCord(Rmean: float, pitchOverCord: float, rTip: float, Tstatic: float, psi: float, phi: float, k: float = 1.4, R: float = 287.05 ):
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
    Z = 22
    Cx = (2 * np.pi * Rmean / Z ) / pitchOverCord 
    return Z, Cx

if __name__ == "__main__":
    #user input
    
    alpha1 = 0
    psi = 0.395
    # phi = 0.563
    
    altitude = 10e3 #[m]
    Minf = 0.78 
    
    M1abs = 0.6  #[-]
    omega = 5000 #[rpm]
    
    hubToTip = 0.3
    #rTip = 0.581 #[m]
    mdot = 80 #[kg/s]
    
    etaIso = 1
    
    
    # ---------------------------
    
    
    
    T02, P02, rho02, Tinf, Pinf, rhoInf = getStagnationInletProperties(altitude,Minf)
    C1mag = M1abs * np.sqrt(1.4 * 287 * Tinf)
    T2, P2, rho2 = getStaticProperties(C1mag,T02,P02)
    #mdot = getMassFlow(hubToTip,rTip,rho2,C1mag)
    rTip = getTipRadiusFromMassFlow(mdot,hubToTip,C1mag,rho2)
    rMean = getMeanLineRadius(rTip,hubToTip)
    
    
    Umean = rMean * omega * 2 * np.pi / 60    
    phi = C1mag[0] / (rMean * omega * 2 * np.pi / 60)
    
    #alpha1, alpha2, beta1, beta2 = computeVelocityTrianglesWithRKnown(psi,phi,DOR)
    alpha2, beta1, beta2, DOR = computeVelocityTrianglesWithAlpha1Known(psi,phi,alpha1)
    
    
    
    MrelHub, MrelMean, MrelTip = getRelativeMachNumbers(hubToTip,rTip,C1mag,omega,np.sqrt(1.4 * 287 * T2)) 
    
    
    T03, P03, rho03, PRTotal , T3, P3, rho3, T0Rotor , P0Rotor, rho0Rotor, PRRotor, Trotor, Protor, rhoRotor, DeltaT, DeltaTis =  getPropertiesAfterStage(
        T0in = T02,
        P0in = P02,
        C1mag = C1mag,
        omega = omega * 2 * np.pi / 60,
        rMean = rMean,
        alpha1 = alpha1,
        alpha2 = alpha2,
        reaction= DOR,
        cp = 1006,
        etaIso = etaIso
    )
    
   
    P02rel = getRelativeTotalPressureOfRotor(T02,C1mag,np.sqrt(C1mag**2 + Umean**2),P2)
    
    PitchOverCordRotor = getPitchOverCord(
        Vm = C1mag[0],
        rho = 0.5 * (rho2 + rhoRotor) ,
        p0 = P02rel,
        p  = P2,
        angleIn=beta1,
        angleOut=beta2
    )
    
    PitchOverCordStator = getPitchOverCord(
        Vm = C1mag[0],
        rho = 0.5 * (rho3 + rhoRotor) ,
        p0 = P0Rotor,
        p  = Protor,
        angleIn=alpha2,
        angleOut=alpha1
    )
    
    Zrotor, CxRotor = getBladeNumberAndAxialCord(rMean,PitchOverCordRotor,rTip,T2,psi,phi)
    Zstator, CxStator = getBladeNumberAndAxialCord(rMean,PitchOverCordStator,rTip, Trotor,psi,phi)
    
    etaTtT, etaTtS = getStageEfficiencies(DeltaT, DeltaTis,C1mag)
    
    solidity = 1.5 #TODO Implement howell & diffusion factor to optimize for solidity

    incidence, deflection = calculate_incidence_deflection(beta1, beta2, solidity, 0.1, "DCA") #TODO actually choose t/c and foil type

    

    #Run Meangen & Stagen
    dirPath = os.path.dirname(os.path.realpath(__file__))

    
    if os.name == "posix":
        path_of_user = os.path.join(dirPath,"multallExecutables","Linux",)
    else:
        path_of_user = os.path.join(dirPath,"multallExecutables","Windows",)
    
    run_meangen(path_of_user, round(P02[0]/1e5,3), round(T02[0],3), DOR, phi, psi, rMean,float(mdot),incidence,deflection,float(CxRotor),float(CxStator))
    # run_stagen(path_of_user)
    #run_multall(path_of_user)
    

    print(f"-------------------------------------------------")
    print(f"rHub          [m]: {rTip*hubToTip:>10.2f}")
    print(f"rMean         [m]: {rMean:>10.2f}")
    print(f"rTip          [m]: {rTip:>10.2f}")
    print(f"mdot       [kg/s]: {mdot:>10.2f}")
    print(f"phi           [-]: {phi:>10.2f}")
    print(f"psi           [-]: {psi:>10.2f}")
    print(f"DOR           [-]: {DOR:>10.2f}")
    print(f"alpha1      [deg]: {np.degrees(alpha1):>10.2f}")
    print(f"beta1       [deg]: {np.degrees(beta1):>10.2f}")
    print(f"alpha2      [deg]: {np.degrees(alpha2):>10.2f}")
    print(f"beta2       [deg]: {np.degrees(beta2):>10.2f}")
    print(f"Tinf          [K]: {Tinf[0]:>10.2f}")
    print(f"Pinf         [Pa]: {Pinf[0]:>10.2f}")
    print(f"rhoInf    [kg/m3]: {rhoInf[0]:>10.2f}")
    print(f"Vm          [m/s]: {C1mag[0]:>10.2f}")
    print(f"T02            [K]: {T02[0]:>10.2f}")
    print(f"P02           [Pa]: {P02[0]:>10.2f}")
    print(f"P02rel        [Pa]: {P02rel[0]:>10.2f}")
    print(f"rho02      [kg/m3]: {rho02[0]:>10.2f}")
    print(f"T2            [K]: {T2[0]:>10.2f}")
    print(f"P2           [Pa]: {P2[0]:>10.2f}")
    print(f"MrelHub       [-]: {MrelHub[0]:>10.2f}")
    print(f"MrelMean      [-]: {MrelMean[0]:>10.2f}")
    print(f"MrelTip       [-]: {MrelTip[0]:>10.2f}")
    print(f"rho2      [kg/m3]: {rho2[0]:>10.2f}")
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
    print(f"1/sigma Rotor [-]: {PitchOverCordRotor[0]:>10.2f}")
    print(f"1/sigma Stat. [-]: {PitchOverCordStator[0]:>10.2f}")
    print(f"sigma Rotor   [-]: { PitchOverCordRotor[0]**-1:>10.2f}")
    print(f"sigma Stat.   [-]: {PitchOverCordStator[0]**-1:>10.2f}")
    print(f"Z Rotor       [-]: {Zrotor:>10.2f}")
    print(f"Z Stator      [-]: {Zstator:>10.2f}")
    print(f"Cx Rotor      [m]: { CxRotor[0]:>10.2f}")
    print(f"Cx Stator     [m]: {CxStator[0]:>10.2f}")
    print(f"-------------------------------------------------")
    #drawVelocityTriangles(alpha1,alpha2,beta1,beta2)