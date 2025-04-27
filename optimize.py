from meanline import *
from loss_models import *
from scipy.optimize import fsolve
import matplotlib.pylab as plt


#-------------------------------------------------------------
Zrotor = 21
Zstator = 35

blockage = 0.00

bttTarget = 1.40
#-------------------------------------------------------------

alpha1 = 0
# phi = 0.563

altitude = 10e3 #[m]
Minf = 0.78 
M1abs = 0.6  #[-]
omega = 5000 #[rpm]

hubToTip = 0.30
mdot = 80 #[kg/s]
thickToCord = 0.06
twist = 0.55


Kv = 0.055 #https://elib.dlr.de/192376/1/ICAS2022_0497_paper.pdf page 20
rhoMat = 4650



#---------------------------------------------------------------------------------------------------------------------------


def getDeterminePsi(psi,bttTarget,optimize=True):
    T02, P02, rho02, Tinf, Pinf, rhoInf = getStagnationInletProperties(altitude,Minf)
    C1mag = M1abs * np.sqrt(1.4 * 287 * Tinf)
    T2, P2, rho2 = getStaticProperties(C1mag,T02,P02)
    rTip = getTipRadiusFromMassFlow(mdot,hubToTip,C1mag,rho2,blockage)
    rMean = getMeanLineRadius(rTip,hubToTip)
    
    phi = C1mag[0] / (rMean * omega * 2 * np.pi / 60)
    
    #alpha1, alpha2, beta1, beta2 = computeVelocityTrianglesWithRKnown(psi,phi,DOR)
    alpha2, beta1, beta2, DOR = computeVelocityTrianglesWithAlpha1Known(psi,phi,alpha1)
    
    etaIso = compute_eta(psi,phi,DOR,thickToCord)
    # etaIso = 1
    T03, P03, rho03, PRTotal , T3, P3, rho3, T0Rotor , P0Rotor, rho0Rotor, PRRotor, Trotor, Protor, rhoRotor, DeltaT, DeltaTis =  getStagePropertiesAfterStage(
        omega = omega * 2 * np.pi / 60,
        rMean=rMean,
        Caxial=C1mag[0],
        alpha1=alpha1,
        alpha2=alpha2,
        psi=psi,
        DOR=DOR,
        eta = etaIso,
        T0inlet=T02,
        p0Inlet=P02,
    )
    if optimize:
        return PRTotal - bttTarget
    return etaIso


#-------------------------------------------------------

btts = list(np.arange(1.4,1.8,0.01))
etas = []
for b in btts:
    psiL = fsolve(getDeterminePsi,0.5,args=(b))[0]
    eta = getDeterminePsi(psiL,b,False)
    etas.append(eta)
    
plt.plot(btts,etas)
    




# ---------------------------

psi = fsolve(getDeterminePsi,0.5,args=(bttTarget))[0]

T02, P02, rho02, Tinf, Pinf, rhoInf = getStagnationInletProperties(altitude,Minf)
C1mag = M1abs * np.sqrt(1.4 * 287 * Tinf)
T2, P2, rho2 = getStaticProperties(C1mag,T02,P02)
#mdot = getMassFlow(hubToTip,rTip,rho2,C1mag)
rTip = getTipRadiusFromMassFlow(mdot,hubToTip,C1mag,rho2,blockage)
rMean = getMeanLineRadius(rTip,hubToTip)


Umean = rMean * omega * 2 * np.pi / 60    
phi = C1mag[0] / (rMean * omega * 2 * np.pi / 60)

#alpha1, alpha2, beta1, beta2 = computeVelocityTrianglesWithRKnown(psi,phi,DOR)
alpha2, beta1, beta2, DOR = computeVelocityTrianglesWithAlpha1Known(psi,phi,alpha1)



MrelHub, MrelMean, MrelTip = getRelativeMachNumbers(hubToTip,rTip,C1mag,omega,np.sqrt(1.4 * 287 * T2)) 

etaIso = compute_eta(psi,phi,DOR,thickToCord)
# etaIso = 1

T03, P03, rho03, PRTotal , T3, P3, rho3, T0Rotor , P0Rotor, rho0Rotor, PRRotor, Trotor, Protor, rhoRotor, etaTtT, etaTtS =  getStagePropertiesAfterStage(
    omega = omega * 2 * np.pi / 60,
    rMean=rMean,
    Caxial=C1mag[0],
    alpha1=alpha1,
    alpha2=alpha2,
    psi=psi,
    DOR=DOR,
    eta = etaIso,
    T0inlet=T02,
    p0Inlet=P02,
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

Zrotor, CxRotor = getBladeNumberAndAxialCord(Zrotor,rMean,PitchOverCordRotor,rTip,T2,psi,phi)
Zstator, CxStator = getBladeNumberAndAxialCord(Zstator,rMean,PitchOverCordStator,rTip, Trotor,psi,phi)

# etaTtT, etaTtS = getStageEfficiencies(DeltaT, DeltaTis,C1mag)

solidityRotor  = float(1/PitchOverCordRotor ) #TODO Implement howell & diffusion factor to optimize for solidity
print("Solidity:::",solidityRotor)
solidityStator = float(1/PitchOverCordStator) #TODO Implement howell & diffusion factor to optimize for solidity

print("beta and alpha:::", beta1*360/(2*np.pi), beta2*360/(2*np.pi), alpha1*360/(2*np.pi), alpha2*360/(2*np.pi))
incidenceRotor, deflectionRotor = calculate_incidence_deflection(beta1*360/(2*np.pi), beta2*360/(2*np.pi), solidityRotor, thickToCord, "DCA") #TODO actually choose t/c and foil type
incidenceStator, deflectionStator = calculate_incidence_deflection(alpha2*360/(2*np.pi), alpha1*360/(2*np.pi), solidityStator, thickToCord, "DCA") #TODO actually choose t/c and foil type
print("Incidence, declection:::",incidenceRotor, deflectionRotor, incidenceStator, deflectionStator)

# print(incidenceRotor,deflectionRotor)

psiOpt = 0.185 * np.sqrt(4*phi**2+1)
psiMax = 0.32 + 0.2*phi



ArotorIn  = mdot/C1mag/rho2/(1-blockage)
AstatorIn   = mdot/C1mag/rhoRotor/(1-blockage) 
AstatorOut   = mdot/C1mag/rho3/(1-blockage) 

aeroForceRotor  = getAeroForces(0.5*(ArotorIn+AstatorIn),Zrotor,P2,Protor,mdot,C1mag,beta1,beta2)
aeroForceStator = getAeroForces(0.5*(AstatorIn+AstatorOut),Zstator,Protor,P3,mdot,C1mag,alpha2,alpha1)


rotorBladeMass = (ArotorIn+AstatorIn) * 0.5 * CxRotor * Kv * rhoMat / Zrotor
statorBladeMass = (AstatorIn+AstatorOut) * 0.5 * CxStator * Kv * rhoMat / Zstator


print( (omega * 2 * np.pi / 60))

centForceRotor = rotorBladeMass * rMean * (omega * 2 * np.pi / 60)**2
# centForceStator = statorBladeMass * rMean * (omega * 2 * np.pi / 60)**2


#Run Meangen & Stagen
dirPath = os.path.dirname(os.path.realpath(__file__))


if os.name == "posix":
    path_of_user = os.path.join(dirPath,"multallExecutables","Linux",)
else:
    path_of_user = os.path.join(dirPath,"multallExecutables","Windows",)



print(f"-------------------------------------------------")
print(f"rHub          [m]: {rTip*hubToTip:>10.2f}")
print(f"rMean         [m]: {rMean:>10.2f}")
print(f"rTip          [m]: {rTip:>10.2f}")
print(f"mdot       [kg/s]: {mdot:>10.2f}")
print(f"phi           [-]: {phi:>10.2f}")
print(f"psi           [-]: {psi:>10.2f}")
print(f"psiOpt        [-]: {psiOpt:>10.2f}")
print(f"psiMax        [-]: {psiMax:>10.2f}")
print(f"DOR           [-]: {DOR:>10.2f}")
print(f"alpha1      [deg]: {np.degrees(alpha1):>10.2f}")
print(f"beta1       [deg]: {np.degrees(beta1):>10.2f}")
print(f"alpha2      [deg]: {np.degrees(alpha2):>10.2f}")
print(f"beta2       [deg]: {np.degrees(beta2):>10.2f}")
print(f"Tinf          [K]: {Tinf[0]:>10.2f}")
print(f"Pinf         [Pa]: {Pinf[0]:>10.2f}")
print(f"rhoInf    [kg/m3]: {rhoInf[0]:>10.2f}")
print(f"Vm          [m/s]: {C1mag[0]:>10.2f}")
print(f"T02           [K]: {T02[0]:>10.2f}")
print(f"P02          [Pa]: {P02[0]:>10.2f}")
print(f"P02rel       [Pa]: {P02rel[0]:>10.2f}")
print(f"rho02     [kg/m3]: {rho02[0]:>10.2f}")
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
print(f"etaIso        [-]: {etaIso:>10.2f}")
print(f"eta TtT       [-]: {etaTtT:>10.2f}")
print(f"eta TtS       [-]: {etaTtS:>10.2f}")
print(f"1/sigma Rotor [-]: {PitchOverCordRotor[0]:>10.2f}")
print(f"1/sigma Stat. [-]: {PitchOverCordStator[0]:>10.2f}")
print(f"sigma Rotor   [-]: { PitchOverCordRotor[0]**-1:>10.2f}")
print(f"sigma Stat.   [-]: {PitchOverCordStator[0]**-1:>10.2f}")
print(f"Z Rotor       [-]: {Zrotor:>10.2f}")
print(f"Z Stator      [-]: {Zstator:>10.2f}")
print(f"Cx Rotor      [m]: { CxRotor[0]:>10.2f}")
print(f"Cx Stator     [m]: {CxStator[0]:>10.2f}")
print(f"Mass Rotor   [kg]: {rotorBladeMass[0]:>10.2f}")
print(f"Mass Stator  [kg]: {statorBladeMass[0]:>10.2f}")
print(f"Faero Rotor   [N]: {aeroForceRotor[0]:>10.2f}")
print(f"Faero Stator  [N]: {aeroForceStator[0]:>10.2f}")
print(f"Fcent Rotor  [kN]: {centForceRotor[0]/1e3:>10.2f}")
# print(f"Fcent Stator [kN]: {centForceStator[0]/1e3:>10.2f}")
print(f"-------------------------------------------------")

input("Give any input to continue with meangen execution")
run_meangen(path_of_user, round(P02[0]/1e5,3), round(T02[0],3), DOR, phi, psi, rMean,float(mdot),incidenceRotor,deflectionRotor,incidenceStator,deflectionStator,float(CxRotor),float(CxStator),float(etaIso),blockage, twist)
input("Give any input to continue with multall execution")
run_stagen(path_of_user)
run_multall(path_of_user)

