import numpy as np
from duty_coefficients import calculate_duty_coefficients

f = 80 #kg/s
rho_ISA = 0.413 # kg/m^3
T_ISA = 223.2 # K
R = 287.05 # J/kgK
P_ISA = rho_ISA * R * T_ISA # Pa
M_inf = 0.78
M_ax = 0.6
M_w_1 = 1.4
gamma = 1.4
beta_tt = 1.6
eta_initial = 0.9
hub_to_tip = 0.3
RPM = 5000 #rpm

"Calculation of the duty coefficients"
psi, phi, DOR = calculate_duty_coefficients(np.sqrt(M_w_1 ** 2 - M_ax ** 2))

P_t0 = P_ISA * (1 + ((gamma-1)/2)*M_inf**2) ** ((gamma-1)/gamma)
T_t0 = T_ISA * (1 + ((gamma-1)/2)*M_inf**2)

P_t1 = P_ISA * (1 + ((gamma-1)/2)*M_ax**2) ** ((gamma-1)/gamma)
T_t1 = T_ISA * (1 + ((gamma-1)/2)*M_ax**2)

dT_t = (T_t1/eta_initial)*(beta_tt**((gamma-1)/gamma) - 1)

#Multall prints:
print("enthalpy = ",dT_t*1005)
print("P_t1 = ",P_t0, P_ISA)
print("T_t1 = ",T_t1)


## Area calculations
rho_1 = rho_ISA * ( (1 + 0.5 * (gamma-1) * M_ax**2) / (1 + 0.5 * (gamma-1) * M_inf**2) )**(-1/(gamma-1))
T_1 = T_ISA * ( (1 + 0.5 * (gamma-1) * M_ax**2) / (1 + 0.5 * (gamma-1) * M_inf**2) )**(-1)
a_1 = np.sqrt(gamma * R * T_1)
A_fan = f/(rho_1 * a_1 * M_ax)
print(A_fan)
print(rho_1)

R_t = np.sqrt(A_fan/(np.pi * (1 -  hub_to_tip**2)))
print(R_t)

##tip mach number
omega = RPM * 2 * np.pi /  60
U_tip = R_t * omega
M_tip = U_tip / a_1
print(M_tip)



