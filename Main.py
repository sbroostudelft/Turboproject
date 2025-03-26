import numpy as np

rho_ISA = 0.413 # kg/m^3
T_ISA = 223.2 # K
R = 287.05 # J/kgK
P_ISA = rho_ISA * R * T_ISA # Pa
M_inf = 0.78
M_ax = 0.6
gamma = 1.4
beta_tt = 1.6
eta_initial = 0.9

P_t0 = P_ISA * (1 + ((gamma-1)/2)*M_inf**2) ** ((gamma-1)/gamma)
T_t0 = T_ISA * (1 + ((gamma-1)/2)*M_inf**2)

P_t1 = P_ISA * (1 + ((gamma-1)/2)*M_ax**2) ** ((gamma-1)/gamma)
T_t1 = T_ISA * (1 + ((gamma-1)/2)*M_ax**2)

dT_t = (T_t1/eta_initial)*(beta_tt**((gamma-1)/gamma) - 1)

print(dT_t*1005)