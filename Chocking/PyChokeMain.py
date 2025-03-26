import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
import PyChokeFunctions as pyf

# Estimation of axial compressor choking point
# Authors: Dr. ir. Andrea Giuffre', Dr. ir. Matteo Pini
# Reference: C. Freeman, N.A. Cumpsty - "A Method for the Prediction of Supersonic Compressor Blade Performance", 1989
# Delft University of Technology - All rights reserved

# compressor design data
beta_blade = np.array([-60.4,-60.4,-60.4])       # hub, mid, tip
R_hub = 0.168
R_tip = 0.3
t_le = 0.0014
Nbl = 48

# operating conditions
Pt = 276000
Tt = 340.0
mass_flow_vec = np.arange(70,90,1)     # range of mass flow to analyze
rpm = 13777

# fluid information (perfect gas)
gamma = 1.4
Rgas = 287.05
cp = 1004.5

# choose hub [0], mid [1] or tip [2]
loc = 1

# preliminary calculations
R_mean = (R_hub + R_tip) / 2
H_in = R_tip - R_hub
beta_blade = np.deg2rad(beta_blade)
Radius = np.array([R_hub, R_mean, R_tip])
throat = 2 * np.pi * Radius / Nbl * np.cos(beta_blade) - t_le
t_th = t_le / throat

incidence_vec = np.array([])
Mach_rel_in_vec = np.array([])
Mach_rel_out_vec = np.array([])
ds_vec = np.array([])
residual_vec = np.array([])

# run main loop
for mass_flow in mass_flow_vec:
    Mach_rel_in, beta, T_in, P_in = pyf.compute_inlet_flow_conditions(mass_flow, Pt, Tt, H_in, rpm, Rgas, gamma, Radius[loc])
    incidence = np.rad2deg(np.abs(beta) - np.abs(beta_blade[loc]))
    result = opt.fsolve(pyf.compute_post_shocks_flow_conditions, (Mach_rel_in),
                        args=(t_th[loc], beta, beta_blade[loc], Mach_rel_in, gamma), full_output=True)
    Mach_rel_out = result[0]
    rho_in = P_in / (Rgas * T_in)
    ds = pyf.compute_entropy_generation(beta_blade[loc], beta, t_th[loc], T_in, rho_in, P_in, Mach_rel_in, Mach_rel_out, Rgas, gamma, cp)

    incidence_vec = np.append(incidence_vec, incidence)
    Mach_rel_in_vec = np.append(Mach_rel_in_vec, Mach_rel_in)
    Mach_rel_out_vec = np.append(Mach_rel_out_vec, Mach_rel_out)
    ds_vec = np.append(ds_vec, ds)
    residual_vec = np.append(residual_vec, result[1]['fvec'])

choking_idx = np.argmax(residual_vec > 1e-6) - 1

# plot charts
plt.figure()
plt.plot(mass_flow_vec, incidence_vec, 'black')
plt.plot(mass_flow_vec[choking_idx], incidence_vec[choking_idx], 'ro')
plt.xlabel('mass flow rate [kg/s]')
plt.ylabel('incidence angle [deg]')
plt.grid(1)

plt.figure()
plt.plot(mass_flow_vec, Mach_rel_in_vec, 'black')
plt.plot(mass_flow_vec[choking_idx], Mach_rel_in_vec[choking_idx], 'ro')
plt.xlabel('mass flow rate [kg/s]')
plt.ylabel('inlet Mach [-]')
plt.grid(1)

plt.figure()
plt.plot(mass_flow_vec, Mach_rel_out_vec, 'black')
plt.plot(mass_flow_vec[choking_idx], Mach_rel_out_vec[choking_idx], 'ro')
plt.xlabel('mass flow rate [kg/s]')
plt.ylabel('post shocks Mach [-]')
plt.grid(1)

plt.figure()
plt.plot(mass_flow_vec, ds_vec, 'black')
plt.plot(mass_flow_vec[choking_idx], ds_vec[choking_idx], 'ro')
plt.xlabel('mass flow rate [kg/s]')
plt.ylabel('entropy rise [J/(kg K)]')
plt.grid(1)

plt.figure()
plt.plot(mass_flow_vec, residual_vec, 'black')
plt.plot(mass_flow_vec[choking_idx], residual_vec[choking_idx], 'ro')
plt.xlabel('mass flow rate [kg/s]')
plt.ylabel('( LHS- RHS ) / LHS [-]')
plt.grid(1)
plt.show()
