"""
Calculates the loss coefficients and the efficiency
Considers 2D loss mechanisms: boundary layer, shock and mixing losses
"""

import numpy as np
from meanline import computeVelocityTrianglesWithRKnown

def compute_denton_bl_losses(psi, phi, R):
    '''
    Denton Boundary Layer loss model implementing
    '''
    alpha_1, alpha_2, beta_1, beta_2 = computeVelocityTrianglesWithRKnown(psi, phi, R)
    C_d = 0.002

    "Old method by Sam"
    #delta_V_bar_v = 0.5 * (1/solidity) * (np.tan(alpha_2) - np.tan(alpha_1))
    #delta_V_bar_v_2 = 0.5 * (1 / solidity) * (np.tan(np.abs(beta_2)) - np.tan(np.abs(beta_1)))

    "New method assuming optimal solidity based on Zweiffel"
    delta_V_bar_v = 1 / np.sqrt(3)
    delta_V_bar_v_2 = 1 / np.sqrt(3)
    zeta_s = C_d * (2 * (1/delta_V_bar_v) + 6 * (delta_V_bar_v)) * (np.tan(alpha_1) + np.tan(alpha_2))
    zeta_r = C_d * (2 * (1/delta_V_bar_v_2) + 6 * (delta_V_bar_v_2)) * (np.tan(np.abs(beta_1)) + np.tan(np.abs(beta_2)))
    return zeta_s, zeta_r

def compute_TE_mixing_losses(t_over_s):
    '''
    Makes use of the incompressible loss model
    '''
    zeta_s = (t_over_s ** 2) + 0.15 * t_over_s
    zeta_r = zeta_s
    return zeta_s, zeta_r

def compute_eta(psi, phi, R, thick_to_chord):

    zeta_s_bl, zeta_r_bl = compute_denton_bl_losses(psi, phi, R)
    zeta_s_te, zeta_r_te = compute_TE_mixing_losses(thick_to_chord)
    zeta_S = zeta_s_bl + zeta_s_te
    zeta_R = zeta_r_bl + zeta_r_te
    term1 = zeta_R * (phi**2 + (R + psi / 2)**2)
    term2 = zeta_S * (phi**2 + (1 - R + psi / 2)**2)
    eta = 1 - (1 / (2 * psi)) * (term1 + term2)
    return eta

if __name__ == "__main__":
    psi = 0.22
    phi = 0.4065
    R = 0.8877

    zeta_s_friction, zeta_r_friction = compute_denton_bl_losses(psi, phi, R)
    zeta_s_wake, zeta_r_wake = compute_TE_mixing_losses(0.1)
    eta = compute_eta(psi, phi, R, 0.1)
    print(f"Stator loss coefficient due to BL: {zeta_s_friction:.6f}")
    print(f"Rotor loss coefficient due to BL: {zeta_r_friction:.6f}")
    print(f"Stator loss coefficient due to TE: {zeta_s_wake:.6f}")
    print(f"Rotor loss coefficient due to TE: {zeta_r_wake:.6f}")
    print(f"Stator loss coefficient: {zeta_s_wake + zeta_s_friction:.6f}")
    print(f"Rotor loss coefficient: {zeta_r_wake + zeta_r_friction:.6f}")
    print(f"Overall stage efficiency: {eta:.4f}")

    psi = 0.45
    phi = 0.4065
    R = 0.8877

    zeta_s_friction, zeta_r_friction = compute_denton_bl_losses(psi, phi, R)
    zeta_s_wake, zeta_r_wake = compute_TE_mixing_losses(0.1)
    eta = compute_eta(psi, phi, R, 0.1)
    print(f"Stator loss coefficient due to BL: {zeta_s_friction:.6f}")
    print(f"Rotor loss coefficient due to BL: {zeta_r_friction:.6f}")
    print(f"Stator loss coefficient due to TE: {zeta_s_wake:.6f}")
    print(f"Rotor loss coefficient due to TE: {zeta_r_wake:.6f}")
    print(f"Stator loss coefficient: {zeta_s_wake + zeta_s_friction:.6f}")
    print(f"Rotor loss coefficient: {zeta_r_wake + zeta_r_friction:.6f}")
    print(f"Overall stage efficiency: {eta:.4f}")