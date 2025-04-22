"""Attempt at implementing the Lieblein model for profile losses
    UNFINISHED
    """


import numpy as np
from meanline import computeVelocityTrianglesWithRKnown

def compute_denton_bl_losses(psi, phi, R, solidity):
    alpha_1, alpha_2, beta_1, beta_2 = computeVelocityTrianglesWithRKnown(psi, phi, R)
    C_d = 0.002
    delta_V_bar_v = 0.5 * (1/solidity) * (np.tan(alpha_2) - np.tan(alpha_1))
    delta_V_bar_v_2 = 0.5 * (1 / solidity) * (np.tan(np.abs(beta_2)) - np.tan(np.abs(beta_1)))
    zeta_s = C_d * (2 * (1/delta_V_bar_v) + 6 * (delta_V_bar_v)) * (np.tan(alpha_1) + np.tan(alpha_2))
    zeta_r = C_d * (2 * (1/delta_V_bar_v_2) + 6 * (delta_V_bar_v+2)) * (np.tan(np.abs(beta_1)) + np.tan(np.abs(beta_2)))
    return zeta_s, zeta_r

def compute_denton_shock_losses(psi, phi, R):
    return

def compute_TE_mixing_losses(t_over_s):
    return

def compute_eta(psi, phi, R, zeta_R, zeta_S):
    term1 = zeta_R * (phi**2 + (R + psi / 2)**2)
    term2 = zeta_S * (phi**2 + (1 - R + psi / 2)**2)
    eta = 1 - (1 / (2 * psi)) * (term1 + term2)
    return eta

if __name__ == "__main__":
    psi = 0.22
    phi = 0.4065
    R = 0.8877

    allah, allah2 = compute_denton_bl_losses(psi, phi, R, 1.5)
    print(allah, allah2)
    print(compute_eta(psi, phi, R, allah2, allah))
