"""Attempt at implementing the Lieblein model for profile losses
    UNFINISHED
    """


import numpy as np

def lieblein(theta_over_c_2, solidity, beta_1, beta_2, H_2=1.08):
    omega = (
        2 * theta_over_c_2 * solidity / np.cos(beta_2) *
        (np.cos(beta_1) / np.cos(beta_2))**2 *
        (2 * H_2 / (3 * H_2 - 1)) *
        (1 - theta_over_c_2 * solidity * H_2 / np.cos(beta_2))**(-3)
    )
    return omega

def compute_eta(psi, phi, R, zeta_R, zeta_S):
    term1 = zeta_R * (phi**2 + (R + psi / 2)**2)
    term2 = zeta_S * (phi**2 + (1 - R + psi / 2)**2)
    eta = 1 - (1 / (2 * psi)) * (term1 + term2)
    return eta

if __name__ == "__main__":
    psi = 0.22
    phi = 0.4065
    R = 0.8877
    M1 = 1.4

    from meanline import computeVelocityTrianglesWithRKnown

    alpha1, alpha2, beta_1, beta_2 = computeVelocityTrianglesWithRKnown(psi, phi, R)
    print(beta_1, beta_2)
    pressure_loss = lieblein(0.02, 1.5,beta_1, beta_2)
    print(pressure_loss)
    "Convert to enthalpy based"
    pressure_loss = 0.015
    zeta_R = pressure_loss * (1 + 0.5 *(1.4 * M1) ** 2)
    print(zeta_R)
    efficiency = compute_eta(psi,phi, R, zeta_R, 0 )
    print(efficiency)