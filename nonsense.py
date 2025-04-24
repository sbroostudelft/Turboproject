import numpy as np
from scipy.optimize import fsolve


def shock_loss_model(beta_b_deg, M_1, i_deg, S_1, R, kappa=1.4):
    """
    Solves the complete shock-loss model system from Schobeiri (1998)

    Parameters:
    beta_b_deg : blade metal angle [degrees]
    M_1 : inlet Mach number
    i_deg : incidence angle [degrees]
    S_1 : blade spacing [m]
    R : curvature radius of mean flow path [m]
    kappa : specific heat ratio (default 1.4)

    Returns:
    dict: Contains all calculated quantities including shock loss coefficient
    """

    # Convert angles to radians
    beta_b = np.radians(beta_b_deg)
    i = np.radians(i_deg)

    # Prandtl-Meyer function (Eq. 5)
    def prandtl_meyer(M, kappa):
        term1 = np.sqrt((kappa + 1) / (kappa - 1))
        term2 = np.arctan(np.sqrt(((kappa - 1) / (kappa + 1)) * (M ** 2 - 1))
        term3 = np.arctan(np.sqrt(M ** 2 - 1))
        return term1 * term2 - term3

    # Function to solve the system of equations
    def equations(vars):
        M_s, theta, gamma = vars

        # Calculate Prandtl-Meyer angles
        nu_1 = prandtl_meyer(M_1, kappa)
        nu_s = prandtl_meyer(M_s, kappa)

        # Eq. 4: incidence angle relation
        eq4 = i - (theta - nu_s + nu_1)

        # Eq. 3: delta angle definition
        beta_s = beta_b + theta
        delta = (np.pi / 2) + beta_s - gamma

        # Eq. 10: shock angle relation
        numerator = np.cos(beta_s) - np.cos(beta_s + theta)
        denominator = -np.sin(beta_s) + np.sin(beta_s + theta) - (S_1 / (2 * R))
        eq10 = np.tan(gamma) - numerator / denominator

        # Eq. 7: momentum balance (simplified)
        lhs = M_1 * np.sqrt(1 + ((kappa - 1) / 2) * M_s ** 2)
        rhs = M_s * np.sqrt(1 + ((kappa - 1) / 2) * M_1 ** 2) * np.cos(beta_b) / np.cos(beta_s)
        eq7 = lhs - rhs

        return [eq4, eq10, eq7]

    # Initial guesses (M_s, theta, gamma)
    initial_guess = [M_1, np.radians(5), np.radians(130)]

    # Solve the system
    solution = fsolve(equations, initial_guess)
    M_s, theta_rad, gamma_rad = solution

    # Calculate remaining quantities
    theta_deg = np.degrees(theta_rad)
    gamma_deg = np.degrees(gamma_rad)
    beta_s = beta_b + theta_rad
    delta = (np.pi / 2) + beta_s - gamma_rad
    delta_deg = np.degrees(delta)

    # Shock loss coefficient (Eq. 12)
    term1 = ((kappa + 1) * M_s ** 2 * np.cos(delta) ** 2) / (2 + (kappa - 1) * M_s ** 2 * np.cos(delta) ** 2)
    term2 = 1 + (2 * kappa / (kappa + 1)) * (M_s ** 2 * np.cos(delta) ** 2 - 1)
    zeta_s = 1 - (term1 ** (kappa / (kappa - 1))) * (term2 ** (-1 / (kappa - 1)))

    return {
        'shock_Mach': M_s,
        'expansion_angle_deg': theta_deg,
        'shock_angle_deg': gamma_deg,
        'delta_angle_deg': delta_deg,
        'shock_loss_coefficient': zeta_s,
        'relative_flow_angle_deg': np.degrees(beta_s)
    }


# Example usage
if __name__ == "__main__":
    # Input parameters (example from paper)
    beta_b_deg = 30  # blade metal angle [deg]
    M_1 = 1.2  # inlet Mach number
    i_deg = 0  # incidence angle [deg]
    S_1 = 0.05  # blade spacing [m]
    R = 0.1  # curvature radius [m]

    results = shock_loss_model(beta_b_deg, M_1, i_deg, S_1, R)

    print("Shock Loss Model Results:")
    for key, value in results.items():
        print(f"{key:25s}: {value:.4f}")