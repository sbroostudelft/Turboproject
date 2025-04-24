import numpy as np
from loss_models import compute_eta

def compute_beta_tt(psi, M_tip, gamma, eta):
    base = 1 + eta * psi * M_tip**2 * (gamma - 1)
    return base ** (1 / ((gamma - 1) / gamma))

def optimize_over_psi_phi(beta_target, psi_range, phi_range, t_c, M_tip, gamma=1.4):
    best_result = None
    min_error = float('inf')

    for psi in psi_range:
        R = -psi / 2 + 1
        for phi in phi_range:
            eta = compute_eta(psi, phi, R, t_c)
            beta_tt = compute_beta_tt(psi, M_tip, gamma, eta)
            error = abs(beta_tt - beta_target)

            if error < min_error:
                min_error = error
                best_result = (psi, phi, R, eta, beta_tt)

    return best_result

if __name__ == "__main__":
    beta_target = 1.6
    M_tip = 0.98
    t_c = 0.1

    psi_values = np.linspace(0.2, 0.5, 50)
    phi_values = np.linspace(0.35, 0.5, 50)

    psi_opt, phi_opt, R_opt, eta_opt, beta_tt = optimize_over_psi_phi(beta_target, psi_values, phi_values, t_c, M_tip)

    print("\n=== Optimal Parameters ===")
    print(f"psi: {psi_opt:.6f}")
    print(f"phi: {phi_opt:.6f}")
    print(f"R: {R_opt:.6f}")
    print(f"eta: {eta_opt:.4f}")
    print(f"beta_tt: {beta_tt:.4f}")
