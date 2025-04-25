import numpy as np
from loss_models import compute_eta

def compute_pressure_ratio(psi, M_tip, gamma, eta_stage):
    base = 1 + eta_stage * psi * M_tip**2 * (gamma - 1)
    return base ** (1 / ((gamma - 1) / gamma))

def iterate_with_R_constraint(beta_target, psi_low, psi_high, phi, t_c, M_tip, gamma=1.4, tol=1e-4, max_iter=100):
    for i in range(max_iter):
        psi = (psi_low + psi_high) / 2
        R = -psi / 2 + 1  # Axial inflow constraint
        eta = compute_eta(psi, phi, R, t_c)
        beta_actual = compute_pressure_ratio(psi, M_tip, gamma, eta)

        print(f"Iter {i}: psi = {psi:.6f}, R = {R:.6f}, eta = {eta:.5f}, beta_tt_actual = {beta_actual:.5f}")

        if abs(beta_actual - beta_target) < tol:
            return psi, R, eta, beta_actual

        if beta_actual > beta_target:
            psi_high = psi
        else:
            psi_low = psi

    raise RuntimeError("Did not converge.")

if __name__ == "__main__":
    beta_target = 1.6
    psi_low = 0.2
    psi_high = 0.5
    phi = 0.4065
    t_c = 0.1
    M_tip = 0.98

    psi_opt, R_opt, eta_opt, beta_tt = iterate_with_R_constraint(beta_target, psi_low, psi_high, phi, t_c, M_tip)

    print("\n=== Final Results ===")
    print(f"psi: {psi_opt:.6f}")
    print(f"R: {R_opt:.6f}")
    print(f"Stage efficiency: {eta_opt:.4f}")
    print(f"Actual beta_tt: {beta_tt:.4f}")
