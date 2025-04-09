import numpy as np


def calculate_incidence_deflection():
    """Calculates duty coefficents based on the compression ratio and tip mach number
    while assuming axial inflow.
    """

    "Calculation of work coefficient with the dimensionless Euler equation"
    exponent: float = (gamma - 1) / gamma
    work_coeff: float = (beta_tt ** exponent - 1) / (M_u ** 2 * (gamma - 1))

    "Calculation of DoR with stage relations assuming axial inflow"
    degree_of_reaction: float  = -(work_coeff/2) + 1

    "Calculation of optimum flow coefficient based on smith charts"
    def negative_eff(phi):
        return -get_repeated_stage_efficiency(phi, 2*work_coeff, degree_of_reaction)

    "Finding optimum efficiency"
    optimum = minimize_scalar(negative_eff, bounds=(0.1, 1.0), method="bounded")
    flow_coeff = optimum.x
    max_eta = -optimum.fun

    print(f"Max efficiency: {max_eta:.4f}")
    print(f"psi: {work_coeff:.4f} phi: {flow_coeff:.4f} DOR: {degree_of_reaction:.4f}")

    return work_coeff, flow_coeff, degree_of_reaction

if __name__ == "__main__":

    #Calculate tip mach number assuming axial inflow
    M_v_1 = 0.6
    M_w_1 = 1.4
    M_u_1 = np.sqrt(M_w_1 ** 2 - M_v_1 ** 2)

    #Calculate work coefficients
    psi, phi, DOR = calculate_duty_coefficients(M_u_1)
