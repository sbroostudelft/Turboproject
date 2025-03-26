import numpy as np
import scipy.optimize as opt


def compute_inlet_flow_conditions(mass_flow, Pt, Tt, H_in, rpm, R, gamma, R_loc):
    """
    Given the operating point specified in terms of mass flow rate and rotational speed and
    the inlet geometry specified in terms of blade height, mean radius and blade angle,
    compute the inlet relative Mach number and the flow angle.
    """
    Area = 2 * np.pi * R_loc * H_in
    omega = rpm * (2 * np.pi) / 60
    U = omega * R_loc
    Mach_abs = opt.fsolve(massflow_isentropic, (0.5), args=(Area, mass_flow, Pt, Tt, R, gamma), full_output=False)
    T = Tt / (1 + (gamma - 1) / 2 * Mach_abs ** 2)
    P = Pt / ((Tt / T) ** (gamma / (gamma - 1)))
    Vm = Mach_abs * np.sqrt(gamma * R * T)
    Wm = Vm
    Wt = 0.0 - U
    W = np.sqrt(Wm ** 2 + Wt ** 2)
    Mach_rel = W / np.sqrt(gamma * R * T)
    beta = np.arctan(Wt / Wm)

    return Mach_rel, beta, T, P


def massflow_isentropic(p, *data):
    """
    Given the annulus area and the mass flow passing through it, compute the absolute Mach number, assuming
    isentropic flow and perfect gas.
    """
    A, m, Pt, Tt, R, gamma = data
    Mach = p

    m_computed = Pt * A / np.sqrt(R * Tt) * Mach * np.sqrt(gamma) * (1 + (gamma - 1) / 2 * Mach ** 2) ** \
                 ((1 + gamma) / (2 * (1 - gamma)))
    res = (m_computed - m) / m

    return res


def freeman_cv(t_th, beta, beta_blade, Mach_rel_in, Mach_rel_out, gamma):
    """
    Compute left-hand and right-hand sides of the equation derived by Freeman and Cumpsty for the estimation
    of choking point in axial compressor cascades, assuming perfect gas.
    Inlet and outlet refer to the boundaries of the control volume, not of the cascade.
    """
    lhs = ((1 + (gamma - 1) / 2 * (Mach_rel_out ** 2)) ** (- 1 / 2)) * \
          (1 + gamma * (Mach_rel_out ** 2) * (1 - t_th)) / (Mach_rel_out * (1 - t_th))
    rhs = ((1 + (gamma - 1) / 2 * (Mach_rel_in ** 2)) ** (- 1 / 2)) * \
          (np.cos(beta_blade) / np.cos(beta) + gamma * (Mach_rel_in ** 2) * np.cos(beta - beta_blade)) / Mach_rel_in
    
    return lhs, rhs


def compute_post_shocks_flow_conditions(p, *data):
    """
    Compute Mach number at the outlet of the control volume defined by Freeman and Cumpsty.
    """
    t_th, beta, beta_blade, Mach_rel_in, gamma = data
    Mach_rel_out = p

    lhs, rhs = freeman_cv(t_th, beta, beta_blade, Mach_rel_in, Mach_rel_out, gamma)
    res = (lhs - rhs) / lhs

    return res


def compute_entropy_generation(beta_blade, beta, t_th, T_in, rho_in, P_in, Mach_rel_in, Mach_rel_out, R, gamma, cp):
    """
    Given inlet state and outlet relative Mach number, compute entropy generation within the control volume,
    assuming perfect gas.
    """
    T_out = T_in * (1 + (gamma - 1) / 2 * Mach_rel_in ** 2) / (1 + (gamma - 1) / 2 * Mach_rel_out ** 2)
    rho_out = rho_in * (Mach_rel_in * np.sqrt(T_in) * np.cos(beta) / np.cos(beta_blade)) / \
              (Mach_rel_out * np.sqrt(T_out) * (1 - t_th))
    P_out = rho_out * R * T_out
    ds = cp * np.log(T_out / T_in) - R * np.log(P_out / P_in)

    return ds

