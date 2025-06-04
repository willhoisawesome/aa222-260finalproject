# sizing.py

# contains aircraft sizing model, engine model, and atmosphere model 

import numpy as np

def aa_260_sizing(turbofan, atmosphere, M_cruise_range, alt_range, W_S_range, passenger_count_range, range_nm, AR=10.5):
    # Constants
    a, b = -0.90, 1.36
    C1, C2, C3, C4, C5 = -0.10, 0.08, 0.05, -0.05, 0.20
    cd0 = 0.022
    e = 0.82
    range_ft = range_nm * 6076

    n_mach = len(M_cruise_range)
    n_alt = len(alt_range)
    n_ws = len(W_S_range)
    n_pax = len(passenger_count_range)

    hp_W0 = 0.036 * turbofan['v'] ** 0.32  # Empirical fit

    w0 = np.zeros((n_mach, n_alt, n_ws, n_pax))
    w_fuel = np.zeros_like(w0)
    fuel_per_passenger = np.zeros_like(w0)

    for i in range(n_ws):
        for j in range(n_mach):
            for k in range(n_alt):
                V = turbofan['v'][j, k]
                TSFC = turbofan['TSFC'][j, k]
                rho = atmosphere['rho'][k]

                coef1 = b * AR**C2 * hp_W0[j, k]**C3 * W_S_range[i]**C4 * V**C5
                coef2 = a

                q = 0.5 * rho / 515.4 * V**2  # dynamic pressure, slug/ft^3
                L_D = 1. / ((q * cd0 / W_S_range[i]) + (W_S_range[i] / (q * np.pi * AR * e)))

                wf_frac = np.exp(-range_ft * TSFC / 3600 / (V * L_D))

                w_frac_total = 0.98 * (0.991 - 0.007 * M_cruise_range[j] - 0.01 * M_cruise_range[j]**2) * wf_frac * 0.992 * 0.995
                wf_w0 = 1.08 * (1 - w_frac_total)

                for m in range(n_pax):
                    w_payload = 220 * passenger_count_range[m]

                    diff = 1.0
                    count = 0
                    w0_guess = 175000.0
                    while diff > 1e-3:
                        if count > 0:
                            w0_guess = 0.2 * (w0_guess - w0_calc) + w0_calc
                        w_empty = ((coef1 * w0_guess**C1 + coef2) * w0_guess)
                        w_f = wf_w0 * w0_guess
                        w0_calc = w_f + w_payload + w_empty
                        diff = abs((w0_guess - w0_calc) / w0_guess)
                        count += 1

                    w0[j, k, i, m] = w0_calc
                    w_fuel[j, k, i, m] = w_f
                    fuel_per_passenger[j, k, i, m] = w_f / passenger_count_range[m]

    return w0, w_fuel, fuel_per_passenger

def aa_260_engine_model(alt_array, mach_array):
    """
    Simulates the engine behavior (thrust, TSFC, fuel flow) over a grid of altitudes and Mach numbers.

    Parameters:
    alt_array (1D np.ndarray): Altitudes in meters
    mach_array (1D np.ndarray): Cruise Mach numbers

    Returns:
    dict: Engine output (thrust T, TSFC, fuel mass flow mdot_f, velocity v)
    dict: Atmosphere properties (T, t_ratio, rho, r_ratio)
    """
    # Constants and engine specs
    TSFC0 = 0.40               # lb/hr/lbf
    n = 0.8                    # TSFC exponent
    ne = 2                     # number of engines
    T0 = 26345                 # max sea-level thrust, lbf
    N1 = 0.8                   # throttle setting (percent)
    m = 0.7                    # density exponent
    a1 = -8.5e-4               # thrust curve coefficients
    a2 = 5.5e-7

    gamma = 1.4
    R = 287

    # Prepare output arrays
    nM = len(mach_array)
    nAlt = len(alt_array)

    T = np.zeros((nM, nAlt))
    TSFC = np.zeros_like(T)
    v = np.zeros_like(T)

    # Compute atmospheric data at each altitude
    T_atm = np.zeros(nAlt)
    t_ratio = np.zeros(nAlt)
    rho = np.zeros(nAlt)
    r_ratio = np.zeros(nAlt)

    for j, alt in enumerate(alt_array):
        atm = isentropic_atmosphere(alt)
        T_atm[j] = atm['T']
        t_ratio[j] = atm['t_ratio']
        rho[j] = atm['rho']
        r_ratio[j] = atm['r_ratio']

    # Calculate thrust, TSFC, fuel flow, velocity
    for j in range(nAlt):
        a = np.sqrt(gamma * R * T_atm[j])  # speed of sound, m/s
        for i in range(nM):
            v[i, j] = mach_array[i] * a * 3.28084  # m/s → ft/s
            TSFC[i, j] = TSFC0 * np.sqrt(t_ratio[j]) * (1 + mach_array[i])**n
            T[i, j] = ne * N1 * T0 * r_ratio[j]**m * (1 + a1 * v[i, j] + a2 * v[i, j]**2)

    mdot_f = TSFC * T  # fuel flow in lb/hr

    output = {
        'T': T,
        'TSFC': TSFC,
        'mdot_f': mdot_f,
        'v': v
    }

    atmosphere = {
        'T': T_atm,
        't_ratio': t_ratio,
        'rho': rho,
        'r_ratio': r_ratio
    }

    return output, atmosphere

def isentropic_atmosphere(alt):
    """
    Compute isentropic atmospheric conditions at a given altitude.

    Parameters:
    alt (float or np.ndarray): Altitude in meters

    Returns:
    dict: Dictionary containing pressure (Pa), temperature (K), and density (kg/m^3),
          as well as their ratios relative to sea level.
    """
    # Constants
    gamma = 1.4  # Ratio of specific heats for air
    p_surf = 101325  # Surface pressure, Pa
    rho_surf = 1.23  # Surface density, kg/m^3
    g = 9.81  # Gravity, m/s^2
    T_surf = 288  # Surface temperature, K
    R = 287  # Specific gas constant for air, J/kg/K
    z_star = R * T_surf / g  # Scale height, m

    # Local pressure
    p = ((1 - ((gamma - 1) / gamma) * (alt / z_star)) ** (gamma / (gamma - 1))) * p_surf
    p_ratio = p / p_surf

    # Local temperature
    T = (1 - ((gamma - 1) / gamma) * (alt / z_star)) * T_surf
    t_ratio = T / T_surf

    # Local density
    rho = ((1 - ((gamma - 1) / gamma) * (alt / z_star)) ** (1 / (gamma - 1))) * rho_surf
    r_ratio = rho / rho_surf

    return {
        'p': p,
        'p_ratio': p_ratio,
        'T': T,
        't_ratio': t_ratio,
        'rho': rho,
        'r_ratio': r_ratio
    }