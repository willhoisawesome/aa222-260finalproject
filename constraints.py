import numpy as np

# --- Constraint Parameters (default values can be overridden) ---
CD0 = 0.022 + 0.02 # from flaps
AR = 10.5
e = 0.82 -
TOP = 175
sigma_takeoff = 0.794
CL_to = 2.0
V_stall = 150 * 1.68781  # convert from knots to ft/s
CL_max_stall = 2.5
CL_max_land = 2.5
sigma_land = 1.0
S_landing = 6500
S_a = 1000
G_climb = 0.0883 # AEO, Appendix F.4

def stall_constraint(V_stall=V_stall, CL_max=CL_max_stall):
    """Returns W/S limit [lb/ft²] based on stall speed and density in slugs/ft³."""
    rhosealevel = 1.225  # Sea level density [kg/m³]
    rho_slug = rhosealevel * 0.00194032  # Convert from kg/m³ to slugs/ft³
    return 0.5 * rho_slug * V_stall**2 * CL_max

def takeoff_constraint(W_S_vals, TOP=TOP, sigma=sigma_takeoff, CL_to=CL_to):
    """Returns required T/W array for each W/S [lb/ft²] to meet takeoff performance."""
    return W_S_vals / (TOP * sigma * CL_to)

def landing_constraint(sigma=sigma_land, CL_max=CL_max_land, S_landing=S_landing, S_a=S_a):
    """Returns max allowable W/S [lb/ft²] to meet landing distance requirements."""
    return sigma * CL_max * (S_landing - S_a) / 80

def cruise_constraint(W_S_vals, rho, V, CD0=CD0, AR=AR, e=e):
    """Returns required T/W array for level cruise based on drag polar and W/S."""
    q = 0.5 * rho * V**2
    CL = W_S_vals / q
    LD = 1 / (CD0 + (CL**2 / (np.pi * AR * e)))
    return 1 / LD

def climb_constraint_vectorized(T_W_grid, rho, V, CD0=CD0, AR=AR, e=e, G=G_climb):
    """
    Returns:
    - W_S_grid: Wing loading grid [lb/ft²] that satisfies climb constraint for given T/W values
    - T_W_min: Minimum T/W needed for real solution (for plotting horizontal bound)
    """
    q = 0.5 * rho * V**2
    K = 1 / (np.pi * AR * e)
    T_W_min = G + np.sqrt(4 * CD0 * K)

    term = (T_W_grid - G)**2 - 4 * CD0 * K
    valid = term >= 0
    W_S_grid = np.full_like(T_W_grid, np.nan)
    W_S_grid[valid] = ((T_W_grid[valid] - G) + np.sqrt(term[valid])) / (2 / (q * np.pi * AR * e)) * 0.0208854

    return W_S_grid, T_W_min
