import numpy as np
import matplotlib.pyplot as plt

# --- Constraint Parameters ---
CD0 = 0.022
AR = 10.5
e = 0.82
TOP = 175
sigma_takeoff = 1.0
CL_to = 2.5
V_stall = 150 * 1.68781  # knots to ft/s
CL_max_stall = 2.5
CL_max_land = 2.5
sigma_land = 1.0
S_landing = 6500
S_a = 1000
G_climb = 0.0883  # AEO climb gradient

# === Constraint Functions ===

def stall_constraint(V_stall=V_stall, CL_max=CL_max_stall):
    rho_sl = 1.225  # kg/m³
    rho_slug = rho_sl * 0.00194032  # slugs/ft³
    return 0.5 * rho_slug * V_stall**2 * CL_max

def landing_constraint(CL_max=CL_max_land, sigma=sigma_land, S=S_landing, S_a=S_a):
    return sigma * CL_max * (S - S_a) / 80

def takeoff_constraint(W_S_vals, TOP=TOP, sigma=sigma_takeoff, CL_to=CL_to):
    return W_S_vals / (TOP * sigma * CL_to)

def cruise_constraint(W_S_vals, rho, V, CD0=CD0, AR=AR, e=e):
    q = 0.5 * rho * V**2
    CL = W_S_vals / q
    LD = 1 / (CD0 + (CL**2 / (np.pi * AR * e)))
    return 1 / LD

def climb_constraint(T_W_vals, V, CD0=CD0, AR=AR, e=e, G=G_climb):
    rho_sl = 1.225  # kg/m³
    rho_slug = rho_sl * 0.00194032  # slugs/ft³
    q = 0.5 * rho_slug * V**2  
    K = 1 / (np.pi * AR * e)
    a = 2 / (q * np.pi * AR * e)
    delta = (T_W_vals - G)**2 - 4 * CD0 * K

    W_S = np.full_like(T_W_vals, np.nan)
    feasible = delta >= 0
    W_S[feasible] = ((T_W_vals[feasible] - G) + np.sqrt(delta[feasible])) / a

    T_W_min = G + np.sqrt(4 * CD0 * K)
    return W_S, T_W_min

# === Plotting Function ===

def plot_constraints(W_S_design, T_W_design, T_W_takeoff_point, CD0, AR, e, V_cruise, rho, G_climb):
    W_S = np.linspace(50, 250, 500)
    T_W_vals = np.linspace(0.05, 0.4, 500)
    q = 0.5 * rho * V_cruise**2

    # Constraints
    W_S_stall_val = stall_constraint()
    W_S_land_val = landing_constraint()
    T_W_takeoff = takeoff_constraint(W_S)
    T_W_cruise = cruise_constraint(W_S, rho, V_cruise, CD0, AR, e)
    W_S_climb, T_W_min = climb_constraint(T_W_vals, q, CD0, AR, e, G_climb)

    # --- Compute feasible region ---
    # Determine T/W lower bounds from max of climb, cruise, takeoff
    T_W_lower = np.maximum.reduce([T_W_takeoff, T_W_cruise])
    W_S_climb_clean = W_S_climb.copy()
    T_W_climb_clean = T_W_vals.copy()
    W_S_climb_clean = W_S_climb_clean[~np.isnan(W_S_climb_clean)]
    T_W_climb_clean = T_W_climb_clean[~np.isnan(W_S_climb)]

    if len(W_S_climb_clean) > 0:
        # Interpolate climb into common W/S space
        T_W_climb_interp = np.interp(W_S, W_S_climb_clean, T_W_climb_clean, left=np.nan, right=np.nan)
        T_W_lower = np.maximum(T_W_lower, T_W_climb_interp)

    W_S_upper = min(W_S_stall_val, W_S_land_val)

    # === Tight hatching bands for constraint direction cues ===
    band_width_w = 5
    band_width_t = 0.01

    # Takeoff (below curve is infeasible)
    plt.fill_between(W_S, T_W_takeoff - band_width_t, T_W_takeoff,
                     color='none', hatch='\\', facecolor='darkorange', alpha=0.1, zorder=1)

    # Cruise (below curve is infeasible)
    plt.fill_between(W_S, T_W_cruise - band_width_t, T_W_cruise,
                     color='none', hatch='//', facecolor='seagreen', alpha=0.1, zorder=1)

    # Climb (left of curve is infeasible)
    if len(W_S_climb_clean) > 0:
        plt.fill_betweenx(T_W_climb_clean,
                          W_S_climb_clean - band_width_w,
                          W_S_climb_clean,
                          color='none'
                          , hatch='//', facecolor='firebrick', alpha=0.1, zorder=1)

    # Climb (below curve is infeasible)
    plt.fill_between(W_S, T_W_min - band_width_t, T_W_min,
                     color='none', hatch='//', facecolor='firebrick', alpha=0.1, zorder=1)
    
    # Stall constraint band (right of line is infeasible)
    plt.axvspan(W_S_stall_val, W_S_stall_val + band_width_w,
                facecolor='rebeccapurple', hatch='//', alpha=0.1, zorder=1)

    # Landing constraint band
    plt.axvspan(W_S_land_val, W_S_land_val + band_width_w,
                facecolor='saddlebrown', hatch='//', alpha=0.1, zorder=1)
                
    # Plot constraint boundaries
    plt.plot(W_S, T_W_takeoff, color='darkorange', label="Takeoff Constraint")
    plt.plot(W_S, T_W_cruise, color='seagreen', label="Cruise Constraint")
    # plt.plot(W_S_climb, T_W_vals, color='firebrick', label="Climb Constraint (Raymer)")
    plt.hlines(T_W_min, xmin=50, xmax=W_S_climb_clean[0], color='firebrick', label="Climb Constraint")

    plt.axvline(W_S_stall_val, color='rebeccapurple', label="Stall Constraint")
    plt.axvline(W_S_land_val, color='saddlebrown', label="Landing Constraint")

    # Design point
    plt.plot(W_S_design, T_W_design, 'ro', label='Cruise Design Point')
    plt.plot(W_S_design, T_W_takeoff_point, 'ro', label='Takeoff Design Point')

    # Formatting
    plt.xlabel('Wing Loading W/S [lb/ft²]')
    plt.ylabel('Thrust-to-Weight T/W')
    plt.title('Constraint Boundaries and Feasible Region')
    plt.grid(True)
    plt.legend()
    plt.xlim(50, 250)
    plt.ylim(0, 0.5)
    plt.tight_layout()
    plt.show()