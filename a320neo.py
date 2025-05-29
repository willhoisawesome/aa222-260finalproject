# a320neo.py

# Model and constraint evaluation for the A320neo configuration

import numpy as np
from sizing import aa_260_engine_model, aa_260_sizing
from constraints import (
    plot_constraints,
    stall_constraint, landing_constraint,
    cruise_constraint, climb_constraint,
    CD0, AR, e as oswald_efficiency, G_climb
)

# --- A320neo Reference Parameters ---
M_cruise = 0.82         # Cruise Mach number
alt = 11900             # Cruise altitude [m]
W_S = 165               # Wing loading [lb/ft²]
pax = 194               # Passenger count
range_nm = 3500         # Nautical mile range

# --- Engine and Atmosphere Evaluation ---
turbofan, atmosphere = aa_260_engine_model(np.array([alt]), np.array([M_cruise]))

# --- Sizing Model ---
w0, w_fuel, fuel_pp = aa_260_sizing(
    turbofan, atmosphere,
    np.array([M_cruise]), np.array([alt]),
    np.array([W_S]), np.array([pax]),
    range_nm=range_nm
)

# --- Extract Key Values ---
rho = atmosphere['rho'][0]
V = turbofan['v'][0, 0]
T = turbofan['T'][0, 0]
w0_val = w0[0, 0, 0, 0]
T_W_actual = T / w0_val
fuel_total = w_fuel[0, 0, 0, 0]
fuel_per_pax = fuel_pp[0, 0, 0, 0]

w0_val_kg = w0_val * 0.45359237  # Convert from lb to kg
fuel_total_kg = fuel_total * 0.45359237  # Convert from lb to kg
wing_area = 1/W_S * 0.092903 * w0_val  # Convert from lb/ft² to m²
range_km = range_nm * 1.852  # Convert from nautical miles to km
T_W_max = 54000 / w0_val  # Max thrust-to-weight ratio at takeoff

# --- Display Sizing Results ---
print("===== A320neo Sizing Results =====")
print(f"Initial Weight [kg]: {w0_val_kg:.2f}")
print(f"Fuel Weight [kg]   : {fuel_total_kg:.2f}")
print(f"Fuel/Passenger [lb]: {fuel_per_pax:.2f}")
print(f"Thrust-to-Weight Cruise   : {T_W_actual:.3f}")
print(f"Thrust-to-Weight Takeoff   : {T_W_max:.3f}")
print(f"Wing Area [m²]   : {wing_area:.2f}")
print(f"Range [km]        : {range_km:.2f}")

# --- Evaluate Constraints ---
W_S_stall = stall_constraint()
W_S_land = landing_constraint()
T_W_cruise = cruise_constraint(np.array([W_S]), rho, V, CD0, AR, oswald_efficiency)[0]
# _, T_W_min_climb = climb_constraint(
#    np.array([[T_W_actual]]), rho, V, CD0, AR, oswald_efficiency, G_climb
# )

# --- Print Results ---
print("===== A320neo Constraint Validation =====")
print(f"Cruise Mach          : {M_cruise:.2f}")
print(f"Cruise Altitude [m]  : {alt}")
print(f"Wing Loading [lb/ft²]: {W_S}")
print(f"Passengers           : {pax}")
print()
print(f"Stall Constraint     : W/S = {W_S:.1f} ≤ {W_S_stall:.1f} → {'OK' if W_S <= W_S_stall else '❌ VIOLATION'}")
print(f"Landing Constraint   : W/S = {W_S:.1f} ≤ {W_S_land:.1f} → {'OK' if W_S <= W_S_land else '❌ VIOLATION'}")
print(f"Cruise Constraint    : T/W = {T_W_actual:.3f} ≥ {T_W_cruise:.3f} → {'OK' if T_W_actual >= T_W_cruise else '❌ VIOLATION'}")
# print(f"Climb Constraint     : T/W = {T_W_actual:.3f} ≥ {T_W_min_climb:.3f} → {'OK' if T_W_actual >= T_W_min_climb else '❌ VIOLATION'}")

# Call plotting
plot_constraints(
    W_S_design=W_S,
    T_W_design=T_W_actual,
    T_W_takeoff_point = T_W_max,
    CD0=CD0,
    AR=AR,
    e=oswald_efficiency,
    V_cruise=V,
    rho=rho,
    G_climb=G_climb
)