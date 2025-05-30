import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import differential_evolution
# from sizing import aa_260_engine_model, aa_260_sizing
from sizingh import aa_260_engine_model, aa_260_sizing # For Hydrogen aircraft
from constraints import (
    stall_constraint,
    landing_constraint,
    takeoff_constraint,
    climb_constraint,
    cruise_constraint,
    plot_constraints,
    CD0, AR, G_climb,
    e as oswald_efficiency
)

# Design variable bounds: [Mach, Altitude (m), Wing Loading (lb/ft²), Passengers]
bounds = [
    (0.2, 0.9),        # Mach
    (1000, 15000),     # Altitude
    (10, 300),         # Wing Loading
    (100, 400)         # Passengers
]

iteration_counter = {'count': 0}

def objective(x):
    Mach, Alt, W_S, pax = x
    Mach = float(Mach)
    Alt = float(Alt)
    W_S = float(W_S)
    pax = int(np.round(pax))

    try:
        turbofan, atmosphere = aa_260_engine_model(np.array([Alt]), np.array([Mach]))
        rho = atmosphere['rho'][0]
        V = turbofan['v'][0, 0]

        w0, w_fuel, fuel_pp = aa_260_sizing(
            turbofan, atmosphere,
            np.array([Mach]), np.array([Alt]),
            np.array([W_S]), np.array([pax]),
            range_nm=1800
        )
        fuel = fuel_pp[0, 0, 0, 0]
        w0_val = w0[0, 0, 0, 0]
        fuel_weight_kg = w_fuel[0, 0, 0, 0] * 0.4536  # lb to kg

    except Exception as e:
        print(f"[Sizing Failure] Mach={Mach:.3f}, Alt={Alt:.1f}, W/S={W_S:.1f}, Pax={pax} → {e}")
        return 1e9

    penalty = 0
    rho = 1e10

    # --- Stall Constraint ---
    W_S_stall = stall_constraint()
    if W_S > W_S_stall:
        delta = W_S - W_S_stall
        penalty += rho * delta**2
        print(f"[Violation] Stall: W/S={W_S:.1f} > {W_S_stall:.1f}")

    # --- Landing Constraint ---
    W_S_land = landing_constraint()
    if W_S > W_S_land:
        delta = W_S - W_S_land
        penalty += rho * delta**2
        print(f"[Violation] Landing: W/S={W_S:.1f} > {W_S_land:.1f}")

    # --- Cruise Constraint ---
    T_W_cruise = cruise_constraint(np.array([W_S]), rho, V, CD0, AR, oswald_efficiency)[0]
    T_W_actual = turbofan['T'][0, 0] / w0_val
    if T_W_actual < T_W_cruise:
        delta = T_W_cruise - T_W_actual
        penalty += rho * delta**2
        print(f"[Violation] Cruise: T/W={T_W_actual:.3f} < {T_W_cruise:.3f}")

    # --- Climb Constraint ---
    _, T_W_min_climb = climb_constraint(
        np.array([[T_W_actual]]), V, CD0, AR, oswald_efficiency, G_climb
    )
    T_W_max = 54000 / w0_val  # estimated max thrust-to-weight ratio
    if T_W_max < T_W_min_climb:
        delta = T_W_min_climb - T_W_max
        penalty += rho * delta**2
        print(f"[Violation] Climb: T/W_max={T_W_max:.3f} < {T_W_min_climb:.3f}")

    # --- Takeoff Constraint ---
    T_W_takeoff_req = takeoff_constraint(np.array([W_S]))[0]
    if T_W_max < T_W_takeoff_req:
        delta = T_W_takeoff_req - T_W_max
        penalty += rho * delta**2
        print(f"[Violation] Takeoff: T/W_max={T_W_max:.3f} < {T_W_takeoff_req:.3f}")

    # --- Wingspan Constraint ---
    S_ft2 = w0_val / W_S  # Wing area in ft²
    b_ft = np.sqrt(AR * S_ft2)  # Wingspan in ft
    b_m = b_ft * 0.3048  # Wingspan in meters

    max_wingspan_m = 36.0
    if b_m > max_wingspan_m:
        delta = b_m - max_wingspan_m
        penalty += rho * delta**2
        print(f"[Violation] Wingspan: b = {b_m:.2f} m > {max_wingspan_m:.2f} m")

    # --- Fuselage Length Constraint ---
    seat_pitch = 0.76  # m
    seats_abreast = 6
    L_galley = 5       # m
    L_cabin = pax * seat_pitch / seats_abreast + L_galley

    R_tank = 2.0       # m
    rho_LH2 = 70.85       # kg/m³
    eta_V = 0.927
    eta_f = 0.85
    V_tank = fuel_weight_kg / (rho_LH2 * eta_V * eta_f)
    L_tank = V_tank / (np.pi * R_tank**2)

    L_total = L_cabin + L_tank + 4 + 4  # add 4m nose and 4m tail
    if L_total > 44.5:
        penalty += rho * (L_total - 44.5)**2
        print(f"[Violation] Fuselage Length: {L_total:.2f} m > 44.5 m")

    if penalty > 0:
        print(f"[Penalty] Mach={Mach:.3f}, Alt={Alt:.1f}, W/S={W_S:.1f}, Pax={pax}, Penalty={penalty:.1e}")

    return fuel + penalty

def callback(xk, convergence):
    iteration_counter['count'] += 1
    Mach, Alt, W_S, pax = xk
    print(f"[Iter {iteration_counter['count']:2d}] "
          f"Mach={Mach:.3f}, Alt={Alt:.1f} m, W/S={W_S:.1f}, Pax={int(np.round(pax))}")

if __name__ == "__main__":
    result = differential_evolution(
        objective,
        bounds,
        strategy='best1bin',
        maxiter=100,
        callback=callback,
        polish=True
    )

    Mach, Alt, W_S, pax = result.x
    pax = int(np.round(pax))

    turbofan, atmosphere = aa_260_engine_model(np.array([Alt]), np.array([Mach]))
    rho = atmosphere['rho'][0]
    V = turbofan['v'][0, 0]
    w0, w_fuel, fuel_pp = aa_260_sizing(
        turbofan, atmosphere,
        np.array([Mach]), np.array([Alt]),
        np.array([W_S]), np.array([pax]),
        range_nm=2000
    )

    w0_val = w0[0, 0, 0, 0]
    fuel_total = w_fuel[0, 0, 0, 0]
    T_W_cruise = turbofan['T'][0, 0] / w0_val
    T_W_takeoff = 54000 / w0_val

    print("\nOptimal Design:")
    print(f"Mach          : {Mach:.3f}")
    print(f"Altitude [m]  : {Alt:.1f}")
    print(f"Wing Loading  : {W_S:.1f} lb/ft²")
    print(f"Passengers    : {pax}")
    print(f"Fuel/Passenger: {fuel_total/pax:.2f} lb")
    print(f"MTOW          : {w0_val:.1f} lb")
    print(f"Fuel Weight   : {fuel_total:.1f} lb")
    print(f"T/W (Cruise)  : {T_W_cruise:.3f}")
    print(f"T/W (Takeoff) : {T_W_takeoff:.3f}")

    # Calculate wingspan
    S_ft2 = w0_val / W_S  # Wing area in ft²
    b_ft = np.sqrt(AR * S_ft2)  # Wingspan in ft
    b_m = b_ft * 0.3048  # Convert to meters

    print(f"Wing Area     : {S_ft2:.1f} ft²")
    print(f"Wingspan      : {b_ft:.1f} ft ({b_m:.2f} m)")

    # Fuselage Length Calculation
    seat_pitch = 0.76  # m
    seats_abreast = 6
    L_galley = 5       # m
    L_cabin = pax * seat_pitch / seats_abreast + L_galley

    R_tank = 2.0       # m
    rho_LH2 = 71       # kg/m³
    eta_V = 0.927
    eta_f = 0.85
    fuel_kg = fuel_total * 0.4536  # Convert lb to kg
    V_tank = fuel_kg / (rho_LH2 * eta_V * eta_f)
    L_tank = V_tank / (np.pi * R_tank**2)

    L_total = L_cabin + L_tank + 4 + 4  # add nose + tail

    print(f"Cabin Length   : {L_cabin:.2f} m")
    print(f"Tank Length    : {L_tank:.2f} m")
    print(f"Fuselage Length: {L_total:.2f} m")

    plot_constraints(
        W_S_design=W_S,
        T_W_design=T_W_cruise,
        T_W_takeoff_point=T_W_takeoff,
        CD0=CD0,
        AR=AR,
        e=oswald_efficiency,
        V_cruise=V,
        rho=rho,
        G_climb=G_climb
    )