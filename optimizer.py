import numpy as np
from scipy.optimize import differential_evolution
from sizing import aa_260_engine_model, aa_260_sizing
from constraints import (
    stall_constraint,
    landing_constraint,
    climb_constraint_vectorized,
    cruise_constraint,
    CD0, AR, G_climb,
    e as oswald_efficiency  # avoid naming conflict
)

# Design variable bounds: [Mach, Alt (m), Wing Loading (lb/ft²), Passengers]
bounds = [
    (0.6, 0.9),       # Mach number
    (9000, 12000),    # Altitude in meters
    (50, 200),        # Wing loading in lb/ft²
    (100, 200)        # Passenger count
]

iteration_counter = {'count': 0}  # mutable counter for callback logging

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
            range_nm=2000
        )
        fuel = fuel_pp[0, 0, 0, 0]
        w0_val = w0[0, 0, 0, 0]

    except Exception as e:
        print(f"[Sizing Failure] Mach={Mach:.3f}, Alt={Alt:.1f}, W/S={W_S:.1f}, Pax={pax} → {e}")
        return 1e9  # large penalty

    # Initialize penalty
    penalty = 0

    # Stall constraint
    W_S_stall = stall_constraint
    if W_S > W_S_stall:
        penalty += 1e6
        print(f"[Constraint Violation] Stall: W/S={W_S:.1f} > {W_S_stall:.1f}")

    # Landing constraint
    W_S_land = landing_constraint()
    if W_S > W_S_land:
        penalty += 1e6
        print(f"[Constraint Violation] Landing: W/S={W_S:.1f} > {W_S_land:.1f}")

    # Cruise constraint
    T_W_cruise = cruise_constraint(np.array([W_S]), rho, V, CD0, AR, oswald_efficiency)[0]
    T_W_actual = turbofan['T'][0, 0] / w0_val
    if T_W_actual < T_W_cruise:
        penalty += 1e6
        print(f"[Constraint Violation] Cruise: T/W={T_W_actual:.3f} < {T_W_cruise:.3f}")

    # Climb constraint
    _, T_W_min_climb = climb_constraint_vectorized(
        np.array([[T_W_actual]]), rho, V, CD0, AR, oswald_efficiency, G_climb
    )
    if T_W_actual < T_W_min_climb:
        penalty += 1e6
        print(f"[Constraint Violation] Climb: T/W={T_W_actual:.3f} < {T_W_min_climb:.3f}")

    if penalty > 0:
        print(f"[Penalty Triggered] Mach={Mach:.3f}, Alt={Alt:.1f}, W/S={W_S:.1f}, Pax={pax}, Penalty={penalty:.1e}")

    return fuel + penalty

def callback(xk, convergence):
    iteration_counter['count'] += 1
    Mach, Alt, W_S, pax = xk
    print(f"[Iter {iteration_counter['count']:2d}] "
          f"Mach={Mach:.3f}, Alt={Alt:.1f} m, W/S={W_S:.1f} lb/ft², Pax={int(np.round(pax))}")

if __name__ == "__main__":
    result = differential_evolution(
        objective,
        bounds,
        strategy='best1bin',
        maxiter=50,
        callback=callback,
        polish=True
    )

    print("\nOptimal Design:")
    print(f"Mach          : {result.x[0]:.3f}")
    print(f"Altitude [m]  : {result.x[1]:.1f}")
    print(f"Wing Loading  : {result.x[2]:.1f} lb/ft²")
    print(f"Passengers    : {int(np.round(result.x[3]))}")
    print(f"Fuel/Passenger: {result.fun:.2f} lb")
