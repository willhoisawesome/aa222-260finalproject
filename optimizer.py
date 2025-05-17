# optimizer.py

import numpy as np
from sizing import aa_260_engine_model, aa_260_sizing

def main():
    # Design variable ranges
    M_cruise = np.arange(0.6, 1.0, 0.1)        # Cruise Mach numbers
    alt = np.arange(9000, 13000, 1000)         # Altitudes in meters
    W_S = np.arange(50, 210, 10)               # Wing loading [lb/ft^2]
    pax = np.arange(100, 210, 10)              # Passenger count

    # Engine and atmospheric modeling
    turbofan, atmosphere = aa_260_engine_model(alt, M_cruise)

    # Sizing model
    w0, w_fuel, fuel_per_passenger = aa_260_sizing(
        turbofan, atmosphere, M_cruise, alt, W_S, pax
    )

if __name__ == "__main__":
    main()
