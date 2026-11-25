# 1) base atmosphere density

import numpy as np

def baseline_density(h, solar_activity): # h is the altitude you want to find the density for
# this can be changed but you need to find the corresponding p0 and Hd value aswell
    ref_alt = 310 
    if solar_activity == "high":
        p0 = 7.65e-12
        Hd = 94
    elif solar_activity == "low":
        p0 = 3.37e-12
        Hd = 80
    else:
        raise ValueError("Choose 'high' or 'low'")
    return p0 * np.exp(-(h - ref_alt) / Hd)
    
rho = baseline_density(400, "high")
print(rho)

# 2) solar factor

def solar_factor(F107, F107_avg, solar_activity):
# F10.7 is the daily solar radio flux
# F10.7 average is the 81 day average
    if solar_activity == "high":
        a1 = 9.43396e-3
        a2 = -2.22615e-6
        P107_ref = 144.7

    elif solar_activity == "low":
        a1 = 2.08690e-2
        a2 = -9.76385e-5
        P107_ref = 79.7

    else:
        raise ValueError("Choose 'high' or 'low'")

    P107 = (F107 + F107_avg) / 2

    delta_P = P107 - P107_ref

    return 1 + a1 * delta_P + a2 * delta_P**2

# 3) geomagnetic factor

def geomagnetic_factor(V_sw, By, Bz, time_hours, solar_activity):

    # Model constants
    tau = 0.5

    if solar_activity == "high":
        m1 = 4.67775e-2
        m2 = 3.35777e-4
        E_m_ref = 1.6

    elif solar_activity == "low":
        m1 = 1.18627e-1
        m2 = -1.36904e-3
        E_m_ref = 1.1

    else:
        raise ValueError("Choose 'high' or 'low'")

    # Step 1: Magnetic field clock angle
    gamma = np.arctan2(np.abs(By), Bz)

    # Step 2: Instantaneous merging electric field
    E_m_prime = (
        (1 / 3000)
        * V_sw**(4/3)
        * (By**2 + Bz**2)**(1/3)
        * np.sin(gamma / 2)**(8/3)
    )

    # Step 3: Weighted average over previous 3 hours
    t = time_hours[-1]

    weights = np.exp((time_hours - t) / tau)

    numerator = np.trapezoid(
        E_m_prime * weights, time_hours
    )

    denominator = np.trapezoid(
        weights, time_hours
    )

    E_m = numerator / denominator

    # Step 4: Difference from reference electric field
    delta_E = E_m - E_m_ref

    # Step 5: Geomagnetic correction factor
    F_geomag = 1 + m1 * delta_E + m2 * delta_E**2

    return F_geomag