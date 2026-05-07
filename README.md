# Dynamic Simulation of a Heated Multi-Component Vessel

> Process modeling project — *Modélisation des Opérations Unitaires* (ENSGTI) 

A FORTRAN simulator that models the dynamic behavior of a closed chemical
vessel containing a non-ideal ternary mixture (chloroform / methanol / acetone),
heated at constant power and equipped with a safety rupture disk.

The model couples mass and energy balances with vapor-liquid equilibrium and
handles the discontinuity introduced by the rupture event, allowing prediction
of pressure peak, rupture time, and the post-rupture quasi-steady regime.

## Physical system

A closed vessel containing a liquid mixture in equilibrium with its vapor is
heated at constant power Q. As the temperature rises, vapor accumulates and
pressure increases until it reaches the rupture threshold P_rup. The disk then
yields, and a vapor outflow E proportional to √(P − P_amb) escapes to the
ambient. The system eventually settles into a quasi-stationary regime where
heat input is exactly balanced by the latent heat of the escaping vapor.

## Mathematical model

The system is formulated as a differential-algebraic equation (DAE) system of
index 1 with N = 2·N_C + 6 = 12 state variables:

- liquid retention U^liq, mole fractions x_i, specific enthalpy h^liq
- vapor retention U^vap, mole fractions y_i
- temperature T, pressure P, leak flow E

Eight equations are implemented:
- partial mass balance (liquid)
- energy balance (liquid)
- vapor mass balance
- summation constraint
- ideal-gas equation of state
- liquid-vapor equilibrium (NRTL + modified Raoult)
- enthalpy constitutive equation
- rupture disk law (E = 0 before, E = K_disc·√(P − P_amb) after)

