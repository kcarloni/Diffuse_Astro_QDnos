

const pCMB_Planck2018 = (
    H0=67.36 * 1km/Mpc*1/s, Ωm=0.3153, ΩΛ=0.6847)

# speed of light in m/s
const c_light = (c/(1*m/s) * m/s)

H(z, pCMB=pCMB_Planck2018) = pCMB.H0 * sqrt( pCMB.Ωm * (1+z)^3 + pCMB.ΩΛ + ( 1 - pCMB.Ωm - pCMB.ΩΛ )*(1+z)^2 )

# ===========================
# distance metrics:

"""
    D_co(z) 

Co-moving distance, from cosmology.
"""
D_co(z) = uconvert(Mpc, c_light * quadgk( (z1 -> 1/H(z1)), 0, z )[1])

"""
    D_co(z) 

Luminosity distance, from cosmology. Related to `D_co(z)` by `D_lumi(z) = (1+z) D_co(z)`
"""
D_lumi(z) = (1 + z) * D_co(z)

"""
    Leff(z, n)

*Effective* distance, related to the total change in QM phase.
`Leff(z, n) = ∫ dz 1 /( H(z) * (1+z)^n )`
"""
Leff(z, n) = uconvert(Mpc, c_light .* quadgk( 
    (z1 -> (H(z1) * (1+z1)^n)^(-1) 
    ), 0, z )[1]) 

# ===========================
# derivatives of distances:

"""
    dD_co_dz(z) 

Derivative of the co-moving distance `D_co` with respect to redshift z. 
"""
dD_co_dz(z) = c_light/H(z)


"""
    dD_lumi_dz(z) 

Derivative of the luminosity distance `D_lumi` with respect to redshift z. 
"""
dD_lumi_dz(z) = (1+z) * c_light/H(z) + c_light * quadgk( z1 -> 1/H(z1), 0, z)[1]

# ===========================
# volume elements:

dV_co_dz(z) = 4π * D_co(z)^2 * c_light/H(z)