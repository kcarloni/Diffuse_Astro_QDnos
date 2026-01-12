
using QuadGK
include( "cosmo.jl" )

# =========================================================

# n = number of sources per differential co-moving volume dV_co
# g = flux function, eg. SPL
struct DiffuseSourcePopulation{F1,F2,Qz}
    n::F1
    g::F2
    zmin::Qz
    zmax::Qz
end
DiffuseSourcePopulation( n, g ) = DiffuseSourcePopulation( 
    n, g, n.xmin, n.xmax )

Base.broadcastable( pop::DiffuseSourcePopulation ) = Ref( pop )

"""
    dϕ_dz( pop::DiffuseSourcePopulation, E, z ) 

Integrand in the expression ϕ_tot = ∫ dz n(z) × dV_co/dz × ϕ(E,z); 
see notes for more details. 
"""
function dϕ_dz( pop::DiffuseSourcePopulation, E, z ) 
    pop.n(z) * dD_co_dz(z) * pop.g( E * (1+z) )
end

function ϕtot( pop::DiffuseSourcePopulation, E )
    quadgk( 
        z -> dϕ_dz( pop, E, z), 
        pop.zmin, 
        pop.zmax; 
        rtol=1e-4
    )[1]
end

# =========================================================

# NB: factor of 3 in C0 corresponds to per-flavor => all-flavor,
# so that osc. effects can then divide again by 1/3

# 2 params 
function flux_SPL( E; γ, ϕ0, E0=100TeV, )
    C0 = 3e-18 * 1/GeV * 1/s * 1/sr * 1/cm^2
    return ϕ0 * C0 * (E / E0)^(-γ)
end

# 4 params
function flux_BPL( E; γ0, γ1, E1, ϕ0, sm=0.1, E0=10TeV )
    C0 = 3e-18 * 1/GeV * 1/s * 1/sr * 1/cm^2

    if (E1 > E0)
        ϕ0b = ϕ0 * (E0/E1)^γ0
    else
        ϕ0b = ϕ0 * (E0/E1)^γ1
    end

    if (E <= E1)
        return C0 * ϕ0b * (E/E1)^(-γ0)
    else
        return C0 * ϕ0b * (E/E1)^(-γ1)
    end

    # return ϕ0 * C0 * (E/E0)^(-γ0) * ( 1 + (E/E1)^(1/sm) )^( -sm * (γ1-γ0) )
end

# 3 params
function flux_SPLcutoff( E; γ, Ecut, ϕ0, E0=100TeV )
    return flux_SPL( E; γ, ϕ0, E0 ) * exp( -E/Ecut )
end

# =========================================================
