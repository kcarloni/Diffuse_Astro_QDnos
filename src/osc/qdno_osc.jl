
calc_λ( E, δmsq ) = 4E/δmsq #unnatural( Gpc, 4E/δmsq )
calc_δmsq( E, λ ) = uconvert( eV^2, 4E/natural(λ) ) 

"""
    calc_oscprob_w_decoherence( z, λ_osc, E )

Calculates the oscillation probability, damped by the decoherence probability. The decoherence effect is often negligible; see `calc_oscprob` for the approximate version.
"""
function calc_oscprob_w_decoherence( z, λ_osc, E )

    λ_coh = λ_osc * natural(E * 1u"fm") * sqrt(2)

    L2 = Leff(z, 2)
    L3 = Leff(z, 3)

    cos_term = cos( 2 * L2/λ_osc )

    # completely negligible because only relevant after oscillations have already averaged out 
    exp_term = exp( - L3^2 / λ_coh^2 )

    return 1/2 * (1 + exp_term * cos_term )
end

"""
    calc_oscprob( z, λ_osc; approx=true )

Calculates cos²( L₂(z)/λ ). If `approx=true`, then high-frequency oscillations, ie. L₂ ≪ λ, are averaged out to 1/2. 
"""
function calc_oscprob( z, λ_osc; approx=true )
    L2 = Leff(z, 2)

    if approx
        cos_term = 1/2
        if λ_osc / 1Gpc > 1/(4 * 2pi * acos(1/sqrt(2)))
            cos_term = cos( L2/λ_osc )^2
        end
    else
        cos_term = cos( L2/λ_osc )^2
    end

    return cos_term
end


function calc_oscprob_ϕtot( E, λ, src_pop )
    quadgk( 
        z -> dϕ_dz( src_pop, E, z ) * calc_oscprob(z, λ),
        src_pop.zmin,
        src_pop.zmax;
        rtol=1e-4
    )[1]
end