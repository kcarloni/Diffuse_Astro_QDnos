
using QuadGK

# base integral
"""
= fac × ∫ dE E⁻² A(E, β) f(E, β)
"""
function calc_int_effA_dE!( arr, binedges, effA, f; fac=1 )

    function f_to_int( u ) 

        E = exp10(u) * GeV
        dE = log(10) * E # * du 

        return dE * sum( 
            effA(E, β) * E^(-2) * f(E, β) for β in (ele, mu, tau)
        )

    end

    log_bine = log10.( binedges / 1GeV )
    inds = 2:length(log_bine)

    for i in inds
        arr[i-1] += quadgk( 
            f_to_int, 
            log_bine[i-1],
            log_bine[i]
        )[1] * fac 
    end
end

# handles composite datasets...
function calc_num_events!( N, dset, f )

    effA = dset.effA
    ΔT = dset.ΔT
    ΔΩ = dset.ΔΩ
    bine = dset.bine

    # expand
    if length(dset.ΔT) == 1 
        effA = ( effA, )
        ΔT = ( ΔT, )
    end

    # mutating...
    for j in eachindex(effA)
        calc_int_effA_dE!( 
            N, 
            bine,
            effA[j],
            f;
            fac = ΔT[j] * ΔΩ
        )
    end
end

# ------------------------

function calc_num_events( dset; effA_int=dset.effA_int )
    N = (
        med = natural.(dset.fdata.med .* effA_int),
        low = natural.(dset.fdata.low .* effA_int),
        up = natural.(dset.fdata.up .* effA_int),
    )
end

function calc_num_events( dset, f )
    N = zeros(
        typeof( 1/GeV * 1cm^2 * 1s * 1sr * f(1GeV, ele) ),
        length(dset.bine)-1
    )
    calc_num_events!( N, dset, f )
    return N
end

# ------------------------

# function convert_num_events_to_flux( N, dset )
#     return N ./ dset.effA_int
# end