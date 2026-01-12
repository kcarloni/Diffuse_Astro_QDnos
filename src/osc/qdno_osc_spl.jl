
using Interpolations
using QuasiDiracOscillations

include( "../flux/z_dists.jl")
include( "../flux/flux.jl")
include( "qdno_osc.jl")


function calc_SPL_posc( logλ, γ; src_zdist )

    src_pop = DiffuseSourcePopulation( 
        src_zdist, 
        (E -> flux_SPL( E; γ, ϕ0=1 )),
        src_zdist.xmin, 
        src_zdist.xmax
    )

    E = 1PeV # factors out, so can just stick a constant here 
    return calc_oscprob_ϕtot( E, exp10(logλ) * 1Gpc, src_pop ) / 
        ϕtot( src_pop, E )
end
calc_SPL_posc( x; src_zdist ) = calc_SPL_posc( x...; src_zdist )

function get_itp_of_SPL_posc(; src_zdist,
    lb = ( -2, 1 ),
    ub = ( 2.5, 5 ),
    steps = ( 0.05, 0.1 )
)

    pts = Tuple( range(l, u; step) for (l, u, step) in zip(lb, ub, steps) )
    iter = Iterators.product( pts... )

    fvals = [ calc_SPL_posc( x; src_zdist ) for x in iter ]

    # integrating ConvolutionInterpolations is really slow... 
    # # interpolation ... 
    # @time citp = convolution_interpolation( pts, fvals,
    #     extrapolation_bc=ConvolutionInterpolations.Line()
    # );

    citp = linear_interpolation(
        pts,
        fvals,
        extrapolation_bc=Interpolations.Flat()
    )

    λ0 = natural(1Gpc )
    f_posc( λ; γ ) = citp( Float64.( (log10( λ / λ0 ), γ) )... )
end




# function get_itp_for_osc( ρ_SFRD )

#     pts_gamma = range(-1.0, 5.0, step=0.2)
#     pts_logλ = range(-2, 2, step=0.005)
#     pts_λ = exp10.( pts_logλ ) * 1Gpc

#     vals_posc = Matrix{Float64}(undef, length(pts_logλ), length(pts_gamma))

#     for (i, γ) in enumerate(pts_gamma)
#         zdist = get_zdist( ρ_SFRD, γ )
#         vals_posc[:, i] .= calc_oscprob_int.( zdist, pts_λ )
#     end

#     itp = linear_interpolation(
#         (pts_logλ, pts_gamma),
#         vals_posc,
#         extrapolation_bc=Line()
#     ) 
#     λ0 = natural(1Gpc)
#     function f_itp( λ, γ )
#         logλ = log10( λ / λ0 )

#         if logλ < -2; return 0.5
#         elseif logλ > 2; return 1
#         else
#             return itp( logλ, γ )
#         end
#     end

#     return f_itp
# end

# calc_oscprob_int( zdist, λ_osc; approx=true ) = quadgk(
#     z -> zdist(z) * calc_oscprob(z, λ_osc; approx), 
#     zdist.zmin, zdist.zmax 
# )[1]
