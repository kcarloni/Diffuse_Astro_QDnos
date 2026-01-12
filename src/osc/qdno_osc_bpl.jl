
include( "../io_funcs.jl" )
include( "../flux/flux.jl" )
include( "qdno_osc.jl" )

using Interpolations

# ====================================================
# calc, save, load BPL *osc_prob*

    """
        calc_BPL_posc( logE, logλ, log_E1, γ0, γ1; src_zdist )

    Calculates the osc_prob ≡ calc_oscprob_ϕtot / ϕtot for energy `logE`, oscillation length `logλ`, and BPL parameters `log_E1, γ0, γ1`.
    """
    function calc_BPL_posc( logE, logλ, log_E1, γ0, γ1; src_zdist )
        src_pop = DiffuseSourcePopulation( 
            src_zdist, 
            (E -> flux_BPL(
                E; γ0, γ1, 
                E1=exp10(log_E1) * 1GeV, 
                ϕ0=1 )
            ),
            src_zdist.xmin, 
            src_zdist.xmax
        )

        calc_oscprob_ϕtot(
            exp10( logE ) * 1GeV,
            exp10( logλ ) * 1Gpc,
            src_pop
        ) / ϕtot( src_pop, exp10( logE ) * 1GeV )
    end
    calc_BPL_posc( x; src_zdist ) = calc_BPL_posc( x...; src_zdist )

    function calc_and_save_BPL_posc( fdir; src_zdist, 
        lb = ( 3, -2, 3, 1, 1 ),
        ub = ( 7, 2.5, 7, 5, 5),
        steps = (0.1, 0.1, 0.5, 0.5, 0.5),
        )

        fdir = fdir * "BPL_posc/"
        mkpath( fdir )

        iter = setup_and_save_iter( fdir; lb, ub, steps )
        func(x) = calc_BPL_posc( x...; src_zdist )
        save_calc_over_iter(iter, func, fdir)
    end

    function load_saved_itp_of_BPL_posc( fdir )

        fdir = fdir * "BPL_posc/"
        pts, fvals = load_saved_calc_over_iter( fdir )

        itp = linear_interpolation( pts, fvals, 
            extrapolation_bc=Interpolations.Flat() )

        λ0 = natural(1Gpc )
        f_posc( E, λ; γ0, γ1, E1 ) = itp(
            log10(E/1GeV),
            log10(λ/λ0),
            log10(E1/1GeV),
            γ0, γ1
        )

        return f_posc
    end

#

# ====================================================
# calc, save, load BPL flux

    """
        calc_BPL_flux( logE, logλ, log_E1, γ0, γ1; src_zdist )
    """
    function calc_BPL_flux( logE, log_E1, γ0, γ1; src_zdist )
        src_pop = DiffuseSourcePopulation( 
            src_zdist, 
            (E -> flux_BPL(
                E; γ0, γ1, 
                E1=exp10(log_E1) * 1GeV, 
                # multiply flux to have correct units for source spectrum
                ϕ0=1e60 * cm^2 )
            ),
            src_zdist.xmin, 
            src_zdist.xmax
        )
        
        ϕtot( src_pop, exp10( logE ) * 1GeV )
    end
    calc_BPL_flux( x; src_zdist ) = calc_BPL_flux( x...; src_zdist )


    """
        calc_and_save_BPL_flux( fdir; lb, ub, steps, src_zdist )
    """
    function calc_and_save_BPL_flux( fdir; src_zdist, 
        lb = (2.5, 2.5, 1, 1),
        ub = (7.5, 7.5, 5, 5),
        steps = (0.1, 0.5, 0.5, 0.5),
        )

        ϕ0 = 1 / GeV / cm^2 / s / sr
        func(x) = log10( calc_BPL_flux( x...; src_zdist ) / ϕ0 )

        fdir = fdir * "BPL_flux/"
        mkpath( fdir )

        iter = setup_and_save_iter( fdir; lb, ub, steps )
        save_calc_over_iter(iter, func, fdir)
    end

    """
        load_saved_itp_of_BPL_flux( fdir )

    Loads the interpolated flux as calculated by `calc_and_save_BPL_flux`.
    """
    function load_saved_itp_of_BPL_flux( fdir )

        fdir = fdir * "BPL_flux/"
        pts, fvals = load_saved_calc_over_iter( fdir )

        itp = linear_interpolation(
            pts, fvals,
            extrapolation_bc=Interpolations.Line()
        )
        function flux(E; γ1, γ0, E1, ϕ0, E0=1TeV ) 
            
            C0 = 3e-18 * 1/GeV * 1/s * 1/sr * 1/cm^2
            
            ϕ_ref = exp10( itp(
                log10(E0/1GeV),
                log10(E1/1GeV), 
                γ0, γ1
            ))

            return ϕ0 * C0 * exp10( itp(
                log10(E/1GeV),
                log10(E1/1GeV), 
                γ0, γ1
            )) / ϕ_ref 

        end

        return flux 
    end

#