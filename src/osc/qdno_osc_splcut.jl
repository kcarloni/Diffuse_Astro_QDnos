
include( "../io_funcs.jl" )
include( "../flux/flux.jl" )
include( "qdno_osc.jl" )

using Interpolations

# ====================================================
# calc, save, load *posc*

    """
        calc_SPLcut_posc( logE, logλ, log_E1, γ0, γ1; src_zdist )

    Calculates the osc_prob ≡ calc_oscprob_ϕtot / ϕtot for 
    """
    function calc_SPLcut_posc( logE, logλ, log_Ecut, γ; src_zdist )
        src_pop = DiffuseSourcePopulation( 
            src_zdist, 
            (E -> flux_SPLcutoff(
                E; γ, 
                Ecut=exp10(log_Ecut) * 1GeV,
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
    calc_SPLcut_posc( x; src_zdist ) = calc_SPLcut_posc( x...; src_zdist )


    """
        calc_and_save_BPL_posc( fdir; lb, ub, steps, src_zdist )

    Calculates the oscillation probabilities for a 4D grid
    """
    function calc_and_save_SPLcut_posc( fdir; src_zdist, 
        lb = ( 3, -2, 3, 1 ),
        ub = ( 7, 2.5, 7, 5 ),
        steps = (0.1, 0.1, 0.5, 0.5),
        )

        func(x) = calc_SPLcut_posc(x; src_zdist)

        fdir = fdir * "SPLcut_posc/"
        mkpath( fdir )

        iter = setup_and_save_iter( fdir; lb, ub, steps )
        save_calc_over_iter(iter, func, fdir)
    end

    """
        load_saved_BPL_posc( fdir )

    Loads the interpolated flux as calculated by `calc_and_save_BPL_posc`.
    """
    function load_saved_itp_of_SPLcut_posc( fdir )

        fdir = fdir * "SPLcut_posc/"
        pts, fvals = load_saved_calc_over_iter( fdir )

        itp = linear_interpolation( pts, fvals, 
            extrapolation_bc=Interpolations.Flat() )

        λ0 = natural(1Gpc )
        f_posc( E, λ; γ, Ecut ) = itp(
            log10(E/1GeV),
            log10(λ/λ0),
            log10(Ecut/1GeV),
            γ
        )

        return f_posc
    end
#

# ====================================================
# calc, save, load *flux* 

    """
        calc_SPLcut_flux( logE, log_Ecut, γ; src_zdist )
    """
    function calc_SPLcut_flux( logE, log_Ecut, γ; src_zdist )
        src_pop = DiffuseSourcePopulation( 
            src_zdist, 
            (E -> flux_SPLcutoff(
                E; γ, 
                Ecut=exp10(log_Ecut) * 1GeV,
                # multiply flux to have correct units for source spectrum
                ϕ0=1e60 * cm^2 )
            ),
            src_zdist.xmin, 
            src_zdist.xmax
        )
        
        ϕtot( src_pop, exp10( logE ) * 1GeV )
    end
    calc_SPLcut_flux( x; src_zdist ) = calc_SPLcut_flux( x...; src_zdist )


    """
        calc_and_save_SPLcut_flux( fdir; lb, ub, steps, src_zdist )
    """
    function calc_and_save_SPLcut_flux( fdir; src_zdist, 
        lb = ( 2.5, 2.5, 1 ),
        ub = ( 7.5, 7.5, 5 ),
        steps = (0.1, 0.5, 0.5 ),
        )

        ϕ0 = 1 / GeV / cm^2 / s / sr
        func(x) = log10( calc_SPLcut_flux( x...; src_zdist ) / ϕ0 )

        fdir = fdir * "SPLcut_flux/"
        mkpath( fdir )

        iter = setup_and_save_iter( fdir; lb, ub, steps )
        save_calc_over_iter(iter, func, fdir)
    end


    """
        load_saved_itp_of_SPLcut_flux( fdir )

    Loads the interpolated flux as calculated by `calc_and_save_SPLcut_flux`.
    """
    function load_saved_itp_of_SPLcut_flux( fdir )
        
        fdir = fdir * "SPLcut_flux/"

        pts, fvals = load_saved_calc_over_iter( fdir )
        itp = linear_interpolation(
            pts, fvals,
            extrapolation_bc=Interpolations.Line()
        )

        function flux(E; γ, Ecut, ϕ0, E0=1TeV ) 
            
            C0 = 3e-18 * 1/GeV * 1/s * 1/sr * 1/cm^2
            
            ϕ_ref = exp10( itp(
                log10(E0/1GeV),
                log10(Ecut/1GeV), 
                γ
            ))

            return ϕ0 * C0 * exp10( itp(
                log10(E/1GeV),
                log10(Ecut/1GeV), 
                γ
            )) / ϕ_ref 

        end

        return flux 
    end

#