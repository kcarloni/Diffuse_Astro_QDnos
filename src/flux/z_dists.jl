
using DelimitedFiles
using Interpolations
using QuadGK
include( "cosmo.jl" )

struct NormalizedFunct{F, Q, QF}
    f::F
    xmin::Q
    xmax::Q
    f_int::QF
    f_int_err::QF
end
Base.broadcastable( F::NormalizedFunct ) = ref( F )
( F::NormalizedFunct )( x ) = F.f(x) / F.f_int

# NormalizedFunct( f, xmin, xmax ) = NormalizedFunct( 
#     f, xmin, xmax, quadgk(f, xmin, xmax)... )

# utility -- integrate SFRD over all space to get total SFR
function int_over_Vco( ρ, zmax=12 )
    return quadgk(
        (z -> ρ(z) * dV_co_dz(z)),
        0, zmax
    )
end
#

# ===========================
# star formation rate density models 

    """
        ρ_YK08( z, a, b, c, ρ0, η, B, C  )

    Star Formation Rate Density ρ̇*(z) parametrization from 
    Yusef + Kuchik 2008, arXiv:0804.4008v2 eq. (5)
    """
    ρ_YK08( z::Float64, a, b, c, ρ0, η, B, C ) = ρ0 * (
        (1 + z)^(a*η) + ( (1+z)/B )^(b*η) + ((1+z)/C)^(c*η)
    )^(1/η)
    ρ_YK08( z::Int, args... ) = ρ_YK08( float(z), args... )

    ρ_YK08(z) = ρ_YK08( z, 3.4, -0.3, -3.5,  0.02 * Msun / u"yr" / Mpc^3, -10, 5000, 9 )

    """
        ρ_EM20( z )

    Star Formation Rate Density ρ̇*(z) parametrization from 
    Elías + Martínez 2020, arXiv:2006.03367v1 eq. (8)
    """
    ρ_EM20(z) = ρ_YK08(z, 3,  -0.94, -4.36,  0.02 * Msun / u"yr" / Mpc^3, -10, 18.27, 6.61)

#

# ===========================
# generic z-distribution functions

    """
        f_PL( z, m )

    Super-generic power law scaling behavior for source evolution, as used by e.g. arXiv:2502.08508, arXiv:1406.2160
    """
    f_PL(z, m) = (1 + z)^m


    """
        f_Capel22(z, p1, p2, zc)

    A three-parameter model for cosmological source evolution from Capel+ (2022), arXiv:2005.02395. 
    """
    f_Capel22(z, p1, p2, zc) = (1 + z)^p1 / (1 + z/zc)^p2 


    f_Capel22_pos(z) = f_Capel22(z, 19.3, 24.9, 1.76)
    f_Capel22_neg(z) = f_Capel22(z, 19.3, 20.5, 1.0)


    """
        f_AhlersHalzen14(z) 

    Approximate SFR distribution as used by Ahlers + Halzen (2014), arXiv:1406.2160
    """
    function f_AhlersHalzen14(z) 

        if z ≤ 1
            return f_PL(z, 3.4)
        elseif 1 < z ≤ 4
            return f_AhlersHalzen14(1) * f_PL(z, -0.3 ) / f_PL(1, -0.3)
        else
            return f_AhlersHalzen14(4) * f_PL(z, -3.5) / f_PL(4, -3.5)
        end

    end
    #

    function f_Groth25( z; m=3, z_br=1.5 )

        if z < z_br;    return f_PL(z, m) #(1+z)^m
        elseif z < 6;   return f_PL(z_br, m) #(1+z_br)^m
        else;           return 0.
        end

    end
#

# ===========================
# load data-picked z-distributions

    zdist_path = (@__DIR__) * "/../../data/zdists/"

    function load_dist_AGNs_Ueda2014(; zdist_path=zdist_path )

        fpath = zdist_path * "agns_ueda2014_fig12/logL_42-43.csv"
        out = readdlm( fpath, ',' )

        # source number density per co-moving volume, as funct of z
        z = out[:, 1]
        n = out[:, 2] * Mpc^(-3)

        itp = linear_interpolation( z, n, extrapolation_bc=Line() )

        # normalize
        return NormalizedFunct( itp, 0, 10,
            int_over_Vco( itp, 10 )...
        )
    end

    function load_dist_FSRQs_Ajello2012(; vers="fig15",  zdist_path=zdist_path )

        # (!!) "fig15" is the version that matches "fig10" from Ajello (2014). 
        # = "number density" plot 

        if vers == "fig14" || vers == "fig15"
            fpath = zdist_path * "fsrqs_ajello2012_$(vers)/"
            out = readdlm( fpath * "fsrqs.csv", ',' )

            # source number density per co-moving volume, as funct of z
            z = out[:, 1]
            n = out[:, 2] # * 1 / Mpc^3

            itp = linear_interpolation( z, log10.(n), extrapolation_bc=Line() )
            func = (z -> exp10.(itp(z)) * u"erg" / s / Mpc^3  )

        elseif vers == "fig2"
            fpath = zdist_path * "fsrqs_ajello2012_fig2/"
            out = readdlm( fpath * "fsrqs.csv", ',' )

            # source number density per co-moving volume, as funct of z
            z = out[:, 1]
            n = out[:, 2] # * 1 / Mpc^3

            itp = linear_interpolation( z, log10.(n), extrapolation_bc=Line() )
            func_dNdz = (z -> exp10(itp(log10(z))) )
            func = (z -> func_dNdz(z) / dV_co_dz(z) )
        end

        # normalize 
        return NormalizedFunct( func, 0, 10,
            int_over_Vco( func, 10 )...
        )

    end

    function load_dist_Ajello2014(; spc="BLLacs", vers="up",  zdist_path=zdist_path )

        fpath = zdist_path * "bllacs_ajello2014_fig10/"

        if spc == "BLLacs"
            out = readdlm( fpath * "bllacs_$(vers).csv", ',' )

        elseif spc == "FSRQs"
            out = readdlm( fpath * "fsrqs_$(vers).csv", ',' )

        end

        # source number density per co-moving volume, as funct of z
        z = out[:, 1]
        n = out[:, 2] # * 1 / Mpc^3

        itp = linear_interpolation( z, log10.(n), extrapolation_bc=Line() )
        func = (z -> exp10.(itp(z)) *1 / Mpc^3  )

        # normalize 
        return NormalizedFunct( func, 0, 10,
            int_over_Vco( func, 10 )...
        )

    end

    function load_dist_RadioGalaxies_Fukazawa2020(; vers="BLLz",  zdist_path=zdist_path )

        fpath = zdist_path * "radio_fukazawa2022_fig10/"
        out = readdlm( fpath * "$(vers).csv", ',' )

        z = out[:, 1]
        dNdz = out[:, 2] 

        itp = linear_interpolation( log10.(z), log10.(dNdz), extrapolation_bc=Line() )
        func_dNdz = (z -> exp10(itp(log10(z))) )

        func_n = (z -> func_dNdz(z) / dV_co_dz(z) )

        return NormalizedFunct( func_n, 0, 10,
            int_over_Vco( func_n, 10 )...
        )
    end

    function load_dist_Groth2025(; spc="BLLacs",  zdist_path=zdist_path  )

        if spc == "BLLacs"
            fpath = zdist_path * "groth2025_fig3/BLLac.csv"
        elseif spc == "FSRQs"
            fpath = zdist_path * "groth2025_fig3/FSRQ.csv"
        elseif spc=="GC-int"
            fpath = zdist_path * "groth2025_fig3/GC-int.csv"
        elseif spc=="RL_AGNs"
            fpath = zdist_path * "groth2025_fig3/RL_AGN.csv"
        elseif spc=="RQ_AGNs"
            fpath = zdist_path * "groth2025_fig3/RQ_AGN.csv"
        elseif spc=="SBGs"
            fpath = zdist_path * "groth2025_fig3/SBG.csv"
        end

        out = readdlm( fpath, ',' )

        # source number density per co-moving volume, as funct of z
        z = out[:, 1]
        n = out[:, 2] #* Mpc^(-3)

        ixs = sortperm(z)

        itp = linear_interpolation( z[ixs], n[ixs], extrapolation_bc=Line() )

        # normalize
        return NormalizedFunct( itp, 0, 10,
            int_over_Vco( itp, 10 )...
        )  
    end

#
