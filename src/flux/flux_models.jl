
using QuasiDiracOscillations
using StaticArrays
using Printf

include("flux.jl")

include("../osc/qdno_osc.jl")
include("../osc/qdno_osc_spl.jl")
include("../osc/qdno_osc_bpl.jl")
include("../osc/qdno_osc_splcut.jl")

# --------

const U_PMNS = QuasiDiracOscillations.get_U_PMNS( 
    QuasiDiracOscillations.PMNS_args_NuFit_noSuperK_NO )

function calc_Posc_given_initial_fraction( β, fα )
    Posc = 0.
    for α in (ele, mu, tau), j in 1:3
        Posc += abs2( U_PMNS[α,j] ) * abs2( U_PMNS[β,j] ) * fα[α]
    end   
    return Posc
end

# ============================

"""
all FluxModels 
    - have property `p` containing parameters to be optimized over
    - are callable `(F::FluxModel)( E, β )` and return F = E²ϕ
"""
abstract type FluxModel end

function Base.show( 
    io::IO, ::MIME"text/plain", FMOD::FluxModel )
    println( io, "$(typeof(FMOD).name.wrapper) with optimizable parameters:" )
    for i in eachindex( FMOD.p )
        str = ( @sprintf "  %s = %.3f" param_names(FMOD)[i] FMOD.p[i] )
        println( io, str )
    end
    println( io, "and other parameters ", [sym for sym in propertynames(FMOD) if sym != :p] )
end

function sprint_params( FMOD::FluxModel )

    str = "$(typeof(FMOD).name.wrapper):\n"
    for i in eachindex( FMOD.p )

        if occursin("ϕ0", param_names(FMOD)[i])
            str *= ( @sprintf "  ϕ0 = %.3f" exp10(FMOD.p[i]) )
            str *= "\n"

        else
            str *= ( @sprintf "  %s = %.3f" param_names(FMOD)[i] FMOD.p[i] )
            str *= "\n"
        end
    end
    return str

end

Base.broadcastable( F::FluxModel ) = Ref(F)

# ===========================
# SPL flux model 

    """
    FluxModel implementing a single power law, with parameters: 
        p[1] = log10( ϕ0/[ 1e-18 (GeV * s * sr * cm²)⁻¹ ] )
        p[2] = γ

    this model is flavor-independent, F(E,β) = F(E).
    """
    struct FMod_SPL{QE} <: FluxModel
        p::MVector{2, Float64}
        # ~ fixed over optimization
        E0::QE
        # for strictly correct flavor-weighting 
        fα::MVector{3, Float64}
    end

    function ( F::FMod_SPL )( E, β::Flavor ) 
        # Posc = calc_Posc_given_initial_fraction( β, F.fα )
        Posc = 1/3
        ϕ0 = exp10(F.p[1]) 
        γ = F.p[2]
        return Posc * E^2 * flux_SPL( E; γ, ϕ0, E0=F.E0)
    end

    param_names( F::FMod_SPL ) = ("log10_ϕ0", "γ")
#

# ===========================
# BPL flux model

    """
    FluxModel implementing a broken power law, with parameters: 
        p[1] = log10( ϕ0/[ 1e-18 (GeV * s * sr * cm²)⁻¹ ] )
        p[2] = γ1
        p[3] = γ2
        p[4] = log10( E_break/1GeV )
    this model is flavor-independent, F(E,β) = F(E).
    """
    struct FMod_BPL{QE} <: FluxModel
        p::MVector{4, Float64}
        # ~ fixed over optimization
        E0::QE
        s::Float64
        # for strictly correct flavor-weighting (?)
        fα::MVector{3, Float64}
    end

    function ( F::FMod_BPL )( E, β::Flavor ) 

        Posc = calc_Posc_given_initial_fraction( β, F.fα )

        ϕ0 = exp10(F.p[1]) 
        E1 = exp10(F.p[4]) * 1GeV 
        γ0, γ1 = F.p[2:3]
        return Posc * E^2 * flux_BPL( 
            E; γ0, γ1, E1, ϕ0, sm=F.s, E0=F.E0 
        )

    end

    param_names( F::FMod_BPL ) = ("log10_ϕ0", "γ0", "γ1", "log10_E1")
# 

# ===========================
# SPL-cutoff flux model 

    """
    FluxModel implementing a single power law, with parameters: 
        p[1] = log10( ϕ0/[ 1e-18 (GeV * s * sr * cm²)⁻¹ ] )
        p[2] = γ
        p[3] = log10( Ecut / 1GeV )

    this model is flavor-independent, F(E,β) = F(E).
    """
    struct FMod_SPLcut{QE} <: FluxModel
        p::MVector{3, Float64}
        # ~ fixed over optimization
        E0::QE
        # for strictly correct flavor-weighting 
        fα::MVector{3, Float64}
    end
    param_names( F::FMod_SPLcut ) = (
        "log10_ϕ0", "γ", "log10_Ecut")

    function ( F::FMod_SPLcut )( E, β::Flavor ) 
        Posc = calc_Posc_given_initial_fraction( β, F.fα )
        ϕ0 = exp10(F.p[1]) #* 1e-18 * (GeV * s * sr * cm^2)^(-1)
        γ = F.p[2]
        Ecut = exp10(F.p[3]) * 1GeV 

        return Posc * E^2 * flux_SPLcutoff( E; γ, Ecut, ϕ0, E0=F.E0)
    end

#


# ===========================
# SPL QD flux model  

    """
    FluxModel implementing a SPL × Quasi-Dirac oscillations,
    with fit parameters: 
        `p[1]` = log10( ϕ0/[ 1e-18 (GeV * s * sr * cm²)⁻¹ ] )
        `p[2]` = γ
    and fixed:
        - mass-squared splitting `δmsq`
        - SPL normalization energy `E0`
        - initial flavor ratio vector fα
    """
    struct FMod_SPL_QDfixed{ QE, QESQ, F } <: FluxModel
        p::MVector{2, Float64}
        # ~ fixed over optimization
        δmsq::MVector{3, QESQ}
        E0::QE
        fα::MVector{3, Float64}
        f_posc::F

        function FMod_SPL_QDfixed( p, dmsq, E0, fα, src_zdist )
            f_itp = get_itp_of_SPL_posc(; src_zdist )
            new{
                typeof(E0),
                eltype(dmsq),
                typeof(f_itp),
            }( p, dmsq, E0, fα, f_itp )
        end
    end
    param_names( F::FMod_SPL_QDfixed ) = ("log10_ϕ0", "γ")

    function ( F::FMod_SPL_QDfixed )( E, β ) 

        γ  = F.p[2]
        ϕ0 = exp10(F.p[1])

        Posc = 0.
        for α in (ele, mu, tau), j in 1:3
            Posc += abs2( U_PMNS[α,j] ) * 
                abs2( U_PMNS[β,j] ) * 
                F.f_posc( calc_λ(E, F.δmsq[j]); γ ) * F.fα[α]
        end   
        return Posc * E^2 * flux_SPL( E; γ, ϕ0, E0=F.E0 )
    end

#

# ===========================
# fittable 1-param SPL QD flux model 

    """
    FluxModel implementing a SPL × Quasi-Dirac oscillations,
    with fit parameters: 
        `p[1]` = log10( ϕ0/[ 1e-18 (GeV * s * sr * cm²)⁻¹ ] )
        `p[2]` = γ
        `p[3]` = log10( δmsq / 1eV² )
    and fixed:
        - SPL normalization energy `E0`
        - initial flavor ratio vector fα
    """
    struct FMod_SPL_QD1Param{ QE, F } <: FluxModel
        p::MVector{3, Float64}
        # ~ fixed over optimization
        E0::QE
        fα::MVector{3, Float64}
        f_posc::F

        function FMod_SPL_QD1Param( p, E0, fα, src_zdist )
            f_itp = get_itp_of_SPL_posc(; src_zdist )
            new{
                typeof(E0),
                typeof(f_itp),
            }( p, E0, fα, f_itp )
        end
    end

    param_names( F::FMod_SPL_QD1Param ) = (
        # "log10_ϕ0", "γ", "log10_δmsq"
        "log10_ϕ0", "log10_δmsq", "γ"
    )

    function ( F::FMod_SPL_QD1Param )( E, β ) 

        ϕ0 = exp10(F.p[1])
        γ = F.p[3]
        δmsq = exp10.(F.p[2]) * 1eV^2

        Posc = 0.
        for α in (ele, mu, tau), j in 1:3
            Posc += abs2( U_PMNS[α,j] ) * 
                abs2( U_PMNS[β,j] ) * 
                F.f_posc( calc_λ(E, δmsq); γ ) * F.fα[α]
        end   

        return Posc * E^2 * flux_SPL( E; γ, ϕ0, E0=F.E0)
    end
#


# ============================
# BPL QD flux model 

    """
    FluxModel implementing a BPL × Quasi-Dirac oscillations,
    with fit parameters: 

    and fixed:
        - mass-squared splitting `δmsq`
        - initial flavor ratio vector fα
        - redshift distribution Dz::RedshiftDistribution
    """
    struct FMod_BPL_QDfixed{ QE, QESQ, F1, F2 } <: FluxModel
        p::MVector{4, Float64}
        # ~ fixed over optimization
        δmsq::MVector{3, QESQ}
        E0::QE
        fα::MVector{3, Float64}
        #
        f_flux::F1
        f_posc::F2

        function FMod_BPL_QDfixed( p, dmsq, E0, fα, src_zdist_str )

            savedir = (@__DIR__) * "/../../saved/calcs/"

            println( "loading flux..." )
            f_flux = load_saved_itp_of_BPL_flux(
                savedir * src_zdist_str * "/" )

            println( "loading posc..." )
            f_posc = load_saved_itp_of_BPL_posc( 
                savedir * src_zdist_str * "/" )

            new{
                typeof(E0),
                eltype(dmsq),
                typeof(f_flux),
                typeof(f_posc)
            }( p, dmsq, E0, fα, f_flux, f_posc )
        end
    end
    param_names( F::FMod_BPL_QDfixed ) = ("log10_ϕ0", "γ0", "γ1", "log10_E1")

    function ( F::FMod_BPL_QDfixed )( E, β ) 

        ϕ0 = exp10(F.p[1])
        E1 = exp10(F.p[4]) * 1GeV 
        γ0, γ1 = F.p[2:3]

        Posc = 0.
        for α in (ele, mu, tau), j in 1:3
            Posc += abs2( U_PMNS[α,j] ) * 
                abs2( U_PMNS[β,j] ) * 
                F.fα[α] * 
                F.f_posc( 
                    E, calc_λ(E, F.δmsq[j]);
                    E1, γ0, γ1,
                ) 
        end   
        return Posc * E^2 * F.f_flux( E; γ0, γ1, E1, ϕ0, E0=F.E0 )
    end

#

# ============================
# SPLcutoff QD flux model 

    """
    FluxModel implementing a BPL × Quasi-Dirac oscillations,
    with fit parameters: 

    and fixed:
        - mass-squared splitting `δmsq`
        - initial flavor ratio vector fα
        - redshift distribution Dz::RedshiftDistribution
    """
    struct FMod_SPLcut_QDfixed{ QE, QESQ, F1, F2 } <: FluxModel
        p::MVector{3, Float64}
        # ~ fixed over optimization
        δmsq::MVector{3, QESQ}
        E0::QE
        fα::MVector{3, Float64}
        #
        f_flux::F1
        f_posc::F2

        function FMod_SPLcut_QDfixed( p, dmsq, E0, fα, src_zdist_str )

            savedir = (@__DIR__) * "/../../saved/calcs/"

            println( "loading flux..." )
            f_flux = load_saved_itp_of_SPLcut_flux(
                savedir * src_zdist_str * "/" )

            println( "loading posc..." )
            f_posc = load_saved_itp_of_SPLcut_posc( 
                savedir * src_zdist_str * "/" )

            new{
                typeof(E0),
                eltype(dmsq),
                typeof(f_flux),
                typeof(f_posc)
            }( p, dmsq, E0, fα, f_flux, f_posc )
        end
    end
    param_names( F::FMod_SPLcut_QDfixed ) = (
        "log10_ϕ0", "γ","log10_Ecut")

    function ( F::FMod_SPLcut_QDfixed )( E, β ) 

        ϕ0 = exp10(F.p[1])
        γ = F.p[2]
        Ecut = exp10(F.p[3]) * 1GeV 

        Posc = 0.
        for α in (ele, mu, tau), j in 1:3
            Posc += abs2( U_PMNS[α,j] ) * 
                abs2( U_PMNS[β,j] ) * 
                F.fα[α] * 
                F.f_posc( 
                    E, calc_λ(E, F.δmsq[j]);
                    Ecut, γ
                ) 
        end   
        return Posc * E^2 * F.f_flux( E; Ecut, γ, ϕ0, E0=F.E0 )
    end

#