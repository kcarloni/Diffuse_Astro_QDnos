
using DelimitedFiles
using Interpolations
using QuasiDiracOscillations: ele, mu, tau

include("../utils.jl")

struct EffectiveArea{QF1,QF2,QF3}
    itp_ele::QF1
    itp_mu::QF2
    itp_tau::QF3
end
function (F::EffectiveArea)( E, β )
    if β == ele
        return exp10( F.itp_ele(log10(E/1GeV)) )*m^2 
    elseif β == mu
        return exp10( F.itp_mu(log10(E/1GeV))  )*m^2 
    elseif β == tau
        return exp10( F.itp_tau(log10(E/1GeV)) )*m^2 
    end
end

"""
    load_effA( dset_label )

returns an EffectiveArea functor `F` which is callable as `F( log10energy, β::Flavor )`

effective areas are available for: 
    - "Cascades 6yr 2020"
    - "NorthernTracks 9.5yr 2022"
    - "ESTES 10.3yr 2024"
    - "CombinedFit 2025"
"""
function load_effA( dset_label; data_path=(@__DIR__) * "/../../data/flux/" )

    if dset_label == "Cascades 6yr 2020"
        return load_effA_Cascade2020(; data_path)

    elseif dset_label == "NorthernTracks 9.5yr 2022"
        # tracks
        # return load_effA_CombinedFit2025(; data_path)[1]
        return load_effA_NorthernTracks2025(; data_path )

    elseif dset_label == "ESTES 10.3yr 2024"
        return load_effA_ESTES2024(; data_path)

    elseif dset_label == "CombinedFit 2025"
        # tracks, cascades
        return load_effA_CombinedFit2025(; data_path)

    end

end

function load_effA_Cascade2020(; data_path=(@__DIR__) * "/../../data/flux/")

    fpath = data_path*"2020_6yr_cascade_effA/"
    u_effA = m^2 
    u_energy = GeV

    nbinedges = 16*4+1
    log_binedges = range(3, 7, nbinedges)
    ΔlogE = step( log_binedges )
    binedges = exp10.( log_binedges )

    effA_objects = []
    flavors = ("nu_e", "nu_mu", "nu_tau")

    for flavor in flavors 

        x = readdlm( fpath * "$flavor.csv", ',' )
        ixs = sortperm( x[:,1] )
        effA = x[ixs,2]
        len_x = length( effA )

        f_itp = linear_interpolation( 
            get_centers(log_binedges[nbinedges-len_x:nbinedges]), 
            log10.(effA), 
            extrapolation_bc=Interpolations.Line() 
        )
        push!( effA_objects, f_itp )
    end

    return EffectiveArea( effA_objects... )
end

function load_effA_NorthernTracks2025(; data_path=(@__DIR__) * "/../../data/flux/")

    fpath = data_path*"2022_9.5yr_NT_effA/"
    u_effA = m^2 
    u_energy = GeV

    nbinedges = (5*6)+1  
    log_binedges = range(2, 8, nbinedges)
    ΔlogE = step( log_binedges )
    binedges = exp10.( log_binedges )

    # tracks
    x = readdlm( fpath * "effA.csv", ',' )
    effA = x[:, 2]
    len_x = length(effA)

    itp = linear_interpolation( 
        get_centers(log_binedges[(nbinedges-len_x):end]), 
        log10.(effA), 
        extrapolation_bc=Interpolations.Line() 
    )

    # an approx...
    effA_t = EffectiveArea(
        (u -> 0.0 * itp(u) ),
        (u -> 1/1.17 * itp(u)),
        (u -> 0.17/1.17 * itp(u))
    )
end

function load_effA_ESTES2024(; data_path=(@__DIR__) * "/../../data/flux/" )

    fpath = data_path*"2024_ESTES_effA/"
    u_effA = m^2 
    u_energy = GeV

    log_binedges = range(3, 6, 16)
    nbinedges = length(log_binedges)
    ΔlogE = step(log_binedges)
    binedges = exp10.( log_binedges )

    effA_objects = []
    flavors = ("nu_e", "nu_mu", "nu_tau")

    for flavor in flavors 

        x = readdlm( fpath * "$flavor.csv", ',' )
        dec_label = x[:,4]
        x_dec_m30 = x[ (" minus30" .== dec_label), 2 ]
        x_dec_m30_to_p30 = x[ (" minus30_to_plus30" .== dec_label), 2 ]
        x_dec_p30 = x[ (" plus30" .== dec_label), 2 ]

        effA =  ( (0.5 * x_dec_m30) .+ 
            (1 * x_dec_m30_to_p30) .+ 
            (0.5 * x_dec_p30) )/2
            

        len_x = length( effA )

        f_itp = linear_interpolation( 
            get_centers(log_binedges[nbinedges-len_x:nbinedges]), 
            log10.(effA), 
            extrapolation_bc=Interpolations.Line() 
            # extrapolation_bc=-Inf
        )
        push!( effA_objects, f_itp )
    end
    return EffectiveArea(effA_objects...)
end

function load_effA_CombinedFit2025(; data_path=(@__DIR__) * "/../../data/flux/")

    fpath = data_path*"2025_combinedfit_effA/"
    u_effA = m^2 
    u_energy = GeV

    nbinedges = 20*3+1  
    log_binedges = range(2, 8, nbinedges)
    ΔlogE = step( log_binedges )
    binedges = exp10.( log_binedges )

    # tracks
    x = readdlm( fpath * "tracks_avg.csv", ',' )
    effA = x[:, 2]
    len_x = length(effA)

    itp = linear_interpolation( 
        get_centers(log_binedges[(nbinedges-len_x):end]), 
        log10.(effA), 
        extrapolation_bc=Interpolations.Line() 
    )
    # an approx...
    effA_t = EffectiveArea(
        (u -> 0.0 * itp(u) ),
        (u -> 1/1.17 * itp(u)),
        (u -> 0.17/1.17 * itp(u))
    )
    
    # cascades
    x = readdlm( fpath * "casc_avg.csv", ',' )
    effA = x[:, 2]
    len_x = length(effA)

    itp = linear_interpolation( 
        get_centers(log_binedges[(nbinedges-len_x):end]), 
        log10.(effA), 
        extrapolation_bc=Interpolations.Line() 
    )
    # an approx...
    effA_c = EffectiveArea(
        (u -> 1.4/3.03 * itp(u)),
        (u -> 0.4/3.03 * itp(u)),
        (u -> 1.23/3.03 * itp(u)),
    )
    
    return (effA_t, effA_c)
end