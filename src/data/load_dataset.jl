
include( "load_data.jl" )
include( "load_effA.jl" )

include("../calcs.jl")

using TypedTables


"""
    load_dataset( dlabel )

Returns a named tuple corresponding to the dataset, including ...
Implemented datasets are: 
    - "Cascades 6yr 2020"
    - "NorthernTracks 9.5yr 2022"
    - "ESTES 10.3yr 2024"
    - "CombinedFit 2025"
"""
function load_dataset( dlabel; color=nothing )

    fdata = load_flux_data( dlabel )
    effA = load_effA( dlabel )
    ΔT = load_lifetime(dlabel)
    ΔΩ = 4π * sr # ... 

    bine = fdata.binedges
    binc = fdata.bincenters
    
    # pre-calculate effA_int = ΔT ΔΩ ∫ dE E⁻² A(E),
    # st. the piecewise event counts can be quickly calculated as 
    # Nᵢ = Fᵢ × effA_intᵢ
    dset = (
        effA=effA,
        bine=bine,
        ΔT=ΔT,
        ΔΩ=ΔΩ,
        fdata = fdata[[:med, :up, :low]],
    )
    effA_int = calc_num_events( 
        dset,
        ( (E,β)->1.0 )
    )
    N_obs = calc_num_events( dset; effA_int )

    return (
        label = dlabel, 
        # -------
        ΔT = ΔT,
        ΔΩ = ΔΩ,
        bine = bine,
        binc = binc,
        # -------
        fdata = Table(fdata[[:med, :up, :low]]),
        effA_int = effA_int,
        N_obs = Table(N_obs),
        # -------
        effA = effA,
        color = color
    )
end

# for plotting... 
function find_limit_pts( dset; islimit=true )
    u_F = 1e-8 * GeV * (cm^2 * sr * s)^(-1)

    if islimit 
        return eachindex( dset.fdata.med )[ 
            dset.fdata.med .== 0 * u_F 
        ]
    else
        return eachindex( dset.fdata.med )[ 
            dset.fdata.med .> 0 * u_F 
        ] 
    end
end
    


function copy_dataset_w_change( dset, nt )

    keys_to_use = [ k for k in keys(dset) if !(k in keys(nt)) ]
    ks = tuple( keys_to_use..., keys(nt)... )

    return NamedTuple{ks}([
        dset[keys_to_use]...,
        nt...
    ])
end
