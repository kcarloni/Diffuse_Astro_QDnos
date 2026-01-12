
using DelimitedFiles

"""
    load_flux_data( dset_label )

Returns the datapoints Fⁱ of the E² scaled per-flavor flux from dataset options:
    - "Cascades 6yr 2020"
    - "NorthernTracks 9.5yr 2022"
    - "ESTES 10.3yr 2024"
    - "CombinedFit 2025"
    - "MESE 10.5yr 2025"
"""
function load_flux_data( dset_label; data_path=(@__DIR__) * "/../../data/flux/" )

    if dset_label == "Cascades 6yr 2020"
        return load_Cascades2020(; data_path)

    elseif dset_label == "ESTES 10.3yr 2024"
        return load_ESTES2024()

    elseif dset_label == "NorthernTracks 9.5yr 2022"
        return load_NorthernTracks2022()

    elseif dset_label == "CombinedFit 2025"
        return load_CombinedFit2025(; data_path)

    elseif dset_label == "MESE 10.5yr 2025"
        return load_MESE2025(; data_path)

    end
end


function load_lifetime( dset_label )

    if dset_label == "Cascades 6yr 2020"
        return 6u"yr"

    elseif dset_label == "NorthernTracks 9.5yr 2022"
        return 9.5u"yr"

    elseif dset_label == "ESTES 10.3yr 2024"
        return 10.3u"yr"

    elseif dset_label == "CombinedFit 2025"
        # tracks, cascades
        return (9.5u"yr", 11u"yr")
    end
end


"""
    load_Cascade2020()

Returns the datapoints Fⁱ of the E² scaled per-flavor flux from the 2020 6.5yr Cascade analysis.
"""
function load_Cascades2020(; data_path=(@__DIR__) * "/../../data/flux/")

    fpath = data_path*"2020_6yr_cascade/"

    # binedges = readdlm( fpath * "cascades_binedges.csv", ',' )
    binedges = exp10.( range( 4-1/3, 7, step=1/3 ) )

    binvalues = readdlm( fpath * "cascades_binvalues.csv", ',' )
    value_label = binvalues[:, 4]
    bin_mean = binvalues[ occursin.("mean", value_label), 2 ]
    bin_up = binvalues[ occursin.("up", value_label), 2 ]
    bin_low = binvalues[ occursin.("down", value_label), 2 ]

    limit_only = [6, 9, 10] #( (bin_up .- bin_low)./bin_mean .< 1e-1 )
    bin_mean[limit_only] .= 0.
    bin_low[limit_only] .= 0. 

    for i in [7]
        # correct cropped errorbars using fixed value
        bin_low[i] = 9e-11
    end

    u_energy = GeV
    u_flux = GeV * cm^(-2) * s^(-1) * sr^(-1)

    return (
        binedges=binedges * u_energy, 
        bincenters=exp10.( get_centers(log10.(binedges)) ) * u_energy,
        med = bin_mean * u_flux,
        up    = bin_up * u_flux,
        low   = bin_low * u_flux,
    )
end

"""
    load_NorthernTracks2022()
    
Returns the datapoints Fⁱ of the E² scaled per-flavor flux from the 2022 NorthernTracks analysis.
"""
function load_NorthernTracks2022()

    binedges = [ 100GeV, 15TeV, 104TeV, 721TeV, 5PeV, 100PeV ]
    bincenters = exp10.( get_centers( log10.(binedges/1GeV)) ) * 1GeV

    bin_mean = [
        0.0
        2.08
        1.14
        0.31
        0.0
    ]
    bin_err_up = [
        3.96
        1.0
        0.38
        0.22
        0.53
    ]
    bin_err_low = [
        0.0
        1.03
        0.40
        0.29
        0.0
    ]
    bin_up = bin_mean .+ bin_err_up
    bin_low = bin_mean .- bin_err_low

    u_flux = 1e-8 * GeV/cm^2/s/sr

    # return (
    #     binedges = binedges[2:5], 
    #     bincenters = bincenters[2:4],
    #     med = bin_mean[2:4] * u_flux,
    #     up  = bin_up[2:4] * u_flux,
    #     low = bin_low[2:4] * u_flux,
    # )

    return (
        binedges = binedges, 
        bincenters = bincenters,
        med = bin_mean * u_flux,
        up  = bin_up * u_flux,
        low = bin_low * u_flux,
    )
end

"""
    load_ESTES2024()
    
Returns the datapoints Fⁱ of the E² scaled per-flavor flux from the 2024 ESTES analysis.
"""
function load_ESTES2024()

    # use tabulated results from msilva's thesis 
    binedges = exp10.(3:0.5:7)
    bincenters = exp10.( get_centers(3:0.5:7) )
    bin_mean = [
        0.0
        13.3
        3.86
        2.60
        0.97
        1.02
        0.00
        0.82
    ]
    bin_err_up = [
        10.2
        3.67
        0.85
        0.43
        0.33
        0.49
        0.28
        1.04
    ]
    bin_err_low = [
        0.0
        3.67
        0.85
        0.43
        0.33
        0.49
        0.0
        0.82
    ]

    bin_up = bin_mean .+ bin_err_up
    bin_low = bin_mean .- bin_err_low

    u_flux = 1e-8 * GeV/cm^2/s/sr
    u_energy = 1GeV

    return (
        binedges = binedges * u_energy, 
        bincenters = bincenters * u_energy,
        med = bin_mean * u_flux,
        up    = bin_up * u_flux,
        low   = bin_low * u_flux,
    )
end  

"""
    load_CombinedFit2025()
    
Returns the datapoints Fⁱ of the E² scaled per-flavor flux from the 2025 CombinedFit analysis.
"""
function load_CombinedFit2025(; data_path=(@__DIR__) * "/../../data/flux/")

    fpath = data_path * "2025_combinedfit/"

    binedges = vec( readdlm( fpath*"/binedges.txt", ',' ) )

    # the table claims the binvalues are actual flux, 
    # but comparing w. the plot I'm pretty sure they are in 1e-8GeV/cm^2/s/sr

    binvalues = readdlm( fpath*"/binvalues.txt", ',' )[:, 1:3] 
    bin_mean = binvalues[:, 1]
    bin_up = bin_mean .+ binvalues[:, 3]
    bin_low = bin_mean .+ binvalues[:, 2] # negative

    u_energy = TeV
    # u_flux = 1e-18 * GeV^(-1) * cm^(-2) * s^(-1) * sr^(-1)
    u_flux = 1e-8 * GeV/cm^2/s/sr

    return (
        binedges=binedges * u_energy, 
        bincenters=exp10.( get_centers(log10.(binedges)) ) * u_energy,
        med = bin_mean * u_flux,
        up    = bin_up * u_flux,
        low   = bin_low * u_flux,
    )

end

function load_MESE2025(; data_path=(@__DIR__) * "/../../data/flux/")

    fpath = data_path * "2025_MESE/"

    binedges = vec( readdlm( fpath*"/binedges.txt", ',' ) )

    # the table claims the binvalues are actual flux, 
    # but comparing w. the plot I'm pretty sure they are in 1e-8GeV/cm^2/s/sr
    binvalues = readdlm( fpath*"/binvalues.txt", ',' )[:, 1:3] 
    bin_mean = binvalues[:, 1]
    bin_up = bin_mean .+ binvalues[:, 3]
    bin_low = bin_mean .+ binvalues[:, 2] # negative

    u_energy = 1TeV
    # u_flux = 1e-18 * GeV^(-1) * cm^(-2) * s^(-1) * sr^(-1)
    u_flux = 1e-8 * GeV/cm^2/s/sr

    return (
        binedges=binedges * u_energy, 
        bincenters=exp10.( get_centers(log10.(binedges)) ) * u_energy,
        med = bin_mean * u_flux,
        up    = bin_up * u_flux,
        low   = bin_low * u_flux,
    )
    
end
