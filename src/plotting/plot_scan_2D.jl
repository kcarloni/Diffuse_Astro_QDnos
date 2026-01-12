
using CSV 
using TypedTables
using StatsBase

# =====================

function setup_ΔTS_scan2D_plot!( g; inds_equal=[1,2] )

    xt = -21:-16
    xtl = CairoMakie.Makie.get_ticks(
        CairoMakie.Makie.LogTicks(xt),
        log10,
        CairoMakie.Makie.Automatic(),
        exp10( minimum(xt) ),
        exp10( maximum(xt) ),
    )[2]
    inds_diff = [ i for i in 1:3 if !(i in inds_equal) ]

    ie1, ie2 = inds_equal
    id = inds_diff[1]

    ax = Axis( g[1,1],
        xticks = (xt, xtl),
        yticks = (xt, xtl),
        xlabel = L"\delta m^2_{%$ie1,%$ie2}\; [\textrm{eV}^2]", 
        ylabel = L"\delta m^2_{%$id} \; [\textrm{eV}^2]",
        limits = (extrema(xt)..., extrema(xt)...),
        aspect = DataAspect()
    )

    return ax 
end

function load_QD2param_scan( fname )

    t = CSV.read( fname, Table  )

    log_dmsq = unique( t.log_dmsq_e )
    dims = ( length(log_dmsq), length(log_dmsq) ) 

    t = Table( NamedTuple{columnnames(t)}([
        reshape( getproperty(t, sym), dims ) for sym in columnnames(t)
    ]))
end