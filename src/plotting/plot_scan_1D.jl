
using DelimitedFiles
using CSV 
using TypedTables
using StatsBase

include("/Users/kiara/harvard_grad/research_projects/load_QD-no_PS.jl")

function load_QD1param_scan( fname )
    t = CSV.read( fname, Table )
    dmsq = exp10.( t.log_dmsq ) * 1eV^2 

    return Table(
        t;
        dmsq=dmsq
    )
    
end

# =====================

function setup_ΔTS_scan_plot()
    fig = Figure(
        # size = (linewidth, 3.5inch), 
        # fontsize = fs
    )
    ax = setup_ΔTS_scan_plot!( fig[1,1] )
    return fig, ax 
end

function setup_ΔTS_scan_plot!( g )

    ax = Axis( g[1,1],
        xlabel=L"\delta m^2 [\textrm{eV}^2]",
        ylabel=L"\Delta \textrm{TS}",
        # -
        xscale=log10,
        xticks=LogTicks(-21:-16),
        xminorticksvisible=true,
        xminorticks=IntervalsBetween(9),
        xgridvisible=false,
        # - 
        yminorticksvisible=true,
        yminorticks=-20:70,
        yticks=[ -25, -16, -9, -4, 0, 4, 9, 16, 25 ],
    )
    ylims!( ax, -5, 19 )
    xlims!( ax, 3e-22, 2e-16 )
    lines!( ax, [3e-22, 2e-16], fill(0, 2); 
        linewidth=1, color=:black )
        
    return ax
end

function setup_chisq_plot!( g )

    ax = Axis( g[1,1],
        xlabel=L"\delta m^2 [\textrm{eV}^2]",
        ylabel=L"\chi^2",
        # -
        xscale=log10,
        xticks=LogTicks(-21:-16),
        xminorticksvisible=true,
        xminorticks=IntervalsBetween(9),
        xgridvisible=false,
        # - 
        yminorticksvisible=true,
        yminorticks=-20:50,
    )
    xlims!( ax, 3e-22, 2e-16 )

    return ax
end

# =====================


function plot_pointsource_ΔTS_scan!( g; kwargs... )

    f_ΔTS = load_ICPS_QD_sensitivity()

    dmsq = exp10.(-22:0.01:-15.5) * 1eV^2
    lines!( ax,
        dmsq/1eV^2,
        f_ΔTS.( -1 * log10.(dmsq/1eV^2) );
        label="IC Point Sources",
        kwargs...
    )

end

function plot_SN1987A_ΔTS_scan!( ax; kwargs... )

    out = readdlm(
        "/Users/kiara/harvard_grad/research_projects/ref_material/QDconstraints/joint_fit.csv",
        ','
    )
    Δchisq = out[:,2] .- out[1,2]
    x_dmsq = out[:,1] * 1eV^2
    
    append!( x_dmsq, 1e-15*eV^2 )
    append!( Δchisq, 0 )

    lines!( ax, x_dmsq[:,1]/1eV^2, Δchisq; 
        label="SN1987A", kwargs... )

    # plot_excluded_region!(
    #     ax, x_dmsq, Δchisq; color=values(kwargs).color
    # )
end


# =====================

function plot_ΔTS_scan!( ax, fname; cost_0=nothing, kwargs... )

    t = load_QD1param_scan( fname )
    cost_pts = t.chi_sq
    dmsq_pts = exp10.( t.log_dmsq ) * 1eV^2 

    cost_0 = isnothing(cost_0) ? cost_pts[1] : cost_0

    ΔTS = cost_pts .- cost_0
    lines!( ax, 
        dmsq_pts/1eV^2, 
        ΔTS; 
        kwargs...
    )

    # plot_excluded_region!( 
    #     ax, dmsq_pts, ΔTS; 
    # )
end


function calc_ΔTS_brazil_bands( files, fname;
    cost_0_vec=nothing, 
    quantiles=[0.34, 0.66],
    
)

    dmsq_pts = exp10.( load_QD1param_scan( files[1] * fname ).log_dmsq ) * 1eV^2

    mat = zeros( length(files), length(dmsq_pts) )
    for (i, fdir) in enumerate( files )
        chi_sq = load_QD1param_scan( fdir * fname ).chi_sq
        if isnothing(cost_0_vec)
            mat[i, :] .= (chi_sq .- chi_sq[1])
        else
            mat[i, :] .= (chi_sq .- cost_0_vec[i])
        end
    end

    d = Dict([ 
        (q => quantile.( eachcol(mat), q)) for q in quantiles 
    ])
    Dict( d..., "dmsq_pts"=>dmsq_pts )
end

function plot_excluded_region!( ax, dmsq, ΔTS; nsig=3, y_val=-4, color=:black )

    ΔTS_val = nsig^2

    i_min = findfirst( 
        ΔTS[1:end-1] .< ΔTS_val  .< ΔTS[2:end] )
    i_max = findfirst( 
        ΔTS[1:end-1] .> ΔTS_val  .> ΔTS[2:end] )

    if !isnothing(i_min) && !isnothing(i_max)
        lines!( ax, 
            dmsq[ [i_min, i_max] ]/1eV^2,
            fill(y_val, 2);
            color, linewidth=2
        )
    end
end
