
# const linewidth = (3+27/64)*inch
# const fullwidth = (7+5/64)*inch
# const fs = 10pt

# ========================================

function setup_flux_ax!( g )
    ax = Axis( g[1,1],
        xscale=log10,
        yscale=log10,
        xgridvisible=false,
        xticks=LogTicks(2:8),
        yticks=LogTicks(-11:-6),
        yminorticksvisible=true, yminorticks=IntervalsBetween(9),
        ylabel=L"E^2\phi\; [\textrm{GeV} \textrm{cm}^{-2} \textrm{s}^{-1} \textrm{sr}^{-1}]", xlabel=L"E [\textrm{GeV}]"
    )
    ylims!( ax, 4e-10, 1e-6 )

    return ax 
end

function plot_obs_Fpts!( ax, dsets; 
    colors=getproperty.(dsets, :color),
    labels=getproperty.(dsets, :label),
    markersize=8, 
    linewidth=1.5,
    kwargs... 
    )

    u_E = 1GeV
    u_F = 1 * GeV * (cm^2 * sr * s)^(-1)

    for (i, dset) in enumerate(dsets)

        color = colors[i]
        label = labels[i]

        bine = dset.bine
        binc = dset.binc

        islim = find_limit_pts( dset; islimit=true)
        notlim = find_limit_pts( dset; islimit=false )

        Fobs = dset.fdata
        Nobs = dset.N_obs

        F_low = copy( Fobs.low )
        F_low[ (Fobs.low .<= 0GeV/cm^2/s/sr) ] .= 1e-11 * GeV/cm^2/s/sr

        # -----------------------------

        lw = linewidth

        # plot Fobs 
            F0 = 1e-12
            scatter!( ax,
                getindex.(Ref(binc), notlim)/u_E,
                getindex.(Ref(Fobs.med), notlim)/u_F;
                color,
                markersize,
                label,
                kwargs...
            )
            rangebars!( ax,
                getindex.(Ref(Fobs.med), notlim)/u_F,
                getindex.(Ref(bine), notlim)/u_E,
                getindex.(Ref(bine), notlim .+ 1)/u_E;
                direction=:x,
                color,
                linewidth=lw,
                label,
            )
            rangebars!( ax,
                getindex.(Ref(binc), notlim)/u_E,
                getindex.(Ref(F_low), notlim)/u_F,
                getindex.(Ref(Fobs.up), notlim)/u_F;
                color,
                linewidth=lw,
            )


        # limit only points:
            rangebars!( ax,
                getindex.(Ref(Fobs.up), islim)/u_F,
                getindex.(Ref(bine), islim)/u_E,
                getindex.(Ref(bine), islim .+ 1)/u_E;
                direction=:x,
                color,
                linewidth=lw,
            )
            rangebars!( ax,
                getindex.(Ref(binc), islim)/u_E,
                [F0 for i in islim],
                getindex.(Ref(Fobs.up), islim)/u_F;
                color,
                linewidth=lw,
            )
        #
    end

end

function plot_pred_Fpts!( ax, dsets, FMOD; 
    method="mixed",
    colors=getproperty.(dsets, :color),
    show_chisq=true,
    show_chisq_pos=3e-7,
    kwargs... )

    u_E = 1GeV
    u_F = 1 * GeV * (cm^2 * sr * s)^(-1)

    tot_cost = 0.
    for (i, dset) in enumerate( dsets )

        # update calc. N_pred 
        N_pred = calc_num_events( dset, FMOD )
        F_pred = N_pred ./ dset.effA_int

        scatter!( ax, 
            dset.binc/u_E,
            F_pred/u_F;
            color=colors[i],
            kwargs...
        )

        cost_per_pt = calc_chisq_per_pt( FMOD, dset; method )

        if show_chisq
            # hack
            ixs = (F_pred/u_F .< 1e-11)
            F_pred[ixs] .= dset.fdata.med[ixs]

            text!( ax,
                dset.binc/u_E,
                # F_pred/u_F;
                fill( show_chisq_pos, length(dset.binc)),
                text = [ (@sprintf "%.1f" x) for x in cost_per_pt ],
                fontsize=8pt,
                align=(:center,:top),
                offset=(0, -10)
            )
        end
        
        tot_cost += sum( cost_per_pt )
    end

    return tot_cost
end

function plot_pred_flux!( ax, FMOD; 
    fβ=(1/3, 1/3, 1/3), kwargs... )

    u_E = 1GeV
    u_F = 1 * GeV * (cm^2 * sr * s)^(-1)
    E = exp10.(3:0.01:8) * 1GeV
    F = sum(
        FMOD.(E, β) * fβ for (β, fβ) in zip( 
            (ele, mu, tau), fβ )
    )
    # F = 1/3 * sum( FMOD.( E, β ) for β in (ele, mu, tau) )
    lines!( ax, E/u_E, F/u_F; kwargs... )

end


