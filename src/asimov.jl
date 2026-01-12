

function get_asimov_rate( FMOD, dset; realization=false )

    d = dset
    N_pred = calc_num_events( d, FMOD ) 

    if realization;
        N_pred = rand.( Poisson.(N_pred) )
    end

    return N_pred
end

function calculate_sys_errors( dset )
    
    stat_err=sqrt.(dset.N_obs.med)

    tot_err_up = (dset.N_obs.up .- dset.N_obs.med)
    tot_err_low = (dset.N_obs.med .- dset.N_obs.low)

    zero_floor(x) = (x >= 0) ? x : 0

    return Table(
        stat_err = stat_err,
        sys_err_up = sqrt.( 
            zero_floor.( tot_err_up.^2 .- stat_err.^2) ),
        sys_err_low = sqrt.( 
            zero_floor.(tot_err_low.^2 .- stat_err.^2) )
    )
end
        
function setup_asimov_dset( N_pred, dset; fac_stat=1, fac_sys=1 )

    # =======================

    N_errs = calculate_sys_errors( dset )

    stat_err = fac_stat * sqrt.( N_pred )
    sys_err_up = fac_sys * N_errs.sys_err_up
    sys_err_low = fac_sys * N_errs.sys_err_low

    N_up = N_pred .+ sqrt.( stat_err.^2 .+ sys_err_up.^2 )
    N_low = N_pred .- sqrt.( stat_err.^2 .+ sys_err_low.^2 )
    N_low[ N_low .< 0. ] .= 0. 

    # =======================

    F_pred = N_pred ./ dset.effA_int
    F_up =  N_up ./ dset.effA_int
    F_low = N_low ./ dset.effA_int

    # =======================

    return copy_dataset_w_change(
        dset,
        ( 
            fdata=Table(
                med=F_pred,
                up=F_up,
                low=F_low
            ),
            N_obs=Table(
                med=N_pred,
                up=N_up,
                low=N_low,
            )
        )
    )
end

function get_asimov_dset( FMOD, dset; fac_stat=1, fac_sys=1, realization=false )
    N_pred = get_asimov_rate( FMOD, dset; realization )
    setup_asimov_dset( N_pred, dset; fac_stat, fac_sys )
end

# ====== 

# save + load
function save_realizations_of_rates( FMOD, dsets, fname; num=100 )

    mkpath(dirname(fname))

    t = Table( NamedTuple{Tuple(Symbol.(getproperty.(dsets, :label)))}([
        [ get_asimov_rate( FMOD, dsets[k]; realization=true ) for i in 1:num ] for k in keys(dsets)
    ]))


    CSV.write( fname, t, append=isfile(fname) )

end