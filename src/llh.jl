
using Minuit2
using Distributions


# =======================
# per-bin likelihoods... 

"""
    m2logL_phi_parabola( f, y, yplus, yminus )

Parabolic approximation for -2logL.
"""
function m2logL_phi_parabola( f, y, yplus, yminus )
    if f < y
        chi_sq = (f - y)^2 / (y - yminus)^2
    else
        chi_sq = (f - y)^2 / (yplus - y)^2
    end

    return chi_sq
end

"""
    m2logL_logphi_parabola( f, y, yplus, yminus )

Parabolic approximation, but in log-phi. 
"""
function m2logL_logphi_parabola( f, y, yplus, yminus )

    u_f = 3e-8 * GeV / cm^2 / s / sr

    # take the log
    f = log10( f / u_f )
    y = (y == 0*u_f) ? -5 : log10( y / u_f )
    yplus = log10( yplus / u_f )
    yminus = (yminus == 0*u_f) ? -5 : log10( yminus / u_f )

    return m2logL_phi_parabola( f, y, yplus, yminus )
end

"""
    m2logL_Poisson( k, N, Nplus, Nminus )

Poisson likelihood.
"""
function m2logL_Poisson( k, N, Nplus, Nminus )

    N = round( Int, N )
    v = -2logpdf( Poisson(k), N )
    v0 = -2logpdf( Poisson(N), N )

    return (v - v0)
end

"""
"""
function m2logL_scaled_Poisson( k, N, Nplus, Nminus )

    # fallback: if 0 < N < 0.5 and 0 < Nminus,
    # then rounding down = treat as 0 will be poor. 
    if (0 < N < 0.5) 
        return m2logL_phi_parabola(
            log10(k),
            log10(N), log10(Nplus), 
            (Nminus == 0) ? -5 : log10(Nminus)
        )
    end

    N = round( Int, N )
    v = -2logpdf( Poisson(k), N )
    v0 = -2logpdf( Poisson(N), N )

    if k > N
        v1 = -2logpdf( Poisson(Nplus), N )
    else
        v1 = -2logpdf( Poisson(Nminus), N )
    end

    # fallback: handle double zeros
    if (N == 0) && (Nminus == 0) && (k <= N)
        return (v - v0)
    elseif (N == 0) && (Nplus == 0) && (k >= N)
        return (v - v0)
    end

    return (v - v0)/(v1 - v0)
end

    # convenience methods for vectors 
        function m2logL_phi_parabola( Xpred::AbstractVector, Xobs )
            return m2logL_phi_parabola.( Xpred, Xobs.med, Xobs.up, Xobs.low )
        end

        function m2logL_logphi_parabola( Xpred::AbstractVector, Xobs )
            return m2logL_logphi_parabola.( Xpred, Xobs.med, Xobs.up, Xobs.low )
        end

        function m2logL_Poisson( Xpred::AbstractVector, Xobs )
            return m2logL_Poisson.( Xpred, Xobs.med, Xobs.up, Xobs.low )
        end

        function m2logL_scaled_Poisson( Xpred::AbstractVector, Xobs )
            return m2logL_scaled_Poisson.( Xpred, Xobs.med, Xobs.up, Xobs.low )
        end
    #

# =======================
# minimizer ...

function run_minuit( FMOD::FluxModel, dsets; 
    method="mixed", limits=Dict(), verbose=false
    )

    # setup ...

        # pre-allocate ...
        N_pred = [
            zeros( length(d.binc) ) for d in dsets
        ]
        F_pred = [
            zero( d.fdata.med ) for d in dsets
        ]

    #

    function cost( p )

        tot_cost = 0.
        FMOD.p .= p

        for i in eachindex(dsets) 

            # update calc. N_pred 
            N_pred[i] .= 0.
            calc_num_events!( N_pred[i], dsets[i], FMOD )

            # update F_pred 
            F_pred[i] .= N_pred[i] ./ dsets[i].effA_int

            tot_cost += sum(
                calc_chisq_per_pt(
                    F_pred[i],
                    N_pred[i],
                    dsets[i];
                    method
                )
            )
        end
        return tot_cost
    end

    m = Minuit( 
        cost, 
        collect( FMOD.p ),
        names = param_names( FMOD )
    )
    for (name, bds) in limits; m.limits[name] = bds; end

    res = migrad!( m )
    verbose && @println res.fval

    return res
end

# =======================
# convenience

function calc_chisq_per_pt(F_pred, N_pred, dset; method="mixed" )

    F_obs = dset.fdata
    N_obs = dset.N_obs

    if method=="phi_parabola"
        return m2logL_phi_parabola( F_pred, F_obs )

    elseif method=="logphi_parabola"
        return m2logL_logphi_parabola( F_pred, F_obs )

    elseif method=="poisson"
        return m2logL_Poisson( N_pred, N_obs )

    elseif method=="scaled_poisson"
        return m2logL_scaled_Poisson( N_pred, N_obs )

    elseif method=="mixed"

        # E_switch = 46415.88833612782GeV
        E_switch = 100TeV
        # E_switch = 215443.46900318822GeV
        # E_switch = 464158.8833612782GeV

        cost = zeros( length(N_pred) )
        j = findfirst( dset.bine .>= E_switch )

        cost[1:(j-1)] .= m2logL_phi_parabola( F_pred, F_obs )[1:(j-1)] 
        
        cost[j:end] .= m2logL_scaled_Poisson( N_pred, N_obs )[j:end]
        # cost[j:end] .= m2logL_Poisson( N_pred, N_obs )[j:end]

        return cost 
    else
        error("method $(method) not implemented!")
    end
    
end

function calc_chisq_per_pt( FMOD::FluxModel, dset; method="mixed" )

    N_pred = calc_num_events( dset, FMOD )
    F_pred = N_pred ./ dset.effA_int

    return calc_chisq_per_pt(
        F_pred,
        N_pred,
        dset;
        method
    )
end

calc_chisq( FMOD::FluxModel, dset; method="mixed" ) = sum( 
    calc_chisq_per_pt( FMOD, dset; method ))
