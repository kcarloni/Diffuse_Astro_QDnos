
proj_dir = (@__DIR__) * "/../"
using Pkg; Pkg.activate( proj_dir )

include( proj_dir * "scripts_setup/setup.jl" )
include( proj_dir * "src/run_scan1D.jl" )

# =============================================

# on asimov = specify 
# - src_zdist_str
# - flux_str
# - dstr 
# + asimov flux hypothesis = ( dstr, flux_str, mode )

function main(; src_zdist_str, dstr, flux_str, asimov_hyp )

    savedir = proj_dir * "results/$(dstr)/$(src_zdist_str)/" * get_asimov_fpath(; asimov_hyp... )
    fname = "fits_$(flux_str).csv"

    dmsq_pts = exp10.( range(-21, -16, step=0.02) ) * eV^2

    FINJ = get_asimov_flux_hyp(; asimov_hyp... )
    dsets = get_asimov_dset.( FINJ, get_dsets(;dstr) )

    FQD = get_QD_model(; src_zdist_str, flux_str )
    limits = get_minuit_limits(; flux_str )
    p0 = Tuple( get_SM_model(; flux_str, dsets ).p  )

    run_QD1param_scan(
        FQD, dmsq_pts, dsets, savedir;
        fname, 
        p0, limits
    )

end

# =============================================

dstr = ARGS[1]
flux_str = ARGS[2]
src_zdist_str = src_zdist_strs[ parse(Int, ARGS[3]) ]
asimov_hyp = NamedTuple{(:flux_str, :dstr, :mode)}(
    split( ARGS[4], "_" )
)

main(; src_zdist_str, dstr, flux_str, asimov_hyp )