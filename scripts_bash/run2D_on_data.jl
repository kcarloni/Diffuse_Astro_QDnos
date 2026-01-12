
proj_dir = (@__DIR__) * "/../"
using Pkg; Pkg.activate( proj_dir )

include( proj_dir * "scripts_setup/setup.jl" )
include( proj_dir * "src/run_scan2D.jl" )

# =============================================

# on data = specify 
# - src_zdist_str
# - flux_str
# - dstr 
# - src_flavor_str = "piondecay" or "muondamped" 

function main(; src_zdist_str, dstr, flux_str, src_flavor_str, inds_equal )

    savedir = proj_dir * "results_2D/$(dstr)/$(src_zdist_str)/$(src_flavor_str)/" * "on_data/"
    fname = "fits_$(flux_str).csv"
    fname = "equal_$(inds_equal[1])$(inds_equal[2])_" * fname

    dmsq_pts = exp10.( range(-21, -16, step=0.05) ) * eV^2

    dsets = get_dsets(; dstr )

    FQD = get_QD_model(; src_zdist_str, flux_str, src_flavor_str )
    limits = get_minuit_limits(; flux_str )
    p0 = Tuple( get_SM_model(; flux_str, dsets ).p  )

    run_QD2param_scan(
        FQD, dmsq_pts, dsets, savedir;
        inds_equal,
        fname, 
        p0, limits
    )
end

# =============================================


dstr = ARGS[1]
flux_str = ARGS[2]
src_flavor_str = ARGS[3]
inds_equal = ( [2,3], [1,3], [1,2] )[parse(Int, ARGS[4])]

src_zdist_str = src_zdist_strs[1]

main(; src_zdist_str, dstr, flux_str, src_flavor_str, inds_equal )