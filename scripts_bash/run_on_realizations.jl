
proj_dir = (@__DIR__) * "/../"
using Pkg; Pkg.activate( proj_dir )

include( proj_dir * "scripts_setup/setup.jl" )
include( proj_dir * "src/run_scan1D.jl" )

# =============================================

# on realization = specify 
# - src_zdist_str
# - flux_str
# - dstr 
# + asimov_hyp
# + realization num

function main( num ; src_zdist_str, dstr, flux_str, asimov_hyp )

    savedir = proj_dir * "results/$(dstr)/$(src_zdist_str)/" * get_asimov_fpath(; asimov_hyp... ) * "run_$(num)/"
    fname = "fits_$(flux_str).csv"

    dmsq_pts = exp10.( range(-21, -16, step=0.02) ) * eV^2

    N_pred = collect( load_realizations(; asimov_hyp, dstr )[num] )
    dsets = setup_asimov_dset.( N_pred, get_dsets(; dstr) )

    FQD = get_QD_model(; src_zdist_str, flux_str )
    limits = get_minuit_limits(; flux_str )
    p0 = Tuple( get_SM_model(; flux_str, dsets ).p  )

    run_QD1param_scan(
        FQD, dmsq_pts, dsets, savedir;
        fname, 
        p0, limits,
        verbose=false
    )

end

# =============================================

dstr = ARGS[1]
flux_str = ARGS[2]
src_zdist_str = src_zdist_strs[ parse(Int, ARGS[3]) ]
asimov_hyp = NamedTuple{(:flux_str, :dstr, :mode)}(
    split( ARGS[4], "_" )
)
num1, num2 = parse.(Int, (ARGS[5], ARGS[6]) )

for num in num1:num2
    main(num; src_zdist_str, dstr, flux_str, asimov_hyp )
end
