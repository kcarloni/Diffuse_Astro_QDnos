
proj_dir = (@__DIR__) * "/../"

using Pkg
Pkg.activate( proj_dir )

include( proj_dir * "src/flux/z_dists.jl" )
include( proj_dir * "src/osc/qdno_osc_bpl.jl" )
include( proj_dir * "src/osc/qdno_osc_splcut.jl" )

include( proj_dir * "scripts_setup/setup_zdists.jl")

# ====================================

zdist_num = parse(Int, ARGS[1])
src_zdist_str = src_zdist_strs[zdist_num]
src_zdist = zdist_dict[src_zdist_str]

println()
@show zdist_num
@show src_zdist_str
println()

# ===========================================

fdir = proj_dir * "saved/calcs/$(src_zdist_str)/"
@time calc_and_save_BPL_posc( 
    fdir; 
    src_zdist, 
    steps=(0.05, 0.02, 0.25, 0.25, 0.25)
)

fdir = proj_dir * "saved/calcs/$(src_zdist_str)/"
@time calc_and_save_BPL_flux(
    fdir;
    src_zdist,
    steps=(0.05, 0.25, 0.25, 0.25)
)

fdir = proj_dir * "saved/calcs/$(src_zdist_str)/"
@time calc_and_save_SPLcut_posc( 
    fdir; 
    src_zdist, 
    steps=(0.05, 0.05, 0.25, 0.25)
)

fdir = proj_dir * "saved/calcs/$(src_zdist_str)/"
@time calc_and_save_SPLcut_flux( 
    fdir; 
    src_zdist, 
    steps=(0.02, 0.25, 0.05)
)
