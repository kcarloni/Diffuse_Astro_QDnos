
project_dir = (@__DIR__) * "/../"
using Pkg; Pkg.activate( project_dir )

include( project_dir * "scripts_setup/setup.jl" )

# =================================================

function main(; dstr, asimov_hyp )

    savedir = proj_dir * "saved/realizations/$(dstr)/" * get_asimov_fpath(; asimov_hyp... )

    FINJ = get_asimov_flux_hyp(; asimov_hyp... )
    dsets = get_dsets(; dstr )

    save_realizations_of_rates( 
        FINJ, dsets, savedir * "Npred.csv",
        num=100
    )

end

# =============================================

dstr = "C20+E"
asimov_hyp = (
    dstr = dstr,
    flux_str="BPL",
    mode = "fit"
)

# asimov_hyp = (
#     dstr = "CF",
#     flux_str="BPL",
#     mode = "pub"
# )

# main(; dstr="CF", asimov_hyp )
main(; dstr="C20+E", asimov_hyp )
# main(; dstr="C20", asimov_hyp )
