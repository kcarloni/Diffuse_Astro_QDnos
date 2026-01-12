
This repository contains the source code for the analysis described in the work "Signatures of quasi-Dirac neutrinos in diffuse high-energy astrophysical neutrinos," by Carloni, Porto, Arguelles, Bhupal, and Jana, available at https://arxiv.org/abs/2503.19960.

To run a series of QD model fits on the CombinedFit (2025) diffuse flux measurements, try:
`julia scripts_bash/run_on_data.jl CF SPL 1`
- The first argument indicates the dataset chosen (either `CF`, `C20`, `E`, or `C20+E` )
- The second argument indicates the model chosen for the flux emitted by the neutrino sources (either `SPL`, `SPE`, or `BPL`)
- The third argument is an integer indicating the redshift evolution function chosen (any integer from 1 to 12 will work, see `scripts_setup/setup_zdists.jl` for more information on the available options.)

Before running fits using `SPE` or `BPL` emission models, you will need to pre-compute a series of splines. This can be done by calling `julia scripts/save_calcs_for_zdist.jl $ZNUM`, where `ZNUM` is an integer indicating the redshift evolution function chosen.

This code is written in the programming language `Julia`; see `https://julialang.org/` for more information.