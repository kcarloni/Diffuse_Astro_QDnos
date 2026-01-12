
This repository contains the source code for the analysis described in the work "Signatures of quasi-Dirac neutrinos in diffuse high-energy astrophysical neutrinos," by Carloni, Porto, Arguelles, Bhupal, and Jana, available at https://arxiv.org/abs/2503.19960.

**To download the necessary Julia packages**, run the following commands in Julia:
```julia
] activate .
resolve
instantiate
```

**To test a purely SM model fit**, try running the following code block in Julia:
```julia

include( "scripts_setup/setup.jl" )

function test(; dstr, flux_str )

    published_fits = get_published_fits(; dstr, flux_str )
    dsets = get_dsets(; dstr )
    FMOD = get_SM_model(; dsets, flux_str )

    println( prod(fill("=", 20))*"\n" )
    @show dstr
    @show flux_str
    println()

    for (i, pname) in enumerate( param_names(FMOD) )

        v = FMOD.p[i]
        if pname == "log10_ϕ0";
            v = exp10(v)
            pname = "ϕ0"
        end

        v0, err_up, err_low = published_fits[i]
        
        if v > v0; nsig = (v - v0)/err_up
        elseif v < v0; nsig = (v - v0)/err_low
        else nsig = 0
        end

        if nsig <= 1; color = :green
        elseif nsig >= 3; color = :red
        else color = :cyan
        end

        println( "target: $(pname) = $(v0) + $(err_up) - $(-1*err_low)" )
        println( "fit $(pname) = $(round(v, sigdigits=4))" )
        printstyled( 
            "DIFF = $(round(nsig, sigdigits=3)) σ"*"\n";
            color
        )
        println()
    end

end

test(; dstr="C20", flux_str="SPL" )
test(; dstr="C20", flux_str="SPLcut" )
test(; dstr="C20", flux_str="BPL" )

test(; dstr="NT", flux_str="SPL" )
test(; dstr="NT", flux_str="SPLcut" )

test(; dstr="E", flux_str="SPL" )
test(; dstr="E", flux_str="BPL" )

test(; dstr="CF", flux_str="SPL" )
test(; dstr="CF", flux_str="SPLcut" )
test(; dstr="CF", flux_str="BPL" )
```

**To run a series of QD model fits** on the CombinedFit (2025) diffuse flux measurements, try:
`julia scripts_bash/run_on_data.jl CF SPL 1`
- The first argument indicates the dataset chosen (either `CF`, `C20`, `E`, or `C20+E` )
- The second argument indicates the model chosen for the flux emitted by the neutrino sources (either `SPL`, `SPE`, or `BPL`)
- The third argument is an integer indicating the redshift evolution function chosen (any integer from 1 to 12 will work, see `scripts_setup/setup_zdists.jl` for more information on the available options.)

Before running fits using `SPE` or `BPL` emission models, you will need to pre-compute a series of splines. This can be done by calling `julia scripts/save_calcs_for_zdist.jl $ZNUM`, where `ZNUM` is an integer indicating the redshift evolution function chosen. The results will be saved to a folder `saved/`. 

This code is written in the programming language `Julia`; see `https://julialang.org/` for more information.