
begin
# setup
proj_dir = (@__DIR__) * "/../"
    using Pkg; Pkg.activate( proj_dir )

    include( proj_dir * "src/data/load_dataset.jl")
    include( proj_dir * "src/flux/flux_models.jl")
    include( proj_dir * "src/llh.jl")
    include( proj_dir * "scripts_setup/setup_zdists.jl" )

    include( proj_dir * "src/asimov.jl" )
end

# ==============================

function get_dsets(; dstr  )

    if dstr == "C20+E" || dstr == "C+E"
        dset_strs = [
            "ESTES 10.3yr 2024",
            "Cascades 6yr 2020",
        ]

    elseif dstr == "CF"
        dset_strs = [
            "CombinedFit 2025",
        ]

    elseif dstr == "E"
        dset_strs = [
            "ESTES 10.3yr 2024",
        ]

    elseif dstr == "NT"
        dset_strs = [
            "NorthernTracks 9.5yr 2022",
        ]

    elseif dstr == "C20"
        dset_strs = [
            "Cascades 6yr 2020",
        ]

    end

    return load_dataset.(dset_strs; color=:red );
end

function get_minuit_limits(; flux_str="BPL" )

    if flux_str == "BPL"
        return Dict(
            "log10_ϕ0" => (-20, 20),
            "γ0" => (0.0, 4),
            "γ1" => (1.5, 4),
            "log10_E1" => (4, 8)
        )
    elseif flux_str=="SPLcut"
        return Dict(
            "log10_ϕ0" => (-20, 20),
            "γ" => (0.0, 4),
            "log10_Ecut" => (5, 10)
        )
    else
        return Dict(
            "log10_ϕ0" => (-20, 20),
            "γ" => (0.0, 4),    
        )
    end
end

# ==============================

function get_published_fits( ; dstr, flux_str )

    if dstr == "C20"
        return Dict(
            "SPL" => [ 
                (1.66, 0.25, -0.27) # ϕ0
                (2.53, 0.07, -0.07) # γ
            ],
            "SPLcut" => [
                (1.83, 0.37, -0.31),
                (2.45, 0.09, -0.11),
                (6.4, 0.9, -0.4)
            ],
            "BPL" => [
                (1.71, 0.64, -0.28),
                (2.11, 0.29, -0.67),
                (2.75, 0.29, -0.14),
                (4.6, 0.5, -0.2)
            ]
        )[flux_str]
    #

    elseif dstr == "NT"
        return Dict(
            "SPL" => [ 
                (1.36, 0.24, -0.25) # ϕ0
                (2.37, 0.08, -0.09) # γ
            ],
            "SPLcut" => [
                (1.64, 0.39, -0.36),
                (2.0, 0.22, -0.28),
                (6.10, 0.38, -0.26)
            ],
        )[flux_str]
    #

    elseif dstr == "E"
        return Dict(
            "SPL" => [ 
                (1.68, 0.19, -0.22) # ϕ0
                (2.58, 0.10, -0.09) # γ
            ],
            "BPL" => [
                (1.7, 0.19, -0.22),
                (2.79, 0.30, -0.50),
                (2.52, 0.10, -0.09),
                (4.36, 1, 1)
            ]
        )[flux_str]
    #

    elseif dstr == "CF"
        return Dict(
            "SPL" => [ 
                (1.8, 0.13, -0.16) # ϕ0
                (2.52, 0.036, -0.038) # γ
            ],
            "SPLcut" => [
                (2.20, 0.30, -0.25),
                (2.386, 0.081, -0.09),
                (6.15, 0.37, -0.24)
            ],
            "BPL" => [
                (1.77, 0.19, -0.18),
                (1.31, 0.51, -1.30),
                (2.735, 0.067, -0.075),
                (4.39, 0.1, -0.1)
            ]
        )[flux_str]
    #

    elseif dstr == "M"
        return Dict(
            "SPL" => [ 
                (2.13, 0.18, -0.17) # ϕ0
                (2.548, 0.039, -0.041) # γ
            ],
            "SPLcut" => [
                (3.9, 1.2, -1.2),
                (2.16, 0.11, -0.16),
                (5.52, 0.39, -0.35)
            ],
            "BPL" => [
                (2.28, 0.22, -0.20),
                (1.72, 0.26, -0.35),
                (2.839, 0.11, -0.091),
                (4.524, 0.097, -0.087)
            ]
        )[flux_str]
    #

    end
end

function get_SM_model_as_published(; flux_str, dstr )

    fα = (1/3, 2/3, 0)
    fits =  get_published_fits(; dstr, flux_str )
    fits = getindex.( fits, 1 )
    fits = [ log10(fits[1]), fits[2:end]... ]

    MOD = Dict(
        "SPL" => FMod_SPL, 
        "SPLcut" => FMod_SPLcut,
        "BPL" => ( (a,b,c) -> FMod_BPL(a,b,0.1,c))
    )[flux_str]

    return MOD(
        MVector(fits...),
        100TeV,
        MVector{3}(fα)
    )
end

# ====================================

function get_fα( src_flavor_str )
    if src_flavor_str == "piondecay"; fα = (1/3, 2/3, 0)
    elseif src_flavor_str == "muondamped"; fα = (0,1,0)
    elseif src_flavor_str == "neutrondom"; fα = (1,0,0)
    end

    return collect(fα)
end

function get_SM_model(; flux_str="BPL", dsets )

    fα = (1/3, 2/3, 0)

    pdefault = Dict(
        "SPL" => (log10(1.66), 2.53),
        "BPL" => (log10(1.71), 2.79, 2.52, 4.36),
        "SPLcut" => (log10(1.68), 2.57, 9.0)
    )[flux_str]
    MOD = Dict(
        "SPL" => FMod_SPL, 
        "SPLcut" => FMod_SPLcut,
        "BPL" => ( (a,b,c) -> FMod_BPL(a,b,0.1,c))
    )[flux_str]

    FSM = MOD(
        MVector(pdefault...),
        100TeV,
        MVector{3}(fα)
    )

    res = run_minuit( FSM, dsets;
        limits=get_minuit_limits(; flux_str),
        verbose=false
    );

    return FSM 
end

function get_QD_model(; 
    flux_str="SPL", 
    src_zdist_str = "SFR_EM20", 
    src_flavor_str="piondecay"
    )

    if src_flavor_str == "piondecay"; fα = (1/3, 2/3, 0)
    elseif src_flavor_str == "muondamped"; fα = (0,1,0)
    elseif src_flavor_str == "neutrondom"; fα = (1,0,0)
    end

    sz = (flux_str == "SPL") ? zdist_dict[src_zdist_str] : src_zdist_str

    pdefault = Dict(
        "SPL" => (log10(1.66), 2.53),
        "BPL" => (log10(1.71), 2.79, 2.52, 4.36),
        "SPLcut" => (log10(1.68), 2.57, 9.0)
    )[flux_str]

    MOD = Dict(
        "SPL" => FMod_SPL_QDfixed, 
        "SPLcut" => FMod_SPLcut_QDfixed,
        "BPL" => FMod_BPL_QDfixed
    )[flux_str]

    return MOD(
        MVector(pdefault...),
        1.91e-19eV^2 * @MVector[1,1,1],
        100TeV,
        MVector{3}(fα),
        sz
    )
end

# ====================================

function get_asimov_flux_hyp(; mode, flux_str, dstr )

    if mode == "fit"
        return get_SM_model(; 
            flux_str,
            dsets = get_dsets(; dstr )
        )

    elseif mode == "pub"
        return get_SM_model_as_published(; flux_str, dstr )
    
    end

end

function get_asimov_fpath(; mode, flux_str, dstr )
    return "hyp=$(flux_str)_$(dstr)_$(mode)/"
end

# ====================================

function load_realizations(; dstr, asimov_hyp )

    parse_vstr( N ) = parse.( Int,
        split(
            strip( N, ( '[', ']' )),
            ","
        )
    )

    savedir = proj_dir * "saved/realizations/$(dstr)/" * get_asimov_fpath(; asimov_hyp... )

    t = CSV.read( savedir * "Npred.csv", Table )
    Table( NamedTuple{columnnames(t)}( map(
        c -> parse_vstr.( getproperty(t,c) ),
        columnnames(t) 
    )) )
end

