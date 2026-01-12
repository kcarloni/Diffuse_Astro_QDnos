
include( "../src/flux/z_dists.jl" )

n_SFR_EM20 = NormalizedFunct( 
    ρ_EM20, 0, 10, 
    int_over_Vco( ρ_EM20, 10 )...
)
n_SFR_YK08 = NormalizedFunct(
    ρ_YK08, 0, 10, 
    int_over_Vco( ρ_YK08, 10 )...
)

n_Capel22_pos = NormalizedFunct( 
    f_Capel22_pos, 0, 10, 
    int_over_Vco( f_Capel22_pos, 10 )...
)
n_Capel22_neg = NormalizedFunct( 
    f_Capel22_neg, 0, 10, 
    int_over_Vco( f_Capel22_neg, 10 )...
)
n_zneg = NormalizedFunct(
    ( z -> f_PL(z,-1) ), 0, 10, 
    int_over_Vco( ( z -> f_PL(z,-1) ), 10 )...
)
n_zflat = NormalizedFunct(
    ( z -> f_PL(z,0) ), 0, 10, 
    int_over_Vco( ( z -> f_PL(z,0) ), 10 )...
)

# n_AGNs_U14 = load_dist_AGNs_Ueda2014()
# # n_FSRQs_A14 = load_dist_Ajello2014(; spc="FSRQs", vers="low")
# n_FSRQS_A12 = load_dist_FSRQs_Ajello2012(; vers="fig15")
# n_BLLacs_A14 = load_dist_Ajello2014(; spc="BLLacs", vers="up")
# n_RadioGalx_F22 = load_dist_RadioGalaxies_Fukazawa2020()

# from Grath 2025
n_BLLacs_G25 = load_dist_Groth2025(; spc="BLLacs" )
n_FSRQs_G25 = load_dist_Groth2025(; spc="FSRQs" )
n_GCint_G25 = load_dist_Groth2025(; spc="GC-int" )
n_RLAGNs_G25 = load_dist_Groth2025(; spc="RL_AGNs" )
n_RQAGNs_G25 = load_dist_Groth2025(; spc="RQ_AGNs" )
n_SBG_G25 = load_dist_Groth2025(; spc="SBGs" )



# ==========================================

zdist_dict = Dict(
    "SFR_EM20" => n_SFR_EM20,
    "SFR_YK08" => n_SFR_YK08,
    "BLLacs_G25" => n_BLLacs_G25,
    "FSRQs_G25" => n_FSRQs_G25,
    "GCint_G25" => n_GCint_G25,
    "RLAGNs_G25" => n_RLAGNs_G25,
    "RQAGNs_G25" => n_RQAGNs_G25,
    "SBG_G25" => n_SBG_G25,
    "n_Capel22_pos" => n_Capel22_pos,
    "n_Capel22_neg" => n_Capel22_neg,
    "n_zflat" => n_zflat,
    "n_zneg" => n_zneg
)

const src_zdist_strs = [
    "SFR_EM20",      # 1 ✓
    "SFR_YK08",      # 2
    "BLLacs_G25",    # 3 ✓
    "FSRQs_G25",     # 4 ✓
    "GCint_G25",     # 5 
    "RLAGNs_G25",    # 6 ✓
    "RQAGNs_G25",    # 7 ✓
    "SBG_G25",       # 8
    "n_Capel22_pos", # 9 ✓
    "n_Capel22_neg", # 10 ✓
    "n_zflat",       # 11 ✓ 
    "n_zneg"         # 12 
]