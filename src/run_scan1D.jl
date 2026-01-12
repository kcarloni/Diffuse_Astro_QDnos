
using TypedTables
using CSV
using StatsBase

# =================

function num_saved_in_csv( fname )
    if !isfile( fname ); return 0
    else; return length( CSV.read(fname, Table) )
    end
end


function run_QD1param_scan( FQD, dmsq_pts, dsets, savedir; 
    p0, limits=Dict(), fname="scan_1D.csv", verbose=true,
    n_per_batch=10
)

    mkpath( savedir )
    num_pts = length(dmsq_pts)

    # pick up from previous save...
    idx0 = num_saved_in_csv( savedir * fname )

    if idx0 < num_pts

        println( "saving to: $(savedir * fname)" )
        
        # setup table 
        pnames = param_names(FQD)
        t = Table( NamedTuple{ Symbol.(pnames) }([ 
            zeros(num_pts) for i in eachindex(FQD.p) ]
        ))
        t = Table( t; 
            chi_sq=zeros(num_pts), 
            log_dmsq=log10.(dmsq_pts/1eV^2) 
        )

        if idx0 > 0
            t[1:idx0] .= CSV.read( savedir * fname, Table )
            for (i, name) in enumerate( Symbol.(param_names(FQD)) )
                FQD.p[i] = t[idx0][name]
            end
            FQD.δmsq .= dmsq_pts[ idx0 ]
        else
            FQD.p .= p0
            FQD.δmsq .= 0*eV^2
        end

        Δt_1pt = @elapsed run_minuit( FQD, dsets; limits, verbose=true )

        println( "num_pts = $(num_pts - idx0)" )
        println( "estimated time = $( (num_pts - idx0) * Δt_1pt )" )
        flush(stdout)

        # ~ takes ~3.5s / pt 
        @time for i in (idx0+1):num_pts

            verbose && println( "minimizing... $i" )

            FQD.δmsq .= dmsq_pts[i]
            res = run_minuit( FQD, dsets; limits, verbose )
            t.chi_sq[i] = res.fval
            for pname in pnames 
                getproperty(t, Symbol(pname))[i] = res.values[pname]
            end

            if i % n_per_batch == 0
                CSV.write( 
                    savedir * fname, 
                    t[ (i-n_per_batch+1):i ]; 
                    append=isfile(savedir*fname) )
            end
        end

        # save remainder
        CSV.write( 
            savedir * fname, 
            t[(1 + num_pts - (num_pts % n_per_batch)):end ]; 
            append=isfile(savedir*fname) 
        )
    end
    
    println("done!")
end

# =================

function run_QD1param_scan_for_all_hypotheses( FQDs, dmsq_pts, dsets, savedir; 
    p0, limits, verbose=false )

    for k in eachindex( FQDs )

        FQD = FQDs[k]        
        fname = "fits_$(k).csv"
        run_QD1param_scan( FQD, dmsq_pts, dsets, savedir; 
            p0=p0[k], limits=limits[k], fname, verbose )

    end

end

