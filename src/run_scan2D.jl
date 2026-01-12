
function num_saved_scan2D( fname )
    if !isfile( fname ); return 0
    else; return length( CSV.read(fname, Table) )
    end
end

function run_QD2param_scan( FQD, dmsq_pts, dsets, savedir;
    inds_equal=[1,2],
    p0, limits=Dict(), fname="scan_2D.csv", verbose=true,
    n_per_batch=5
)

    mkpath( savedir )
    num_pts = length(dmsq_pts)^2

    # pick up from previous save...
    idx0 = num_saved_scan2D( savedir * fname )

    if idx0 < num_pts

        inds_diff = [ i for i in 1:3 if !(i in inds_equal) ]

        function update_dmsq!( FQD, ldm_e, ldm_d )
            FQD.δmsq[inds_equal] .= exp10.(ldm_e)*1eV^2
            FQD.δmsq[inds_diff] .= exp10.(ldm_d)*1eV^2
        end

        function update_result!( t, idx, res )

            t.chi_sq[idx] = res.fval
            for pname in pnames 
                getproperty(t, Symbol(pname))[idx] = res.values[pname]
            end

        end

        log_dmsq = log10.( dmsq_pts/1eV^2 )
        iter = Iterators.product( log_dmsq, log_dmsq )

        # set up table 
        pnames = param_names(FQD)
        t = Table( NamedTuple{ Symbol.(pnames) }([ 
            zeros(num_pts) for i in eachindex(FQD.p) ]
        ))
        t = Table( t; 
            chi_sq=zeros(num_pts), 
            log_dmsq_e = vec( getindex.(iter, 1) ),
            log_dmsq_d = vec( getindex.(iter, 2) )
        )

        if idx0 > 0;
            t[1:idx0] .= CSV.read( savedir * fname, Table )
            
            for (i, name) in enumerate( Symbol.(param_names(FQD)) )
                FQD.p[i] = t[idx0][name]
            end
            update_dmsq!( FQD, t.log_dmsq_e[idx0], t.log_dmsq_d[idx0])
        else
            FQD.p .= p0
            update_dmsq!( FQD, t.log_dmsq_e[1], t.log_dmsq_d[1] )
        end

        @println Tuple( FQD.p )

        Δt_1pt = @elapsed run_minuit( 
            FQD, dsets; limits, verbose=true )
        @println "num_pts = $(num_pts - idx0)"
        @println "estimated time = $( (num_pts - idx0) * Δt_1pt )"
        flush(stdout)

        # ------

        @time for i in (idx0+1):length(t)

            verbose && @println "minimizing... $i"

            # reset, for speed
            if (t.log_dmsq_e[i] <= -21) || (t.log_dmsq_d[i] <= -21)
                FQD.p .= p0
            end

            update_dmsq!( FQD, t.log_dmsq_e[i], t.log_dmsq_d[i] )
            res = run_minuit( FQD, dsets; limits, verbose )
            update_result!( t, i, res )

            if i % n_per_batch == 0
                CSV.write( 
                    savedir * fname, 
                    t[ (i-n_per_batch+1):i ]; 
                    append=isfile(savedir*fname) )
            end

        end

        CSV.write( 
            savedir * fname, 
            t[(1 + num_pts - (num_pts % n_per_batch)):end ]; 
            append=isfile(savedir*fname) 
        )

    end
    println("done!")
end