
using TypedTables
using CSV

function num_saved( fname ) 
    filesize( fname ) ÷ sizeof( Float64 )
end

function setup_and_save_iter( fdir ; lb, ub, steps )
    t = Table(
        lb = collect(lb), 
        ub = collect(ub), 
        step = collect(steps) 
    )

    pts = Tuple( range(l, u; step) for (l, u, step) in zip(lb, ub, steps) )
    iter = Iterators.product( pts... )

    fname = "pts.csv"
    if !isfile( fdir * fname )
        CSV.write( fdir*fname, t )
    end

    return iter 
end

function save_calc_over_iter( 
    iter, func, fdir; 
    n_per_batch=50_000, fname="fvals.dat" )

    nsaved = 0
    if isfile( fdir * fname )
        nsaved = num_saved( fdir * fname )
    end

    println( "running calc for... $(length(iter)) pts = $(length(iter)÷n_per_batch + 1) batches")
    println( "already calculated... $(nsaved) ")

    t0 = time()
    open( fdir * fname, "a" ) do io 
        for (i, x) in enumerate( iter )
            if i > nsaved
                v = func( x )
                write( io, v )
            end

            if i % n_per_batch == 0
                str = (@sprintf "finished batch %i... in %.2f s" (i ÷ n_per_batch) time() - t0)
                println( str )
                t0 = time()
            end

        end
    end

end

function read_saved_calc( fdir; fname="fvals.dat" )

    nvals = num_saved( fdir * fname )
    fvals = zeros( nvals )

    open( fdir * fname, "r" ) do io
        for i in 1:nvals
            fvals[i] = read(io, Float64)
        end
    end

    return fvals
end

# ==============================================

function load_saved_calc_over_iter( fdir; fname="fvals.dat" )

    t = CSV.read( fdir * "pts.csv", Table )
    pts = Tuple( range( r.lb, r.ub, step=r.step ) for r in t )

    nvals = prod( length.(pts) )
    fvals = zeros( nvals )

    nsaved = num_saved( fdir * fname )
    fvals[1:nsaved] .= read_saved_calc( fdir; fname )

    fvals = reshape( fvals, Tuple(length.(pts)) )
    return pts, fvals
end