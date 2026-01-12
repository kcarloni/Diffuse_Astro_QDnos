
get_centers( x ) = 0.5 .* (x[1:end-1] .+ x[2:end])
get_centers( x::AbstractRange ) = range( first(x)+step(x)/2, last(x)-step(x)/2, step=step(x) )