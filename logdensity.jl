using Statistics, SpecialFunctions, QuadGK, LinearAlgebra



g(u::Float64, j::Int64; zl = 0.0, zr = 1.0) =  - ( u + zl )^ j * (zr -u) 
dg(u::Float64, j::Int64; zl = 0.0, zr = 1.0) = ( u + zl )^(j-1) * ( (j+1) * u + zl - j * zr ) 
epanechnikov(x::Float64) =  (abs(x) <= 1.0) ? 0.75 * (1.0-x*x) : 0.0
gaussian(x::Float64) = exp( -0.5 * x * x) * 0.3989422804014327
triweight(x::Float64) = (abs(x)<=1.0) ? 1.09375 * (1.0-x*x)^3  : 0.0
uniform(x::Float64) = (abs(x)<=1.0) ? 0.5 : 0.0
cosinus(x::Float64) = (abs(x)<=1.0) ?   0.7853981633974483    * cos( 1.5707963267948966 * x ) : 0.0
quartic(x::Float64) = (abs(x)<=1.0) ? (0.9375 * (1.0 - x*x)^2) : 0.0
jassert(b::Bool, s::String) = b ? nothing : throw(s)
function isinvertibleproper(A::Matrix)      # clunky and slow, I hate it
    B = try inv(A)
        catch e
            return false
        end
    true
end

function logdensity(
    X::Vector{Float64},               # data
    x::Vector{Float64},               # evaluation points
    h::Float64;                       # bandwidth
    g = g,                            # g-function in paper
    dg = dg,                          # its derivative                                        
    S = 1,                            # degree 
    minx = 0.0,                       # bottom of support
    maxx = Inf,                       # top of support
    logf = true,                      # whether to provide log density function estimates
    mz::Function = epanechnikov,      # kernel used to provide estimates of f
    anal = true                       # whether to be fussy about checking arguments
    )
    
    n = length(x)
    N = length(X)
    # run some minimal sanity checks
    jassert(1 <= S <= N, "polynomial order $S must be at least one and no greater than the sample size")
    jassert(0.0 < h < Inf, "bad bandwidth $h")
    jassert( all(minx .<= X .<= maxx), "not all data points are in the support [$minx, $maxx]")
    jassert( all(minx .<= x .<= maxx), "not all evaluation points are in the support [$minx, $maxx]")
  
    X .-= minx
    x .-= minx
    maxx -= minx
    Xh = X / h
    maxxh = maxx / h
    β = Array{Float64}(undef, n, S)
    xh = x / h
    z = min.(xh, 1.0)
    zr = min.(maxxh .- xh, 1.0)
    if logf
        num = Array{Float64}(undef, n)
    end
#~     @Threads.threads for t = 1 : n                                   # use this if you want the multithreaded version
    for t = 1 : n
        v = Xh .- xh[t]
        u = v[ -z[t] .<= v .<= zr[t] ]
        ii = length(u)
        if anal
            for j = 1:S
                jassert( g( -z[t], j; zl = z[t], zr = zr[t] ) == 0.0, "function g used is not zero at lower boundary" )
                jassert( g( zr[t], j; zl = z[t], zr = zr[t] ) == 0.0, "function g used is not zero at upper boundary" )
            end
        end
        ℓ = [ mean( dg(u[i], j; zl = z[t], zr = zr[t]) for i = 1:ii ) for j = 1 : S ] / h
        R = [ mean( g( u[i], j; zl=z[t], zr = zr[t]) * u[i]^(s-1)/ factorial(s-1) for i = 1:ii ) * h^(s-1) for j = 1 : S, s = 1 : S ]
        if anal
            if !isinvertible(R) 
                throw("right hand side matrix not invertible")
            end
        end
        β[t,:] = - R \ ℓ
        if logf
            num[t] = sum( mz( u[i] ) for i = 1 : ii ) / (N*h)
        end
    end
    if !logf return β end
    denfunc(x::Float64, c::Float64) = 0.75 * ((c == 0.0) ? x - x^3 / 3.0 : exp( c * x ) * ( c*c * (1.0-x*x) + 2.0 * c * x - 2.0 ) / c^3)
    if S == 1 && mz == epanechnikov
        denom = [ denfunc( zr[t],  β[t,1] * h ) - denfunc( -z[t], β[t,1] * h ) for  t = 1 : n ]
    elseif S == 2 && mz == epanechnikov
        function denfunc(x::Float64, a::Float64, b::Float64) 
            if b==0.0 return denfunc(x,a) end
            if b > 0.0
                sqrtb = sqrt(b)
                γ = 0.5 * (a + 2.0*b*x) / sqrtb
                return 0.75 * (- 0.125 * exp(-0.25*a*a/b) * (1.7724538509055159 * (a*a - 2.0 * b * (2.0*b+1.0) ) * erfi(γ) -
                    2 * sqrtb * exp(γ*γ) * (a-2.0*b * x) ) / sqrtb^5)
            end
            b = -b
            sqrtb = sqrt(b)
            0.75 * (exp( x * (a-b*x)) * (a+2.0*b*x) *0.25 / (b*b) -
              1.7724538509055159  * (a*a + 2.0 * (1.0-2.0*b) * b) * exp( a*a *0.25/b) * erf( (b*x -0.5*a) / sqrtb ) * 0.125 / sqrtb^5)
        end
        denom = [ denfunc( zr[t],  β[t,1] * h, β[t,2] * h^2 * 0.5 ) - denfunc( -z[t], β[t,1] * h, β[t,2] * h^2 * 0.5 ) for  t = 1 : n ]
    else
        denom  = zeros(n)
        for t = 1:n
            function numint(x::Float64)
                 mz(x) * exp( sum(β[t,s] * (x * h)^s / factorial(s) for s = 1:S))
            end
            denom[t] = try quadgk(numint, -z[t], zr[t], order = 15)[1]
                       catch e  
                          throw("integration failure in denominator")
                       end
        end

    end
    ( log.( num ./ denom), β )
end
