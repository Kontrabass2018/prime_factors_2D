function load_prime_factors_data(filename)
    infile = h5open(filename, "r")
    MM = infile["data"][:,:]
    close(infile)
    return MM
end 

struct NPrimes{T} P::T end
NPrimes(N::Integer = 10^8) = NPrimes(primes(N))
function _nprimes(x::T, P) where T
    x < P[1] && return zero(T)
    x == P[end] && return length(P)
    x > P[end] && return typemax(T)
    i, j = 1, length(P)
    while i != j-1
        k = i + (j - i) ÷ 2
        if P[k] > x
            j = k
        elseif P[k] ≤ x
            i = k
        end
    end
    return T(i)
end
(f::NPrimes)(x) = _nprimes(x, f.P)
nprimes = NPrimes()

function prime_factor_div(X, PF)
    return [X % div == 0 for div in PF]
end 

function load_prime_factors_data(filename)
    infile = h5open(filename, "r")
    MM = infile["data"][:,:]
    close(infile)
    return MM
end 

function get_prime_matrix(;n=100_000, outdir)

    PF = NPrimes(n).P
    d = length(PF)
    MM = zeros(n, d)

    for i in 1:n
        MM[i,:] = prime_factor_div(i, PF)
        if i % 10_000 == 0
            println("completed : $(i * 100 / n)% ")
        end 
    end 
    # save number matrix for later 
    outfile = h5open("$(outdir)/prime_divisibility_$(n).h5", "w")
    outfile["data"] = MM
    close(outfile)
    return MM, PF
end 