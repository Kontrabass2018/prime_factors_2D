include("engines/init.jl")

using Primes
using Distances

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

# n = 1_000_000
n = 100_000
PF = NPrimes(n).P
d = length(PF)
# MM = Array{Int64, 2}(undef, (n,d))
MM = zeros(n, d)
## 
for i in 1:n
    MM[i,:] = prime_factor_div(i, PF)
    if i % 10_000 == 0
        println("completed : $(i * 100 / n)% ")
    end 
end 
outfile = h5open("Data/NUMBERS/prime_divisibility_$(n).h5", "w")
outfile["data"] = MM
close(outfile)

@time prime_umap = UMAP_(Matrix(MM'), 2;min_dist = 0.7, n_neighbors = 21);
## save for later
outfile = h5open("Data/NUMBERS/UMAP_prime_divisibility_$(n).h5", "w")
outfile["data"] = prime_umap.embedding
close(outfile)
## plot 
fig = Figure(size = (1024,1024));
ax = Axis(fig[1,1], backgroundcolor = :black);
scatter!(ax, prime_umap.embedding[1,:], prime_umap.embedding[2,:], color = collect(1:n), markersize = 1, colormap=:viridis)
hidedecorations!(ax)
hidespines!(ax)
fig
CairoMakie.save("figures/UMAP_prime_divisibility_$(n).pdf",fig)
CairoMakie.save("figures/UMAP_prime_divisibility_$(n).png",fig)






