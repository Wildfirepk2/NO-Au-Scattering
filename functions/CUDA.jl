using CUDA 
using BenchmarkTools

n = 2^15
THREADS_PER_BLOCK = 256
threads = THREADS_PER_BLOCK # 256
blocks = n ÷ threads # 4096

function add!(c, a, b)
    x = (blockIdx().x - 1) * blockDim().x + threadIdx().x 
    @inbounds c[x] = a[x] + b[x]
    return
end

function cpu_add!(f, d, e)
    for i in eachindex(f)
        f[i] = d[i] + e[i]
    end
    return
end

a = [ones(3) for _ in 1:528]
b = reduce(hcat, a) # ~2 μs
c = CuArray(b) # ~14.7 μs
@btime reduce(+, b; dims=2)

# print("ADDING\n")
# @btime @cuda threads=threads blocks=blocks add!(C, A, B)
# @btime cpu_add!(F, D, E)

# @btime reduce(+, D)
# @btime reduce(+, A)