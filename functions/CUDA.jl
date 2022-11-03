using CUDA 

const n = 2^20 # 1048576, number of elements in 1D arrays
const THREADS_PER_BLOCK = 256

# create a vector [0.0, 0.0, 0.0...], and send to GPU
C = zeros(Float32, n) |> cu 

# create two vectors, fill them with 1s and 2s, and send to GPU
A = fill(1.0f0, n) |> cu 
B = fill(2.0f0, n) |> cu

function add!(c, a, b)
    # compute the thread_id
    x = (blockIdx().x - 1) * blockDim().x + threadIdx().x 

    # i'th thread should get the i'th elements of a and b 
    # sum them, and then store them in i'th element of c
    @inbounds c[x] = a[x] + b[x]
    return
end

threads = THREADS_PER_BLOCK # 256
blocks = n ÷ threads # 4096

# launch the kernel with 4096 blocks, of 256 threads each
@cuda threads=threads blocks=blocks add!(C, A, B)

println(C[1])
