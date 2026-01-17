""" Julia bindings for `thor`
"""

using Plots, Printf
previous_dir = pwd()
cd(@__DIR__)
thorlib = "../build/thorlib.so"

@enum IntegrationMethod direct=1 octree=2

""" 
    Function to call the full-summation Biot Savart code
"""
function bfield(
    source_pts::Matrix{<:Real}, vol::Vector{<:Real}, current_density::Matrix{<:Real}, 
    target_pts::Matrix{<:Real}, method::IntegrationMethod; nthreads::Integer=1, phi=1e-3
)
    # Need size check on input matrices
    n = size(source_pts)[1]
    m = size(target_pts)[1]
    bx = zeros(Float64, m)
    by = zeros(Float64, m)
    bz = zeros(Float64, m)

    if (size(source_pts) != size(current_density)) && (n != length(vol))
        return 0.0      # cause an error, bad practice...
    end

    xs = Base.unsafe_convert(Ptr{Float64}, source_pts[:,1])
    ys = Base.unsafe_convert(Ptr{Float64}, source_pts[:,2])
    zs = Base.unsafe_convert(Ptr{Float64}, source_pts[:,3])
    _vol = Base.unsafe_convert(Ptr{Float64}, vol) 
    jx = Base.unsafe_convert(Ptr{Float64}, current_density[:,1])
    jy = Base.unsafe_convert(Ptr{Float64}, current_density[:,2])
    jz = Base.unsafe_convert(Ptr{Float64}, current_density[:,3])
    xt = Base.unsafe_convert(Ptr{Float64}, target_pts[:,1])
    yt = Base.unsafe_convert(Ptr{Float64}, target_pts[:,2])
    zt = Base.unsafe_convert(Ptr{Float64}, target_pts[:,3])
    _bx = Base.unsafe_convert(Ptr{Float64}, bx)
    _by = Base.unsafe_convert(Ptr{Float64}, by)
    _bz = Base.unsafe_convert(Ptr{Float64}, bz)
    nthreads = convert(Int32, nthreads)

    if method == direct
        check = @ccall thorlib.bfield_direct(
            xs::Ptr{Float64}, ys::Ptr{Float64}, zs::Ptr{Float64}, 
            vol::Ptr{Float64}, 
            jx::Ptr{Float64}, jy::Ptr{Float64}, jz::Ptr{Float64}, 
            n::Int64, 
            xt::Ptr{Float64}, yt::Ptr{Float64}, zt::Ptr{Float64},
            m::Int64, 
            _bx::Ptr{Float64}, _by::Ptr{Float64}, _bz::Ptr{Float64}, 
            nthreads::Int32
        )::Cint
    elseif method == octree
        check = @ccall thorlib.bfield_octree(
            xs::Ptr{Float64}, ys::Ptr{Float64}, zs::Ptr{Float64}, 
            _vol::Ptr{Float64}, 
            jx::Ptr{Float64}, jy::Ptr{Float64}, jz::Ptr{Float64}, 
            n::Int64, 
            xt::Ptr{Float64}, yt::Ptr{Float64}, zt::Ptr{Float64},
            m::Int64, 
            _bx::Ptr{Float64}, _by::Ptr{Float64}, _bz::Ptr{Float64}, 
            nthreads::Int32, phi::Float64
        )::Cint
    end

    return hcat(bx, by, bz)
end


"""
    Calculate the field along the axis of a current-carrying loop using the 
    exact (analytical) solution.
"""
function bfield_loop_axis(I::Float64, R::Float64, z::Vector{<:Real})

    N = length(z) 
    Bz = zeros(N) 
    _z = Base.unsafe_convert(Ptr{Float64}, z)
    _Bz = Base.unsafe_convert(Ptr{Float64}, Bz) 

    @ccall thorlib.bfield_loop_axis(_z::Ptr{Float64}, N::Int64, I::Cdouble, R::Cdouble, _Bz::Ptr{Float64})::Cint
    return Bz

end

cd(previous_dir)