import Base:show, ^
using .Threads 

include("randomwalks.jl")
include("langevin.jl")
include("functionals.jl")

abstract type StochasticProcess end

struct Trajectory{SP<:StochasticProcess, T}
    sp::SP
    T::Real
    function Trajectory(sp::SP, T) where {SP<:StochasticProcess}
        @assert T > 0
        new{SP, T}(sp, T)
    end
end

show(io::IO, traj::Trajectory) = print(io, "$(traj.sp) 长为 $(traj.T) 的路径")

(sp::StochasticProcess)(T) = Trajectory(sp, T)

struct TelomereRW <: StochasticProcess
    args::NTuple{5, Real}
    method
    function TelomereRW(args)
        _, α, λ, __, σ = args
        @assert 0 < α < 1 || 1 < α < 2
        @assert λ > 0 && σ > 0
        new(args, telomere_randomwalk)
    end
end

function show(io::IO, l::TelomereRW)
    str = "端粒的随机游走模型 (初始长度 $(l.args[1]), α = $(l.args[2]), λ = $(l.args[3]), μ = $(l.args[4]), σ = $(l.args[5]))"
    print(io, str)
end

struct TelomereTTRW <: StochasticProcess
    args::NTuple{6, Real}
    method
    function TelomereTTRW(args)
        _, α, γ, λ, __, σ = args
        @assert 0 < α < 1 && γ > 0
        @assert  λ > 0 && σ > 0
        new(args, telomere_temporal_tempered_randomwalk)
    end
end

function show(io::IO, l::TelomereTTRW)
    str = "端粒的时间回火随机游走模型 (初始长度 $(l.args[1]), α = $(l.args[2]), γ = $(l.args[3]), λ = $(l.args[4]), μ = $(l.args[5]), σ = $(l.args[6]))"
    print(io, str)
end

struct TelomereSTRW <: StochasticProcess
    args::NTuple{7, Real}
    method
    function TelomereSTRW(args)
        _, α, β, γ, λ, __, σ = args
        @assert α > 0  && γ > 0
        @assert β > 0
        @assert λ > 0 && σ > 0
        new(args, telomere_spatial_tempered_randomwalk)
    end
end

function show(io::IO, l::TelomereSTRW)
    str = "端粒的空间回火随机游走模型 (初始长度 $(l.args[1]), α = $(l.args[2]), β = $(l.args[3])， γ = $(l.args[4]), λ = $(l.args[5]), μ = $(l.args[6]), σ = $(l.args[7]))"
    print(io, str)
end


struct TelomereLangevin <: StochasticProcess
    args::NTuple{5, Real}
    method
    function TelomereLangevin(args)
        _, α, λ, __, σ = args
        @assert 0 < α < 1 || 1 < α < 2
        @assert  λ > 0 && σ > 0
        new(args, telomere_langevin)
    end
end

function show(io::IO, l::TelomereLangevin)
    str = "端粒的 Lagevin 方程模型 (初始长度 $(l.args[1]), α = $(l.args[2]), λ = $(l.args[3]), μ = $(l.args[4]), σ = $(l.args[5]))"
    print(io, str)
end

simulate(traj::Trajectory, τ = 1e-2) = traj.sp.method(traj.T, τ, traj.sp.args...)

struct PowerTrajectory
    traj::Trajectory
    order::Int
end

^(traj::Trajectory, order::Int) = PowerTrajectory(traj, order)

function moments(traj::Trajectory, N::Int; τ=0.01, order::Int=1)
    moment = zeros(nthreads())
    @threads for _ in 1:N
        __, x = simulate(traj, τ)
        @inbounds moment[threadid()] += x[end]^order
    end
    sum(moment)/N
end

# abstract type Functional end

struct FPT
    domain::NTuple{2, Real}
    sp::StochasticProcess
    function FPT(domain, sp)
        @assert domain[1] < domain[2]
        new(domain, sp)
    end
end

show(io::IO, f::FPT) = print(io, "$(f.sp) 关于区间 $(f.domain) 的首次通过时间")

struct OccupationTime <: StochasticProcess
    domain::NTuple{2, Real}
    sp::StochasticProcess
    method
    args
    function OccupationTime(sp)
        @assert domain[1] < domain[2]
        new(domain, sp, occupationtime, sp.args)
    end
end

show(io::IO, ot::OccupationTime) = print(io, "$(ot.sp) 在区间 $(ot.domain) 内的占据时间")


struct PowerFPT 
    fpt::FPT
    order::Int
end 

^(f::FPT, order::Int) = PowerFPT(f, order)

simulate(f::FPT, τ=1e-2) = firstpassagetime(f.domain, f.sp.method, τ, f.sp.args...)

function moments(f::FPT, N::Int; τ=1e-2, order::Int=1)
    moment = zeros(nthreads())
    @threads for _ in 1:N
        x = simulate(f, τ)
        @inbounds moment[threadid()] += x^order
    end
    sum(moment)/N
end

𝔼(traj::Trajectory; N::Int=100_000, τ=0.01) = moments(traj, N; τ=τ)
𝔼(ptraj::PowerTrajectory; N::Int=100_000, τ=0.01) = moments(ptraj.traj, N; τ=τ, order=ptraj.order)
𝔼(f::FPT; τ=1e-2, N::Int=100_000) = moments(f, N; τ=τ)
𝔼(pf::PowerFPT; τ=1e-2, N::Int=100_000) = moments(pf.fpt, N; τ=τ, order=pf.order)

