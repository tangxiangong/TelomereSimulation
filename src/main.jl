import Base:show, ^
using .Threads 
include("randoms.jl")
include("randomwalks.jl")
include("subordinators.jl")
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

abstract type Functional end

struct FunctionalPower 
    functional::Functional
    order::Int
end

struct FPT <: Functional
    domain::NTuple{2, Real}
    sp::StochasticProcess
    function FPT(domain, sp)
        @assert domain[1] < domain[2]
        new(domain, sp)
    end
end

show(io::IO, f::FPT) = print(io, "$(f.sp) 关于区间 $(f.domain) 的首次通过时间")

^(functional::Functional, order::Int) = FunctionalPower(functional, order)

simulate(f::FPT, τ=1e-2) = firstpassagetime(f.domain, f.sp.method, τ, f.sp.args...)

function moments(functional::Functional, N::Int; τ=1e-2, order::Int=1)
    moment = zeros(nthreads())
    @threads for _ in 1:N
        x = simulate(functional, τ)
        @inbounds moment[threadid()] += x^order
    end
    sum(moment)/N
end

struct OccupationTime <: Functional
    T::Real
    domain::NTuple{2, Real}
    sp::StochasticProcess
    function OccupationTime(T, domain, sp)
        @assert domain[1] < domain[2] && T > 0
        new(T, domain, sp)
    end
end

show(io::IO, ot::OccupationTime) = print(io, "$(ot.sp) 在 [0, $(ot.T)] 内逗留在 $(ot.domain) 的时间")

simulate(oc::OccupationTime, τ=1e-2) = occupationtime(oc.domain, oc.sp.method, oc.T, τ, oc.sp.args...)


𝔼(traj::Trajectory; N::Int=100_000, τ=0.01) = moments(traj, N; τ=τ)
𝔼(ptraj::PowerTrajectory; N::Int=100_000, τ=0.01) = moments(ptraj.traj, N; τ=τ, order=ptraj.order)
𝔼(functional::Functional; τ=1e-2, N::Int=100_000) = moments(functional, N; τ=τ)
𝔼(fp::FunctionalPower; τ=1e-2, N::Int=100_000) = moments(fp.functional, N; τ=τ, order=fp.order)