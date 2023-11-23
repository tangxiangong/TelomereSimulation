include("./src/main.jl")
using Plots
pgfplotsx()

#= 模型参数 =#
begin
    l₀ = 0
    α = 0.7
    β = 1.3
    γ = 1
    λ = 1
    μ = 0
    σ = 1
    T = 100
    domain = (-10, 10)
end
#= 1. 随机游走 =#

rw = TelomereRW((l₀, α, λ, μ, σ)) 
# 相对应的路径
rwₜ = rw(T)
# 模拟路径
simulate(rwₜ, 0.1)
# 阶矩
𝔼(rwₜ,τ=0.1), 𝔼(rwₜ^2)
# 随机变量：首次通过时间
ot = OccupationTime(T, domain, rw)
simulate(ot)
fpt = FPT(domain, rw)
# 模拟
simulate(fpt)
# 阶矩
𝔼(fpt), 𝔼(fpt^2)

#= 2. Lagevin 方程 =#

langevin = TelomereLangevin((l₀, α, λ, μ, σ)) 
# 相对应的路径
langevinₜ = langevin(T)
# 模拟路径
simulate(langevinₜ, 0.1)
# 阶矩
𝔼(langevinₜ), 𝔼(langevinₜ^2)
# 随机变量：首次通过时间
fptl = FPT(domain, langevin)
# 模拟
simulate(fptl)
# 阶矩
𝔼(fptl), 𝔼(fptl^2)

#= 3. 时间回火随机游走 =#
ttrw = TelomereTTRW((l₀, α, γ, λ, μ, σ)) 
# 相对应的路径
ttrwₜ = ttrw(T)
# 模拟路径
simulate(ttrwₜ)
# 阶矩
𝔼(ttrwₜ), 𝔼(ttrwₜ^2)
# 随机变量：首次通过时间
fpttt = FPT(domain, rw)
# 模拟
simulate(fpttt)
# 阶矩
𝔼(fpttt), 𝔼(fpttt^2)

#= 4. 空间回火随机游走 =#
strw = TelomereSTRW((l₀, α, β, γ, λ, μ, σ)) 
# 相对应的路径
strwₜ = strw(T)
# 模拟路径
simulate(strwₜ)
# 阶矩
𝔼(strwₜ), 𝔼(strwₜ^2)
# 随机变量：首次通过时间
fptst = FPT(domain, strw)
# 模拟
simulate(fptst)
# 阶矩
𝔼(fptst; N=1_000)
𝔼(fptst^2)
