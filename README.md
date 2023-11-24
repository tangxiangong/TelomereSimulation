# Telomere Simulation in Julia
## 例子
给定模型参数.
```julia
include("./src/main.jl")
l₀ = 0
α = 0.7
λ = 1
μ = 0
σ = 1
T = 10
Δ = 1
domain = (-10, 10)
```
实例化相应的模型与路径， 通过函数 `simulate` 模拟路径. 
```julia
x = TelomereLangevin((l₀, α, λ, μ, σ))
xₜ = x(T)
simulate(xₜ, 0.1)
```
使用函数 `𝔼` 计算相应的阶矩.
```julia
𝔼(xₜ), 𝔼(xₜ^2)
```
模拟相应模型的首次通过时间和占据时间并计算其阶矩.
```julia
fpt = FPT(domain, x)
occuptime = OccupationTime(T, domain, x)
simulate(fpt) 
simulate(occuptime)
𝔼(fpt)
𝔼(occuptime)
```
计算时间平均均方位移.
```julia
𝔼(δ̄²(x, T, Δ))
```
对于一族时刻, 计算其阶矩.
```julia
t = collect(10:10:100)
@. 𝔼(x(t))
@. 𝔼((x(t))^2)
```

## 第三方依赖
```
PoissonRandom
FastGaussQuadrature
```

## API
### 自定义类型
- 抽象类型 `StochasticProcess`
- `Trajectory(::StochasticProcess, ::Real)`  路径类型, 参数为具体的随机过程实例及路径长度; 还可由 `(::StochasticProcess)(::Real)` 实例化
- `TelomereRW(::NTuple{5, Real}) <: StochasticProcess` 随机游走类型, 参数为 `(x₀, α, λ, μ, σ)`
- `TelomereTTRW(::NTuple{6, Real}) <: StochasticProcess` 时间回火随机游走类型, 参数为 `(x₀, α, γ, λ, μ, σ)`
- `TelomereSTRW(::NTuple{7, Real}) <: StochasticProcess` 时间回火随机游走类型, 参数为 `(x₀, α, β, γ, λ, μ, σ)`
- `TelomereLangevin(::NTuple{5, Real}) <: StochasticProcess` Langevin 方程类型, 参数为 `(x₀, α, λ, μ, σ)`
- 抽象类型 `Functional`
- `FPT(::NTuple{2, Real}, ::StochasticProcess) <: Functional` 首次通过时间
- `OccupationTime(::NTuple{2, Real}, ::StochasticProcess) <: Functional` 占据时间
- `TimeAverage(::StochasticProcess, ::Real, ::Real)` 时间平均, alias `δ̄²`

### 函数/方法
模拟函数, 默认参数为模拟步长(对于 Langevin 方程型随机过程)
- `simulate(::Trajectory, ::Float64=1e-2)`
- `simulate(::FPT, ::Float64=1e-2)`
- `simulate(::OccupationTime, ::Float64=1e-2)`

期望函数， 默认关键字为蒙特卡罗模拟次数 `N` 以及模拟步长 `τ`
- `𝔼(::Trajectory; N::Int=100_000, τ::Float64=0.01)`
- `𝔼(::Functional; N::Int=100_000, τ::Float64=0.01)`
- `𝔼(::TimeAverage; N::Int=100_000, τ::Float64=0.01)`

