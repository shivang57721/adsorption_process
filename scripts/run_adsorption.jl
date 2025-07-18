using Revise

include("../src/AdsorptionModel.jl"); Revise.track("src/AdsorptionModel.jl")

include("../src/util.jl"); Revise.track("src/util.jl")

using Plots
using DifferentialEquations
using BenchmarkTools
using Sundials
using Random
using ODEInterfaceDiffEq
using Profile
Random.seed!(1) 

# Load params file
include("../simulations/default_params.jl")

# Set initial conditions
N = params.N
y0_CO₂  = fill(0.0, N)
y0_N₂   = fill(1.0, N)
y0_H₂O  = fill(0.0, N)
x0_CO₂  = fill(0, N)
# x0_N₂ removed since it's always zero
x0_H₂O  = fill(0, N)

# v0(z) = exp(-100 * z^2)
P̅0      = fill(1.0, N)
T̅0      = fill(0.95, N)
T̅_wall0 = fill(1.0, N)

u0 = vcat(y0_CO₂, y0_N₂, y0_H₂O,
            x0_CO₂, x0_H₂O,  # x0_N₂ removed
            P̅0, T̅0, T̅_wall0)
tspan = (0.0, 200.0)

prob = ODEProblem(adsorption_equations!, u0, tspan, params)

@time sol = solve(prob, Rodas5P(autodiff = AutoFiniteDiff());
            reltol = 1e-6,
            abstol = 1e-6,
            saveat = 0.01,
            verbose = true,
            dtmin = 1e-12,
            isoutofdomain = (u,p,t) -> minimum(u) < -1e-12)

# Unpack sol
y_CO₂  = t -> sol(t)[1:N]
y_N₂   = t -> sol(t)[N+1:2N]
y_H₂O  = t -> sol(t)[2N+1:3N]
x_CO₂  = t -> sol(t)[3N+1:4N]
x_H₂O  = t -> sol(t)[4N+1:5N]
P̅      = t -> sol(t)[5N+1:6N]
T̅      = t -> sol(t)[6N+1:7N]
T̅_wall = t -> sol(t)[7N+1:8N]
# x_N₂ is implicitly zero everywhere
x_N₂   = t -> zeros(N)
v̅_zf   = t -> begin
    buffer = zeros(N+1)
    compute_velocity!(buffer, P̅(t); y_CO₂ = y_CO₂(t), y_N₂ = y_N₂(t), y_H₂O = y_H₂O(t), params, t)
    buffer
end

plot(P̅(10), title="pressure")

plot(v̅_zf(10), title="velocity")

plot(T̅(100), title="temperature")

plot(y_CO₂(1), title="y_CO₂")

plot(y_H₂O(1), title="y_H₂O")
plot(y_N₂(1), title="y_N₂")

plot(T̅_wall(200), title="temperature wall")
