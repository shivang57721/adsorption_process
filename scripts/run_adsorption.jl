using Revise

include("../src/AdsorptionModel.jl"); Revise.track("src/AdsorptionModel.jl")
using .AdsorptionModel

include("../src/util.jl")

using Plots
using DifferentialEquations
using Sundials
using Random
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
P̅0      = fill(1.0, N)
T̅0      = fill(1.0, N)
T̅_wall0 = fill(1.0, N)

u0 = vcat(y0_CO₂, y0_N₂, y0_H₂O,
            x0_CO₂, x0_H₂O,  # x0_N₂ removed
            P̅0, T̅0, T̅_wall0)
tspan = (0.0, 200.0)

prob = ODEProblem(adsorption_equations!, u0, tspan, params)

@time sol = solve(prob, QNDF(autodiff = AutoFiniteDiff());
            abstol = 1e-6,
            reltol = 1e-5,
            saveat = 0.01,
            verbose = true,
            dtmin = 1e-15,
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
v̅_zf   = t -> compute_velocity(P̅(t); y_CO₂ = y_CO₂(t), y_N₂ = y_N₂(t), y_H₂O = y_H₂O(t), params)

plot(P̅(sol.t[end]), title="pressure")
plot(v̅_zf(200), title="velocity")
plot(T̅(200), title="temperature")
plot(y_CO₂(200), title="y_CO₂")
plot(y_H₂O(200), title="y_CO₂")
plot(y_N₂(200), title="y_N₂")

plot(T̅_wall(200), title="temperature wall")
