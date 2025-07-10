using Revise
# using DiffEqDiffTools
include("../src/AdsorptionModel.jl"); Revise.track("src/AdsorptionModel.jl")
using .AdsorptionModel

include("../src/IsothermModel.jl"); Revise.track("src/IsothermModel.jl")
using .IsothermModel

using Plots
using DifferentialEquations
using Sundials
using Random
Random.seed!(1) 

isotherm_params = IsothermParams(
    T₀ = 296,
    nₛ₀ = 2.38,
    b₀ = 70.74e4,
    t₀ = 0.4148,
    ΔH₀ = −57.047,
    χ = 0,
    γ_iso = 0.0061,
    β = 28.907,
    c_G = 0.1489,
    K_ads = 0.5751,
    cₘ = 36.48,
    R = 8.314,
    α = −1.606
)

params = AdsorptionParams(
    N = 200,
    L = 1,
    r_in = 0.1445,
    r_out = 0.1620,
    ε = 0.37,
    εₚ = 0.35,
    rₚ = 1e-3,
    τ = 3.0,
    ρₛ = 1130,
    ρ_wall = 7800,
    Cₚ_gas = 30.7,
    Cₚ_ads = 30.7,
    Cₚₛ = 1070,
    Cₚ_wall =  502,
    μ_N₂ = 1.789e-5,
    μ_CO₂ = 1.5003e-5,
    μ_H₂O = 1.0e-5,
    Dₘ = 1.60e-5,
    γ = 1.4,
    K_wall = 16,
    K_gas = 0.09,
    h_in = 8.6,
    h_out = 2.5,
    R = 8.314,
    y_feed = Dict("CO₂" => 0.0004, "N₂" => 0.8996, "H₂O" => 0.1), # DAC
    v_feed = 1,
    T_feed = 298.15,
    Tₐ = 298.15,
    Pₕ = 1.0e5,
    Pᵢ = 0.2e5,
    Pₗ = 0.1e5,
    qₛ₀ = 1,
    q_star_CO₂ = (T, P_CO₂, q_H₂O) -> IsothermModel.q_star_CO₂(T, P_CO₂, q_H₂O, isotherm_params), 
    q_star_H₂O = (T, P_H₂O) -> IsothermModel.q_star_H₂O(T, P_H₂O, isotherm_params),
    k_CO₂ = 0.0002,
    k_H₂O = 0.002,
    # k_CO₂ = 0,
    # k_H₂O = 0,
    ΔH_CO₂ = -57000,
    ΔH_H₂O = -49000
)

# Set initial conditions
N = params.N
y0_CO₂  = fill(0.0, N)
y0_N₂   = fill(1.0, N)
y0_H₂O  = fill(0.0, N)
x0_CO₂  = fill(0, N)
x0_N₂   = fill(0, N)
x0_H₂O  = fill(0, N)
P̅0      = fill(1.0, N)
T̅0      = fill(0.9, N)
T̅_wall0 = fill(1.0, N)

u0 = vcat(y0_CO₂, y0_N₂, y0_H₂O,
            x0_CO₂, x0_N₂, x0_H₂O,
            P̅0, T̅0, T̅_wall0)
tspan = (0.0, 1000.0)

prob = ODEProblem(adsorption_equations!, u0, tspan, params)

sol = solve(prob, TRBDF2(autodiff = AutoFiniteDiff()); saveat=0.01, verbose = true)

# Unpack sol
y_CO₂  = t -> sol(t)[1:N]
y_N₂   = t -> sol(t)[N+1:2N]
y_H₂O  = t -> sol(t)[2N+1:3N]
y      = t -> Dict("CO₂" => y_CO₂(t), "N₂" => y_N₂(t), "H₂O" => y_H₂O(t))
x_CO₂  = t -> sol(t)[3N+1:4N]
x_N₂   = t -> sol(t)[4N+1:5N]
x_H₂O  = t -> sol(t)[5N+1:6N]
P̅      = t -> sol(t)[6N+1:7N]
T̅      = t -> sol(t)[7N+1:8N]
T̅_wall = t -> sol(t)[8N+1:9N]
v̅_zf   = t -> compute_velocity(P̅(t), y(t), params)

plot(P̅(sol.t[end]), title="pressure")
plot(v̅_zf(1000), title="velocity")
plot(T̅(1000), title="temperature")
plot(T̅_wall(1000), title="temperature wall")


# plot(compute_velocity(P̅(sol.t[end]), params))