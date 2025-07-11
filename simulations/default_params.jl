# simulations/sim1_params.jl
include("../src/params.jl")
include("../src/IsothermModel.jl")

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

N = 100

params = AdsorptionParams(
    N = N,
    L = 0.1,
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
    y_feed_CO₂ = 0.0004,
    y_feed_N₂ = 0.8996,
    y_feed_H₂O = 0.1,
    v_feed = 1,
    T_feed = 298.15,
    Tₐ = 298.15,
    Pₕ = 1.0e5,
    Pᵢ = 0.2e5,
    Pₗ = 0.1e5,
    qₛ₀ = 1,
    q_star_CO₂ = (T, P_CO₂, q_H₂O) -> q_star_CO₂(T, P_CO₂, q_H₂O, isotherm_params), 
    q_star_H₂O = (T, P_H₂O) -> q_star_H₂O(T, P_H₂O, isotherm_params),
    # k_CO₂ = 0.0002,
    # k_H₂O = 0.002,
    k_CO₂ = 0,
    k_H₂O = 0,
    ΔH_CO₂ = -57000,
    ΔH_H₂O = -49000
)

# Buffers
params.buffers["μ"] = zeros(N)
params.buffers["ρ_gas"] = zeros(N)
params.buffers["Σxᵢ"] = zeros(N)
params.buffers["Ω₁"] = zeros(N)
params.buffers["Ω₂"] = zeros(N)
params.buffers["Ω₃"] = zeros(N)
params.buffers["Ω₄"] = zeros(N)
params.buffers["σ_CO₂"] = zeros(N)
params.buffers["σ_H₂O"] = zeros(N)
params.buffers["σ_N₂"] = zeros(N)

params.buffers["P̅_zf"] = zeros(N+1)
params.buffers["T̅_zf"] = zeros(N+1)
params.buffers["T̅_flux"] = zeros(N+1)
params.buffers["P̅_flux"] = zeros(N+1)
params.buffers["y_flux"] = zeros(N+1)