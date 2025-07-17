"""
Simulation of adsorption equations with ModelingToolkit.jl
dU/dt + T/P * ∂/∂Z (P/T * v) = 0
dV/dt - Ω₁ * ∂²T/∂Z² + Ω₂ ∂/∂Z (v * P) = 0
U = log(P) - log(T)
V = T + Ω₂P
∂P/∂Z + F(v) = 0 (Ergun's equation)
"""

using ModelingToolkit, MethodOfLines, OrdinaryDiffEq, DomainSets

@constants begin
    L = 1,
    r_in = 0.1445,
    r_out = 0.1620,
    ε = 0.37,
    rₚ = 1e-3,
    rₚ² = rₚ^2
    τ = 3.0,
    ρₛ = 1130,
    ρ_wall = 7800,
    Cₚ_gas = 30.7,
    Cₚ_ads = 30.7,
    Cₚₛ = 1070,
    Cₚ_wall = 502,
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

    P₀ = Pₕ,
    T₀ = T_feed,
    v₀ = v_feed,

    ρ_gas = P₀ / (R * T_feed) * 1.0 * 28.0134 * 1e-3

    Dₗ       = 0.7Dₘ + 0.5v₀ * 2rₚ     # Axial dispersion
    Pe   = v₀ * L / Dₗ                 # Peclet number
    Peₕ  = (ε * v₀ * L) / (K_gas / (ρ_gas * Cₚ_gas))      # Heat Peclet number

    Ω₁ = (K_gas / (v₀ * ε * L)) / ((1 - ε) / ε * (ρₛ * Cₚₛ))
    Ω₂ = (Cₚ_gas / R) * (P₀ / T₀) / ((1 - ε) / ε * (ρₛ * Cₚₛ))
end

@parameters t Z
@variables U(..) V(..) T(..) P(..) v(..)
Dt = Differential(t)
Dz = Differential(z)
Dzz = Differential(z)^2

eq = [
    # dU/dt + T/P * ∂/∂Z (P/T * v) = 0
    Dt(U(t, z)) ~ - T(t, z) / P(t, z) * Dz(P(t, z) / T(t, z) * v(t, z)),
    # dV/dt - Ω₁ * ∂²T/∂Z² + Ω₂ * ∂/∂Z (v * P) = 0
    Dt(V(t, z)) ~ Ω₁ * Dzz(T) - Ω₂ * Dz(v(t, z) * P(t, z)),
    # U = log(P) - log(T)
    0 ~ U(t, z) - log(P(t,z)) + log(T(t,z)),
    # V = T + Ω₂P
    0 ~ V(t,z) - T(t,z) - Ω₂ * P(t,z),
    # Ergun
    0 ~ Dz(P(t,z)) + (150 * μ_N₂ / (4rₚ²) * (1 - ε)^2 / ε^2 * v(t, z) + 
                    1.75 * ρ_gas[1] / (2rₚ) * (1 - ε) / ε * v(t, z)^2 ) / P₀
]

domains = [
    z ∈ Interval(0.0, 1.0)
    t ∈ Interval(0.0, 10.0)
]

# IC and BC
bcs = []