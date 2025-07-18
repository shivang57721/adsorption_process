"""
Simulation of adsorption equations with ModelingToolkit.jl
dU/dt + T/P * ∂/∂Z (P/T * v) = 0
dV/dt - Ω₁ * ∂²T/∂Z² + Ω₂ ∂/∂Z (v * P) = 0
U = log(P) - log(T)
V = T + Ω₂P
∂P/∂Z + F(v) = 0 (Ergun's equation)
"""

using ModelingToolkit, MethodOfLines, OrdinaryDiffEq, DomainSets

# Constants 
L = 1
r_in = 0.1445
r_out = 0.1620
ε = 0.37
rₚ = 1e-3
rₚ² = rₚ^2
τ = 3.0
ρₛ = 1130
ρ_wall = 7800
Cₚ_gas = 30.7
Cₚ_ads = 30.7
Cₚₛ = 1070
Cₚ_wall = 502
μ_N₂ = 1.789e-5
μ_CO₂ = 1.5003e-5
μ_H₂O = 1.0e-5
Dₘ = 1.60e-5
γ = 1.4
K_wall = 16
K_gas = 0.09
h_in = 8.6
h_out = 2.5
R = 8.314
y_feed_CO₂ = 0.0004
y_feed_N₂ = 0.8996
y_feed_H₂O = 0.1
v_feed = 1
T_feed = 298.15
Tₐ = 298.15
Pₕ = 1.0e5

P₀ = Pₕ
T₀ = T_feed
v₀ = v_feed

ρ_gas = P₀ / (R * T_feed) * 1.0 * 28.0134 * 1e-3

Dₗ       = 0.7Dₘ + 0.5v₀ * 2rₚ     # Axial dispersion
Pe   = v₀ * L / Dₗ                 # Peclet number
Peₕ  = (ε * v₀ * L) / (K_gas / (ρ_gas * Cₚ_gas))      # Heat Peclet number

Ω₁ = (K_gas / (v₀ * ε * L)) / ((1 - ε) / ε * (ρₛ * Cₚₛ))
Ω₂ = (Cₚ_gas / R) * (P₀ / T₀) / ((1 - ε) / ε * (ρₛ * Cₚₛ))

@parameters t z
@variables U(..) V(..) T(..) P(..) v(..)
Dt = Differential(t)
Dz = Differential(z)
Dzz = Differential(z)^2

a = - 150 * μ_N₂ / (4rₚ²) * (1 - ε)^2 / ε^2 / P₀
b = - 1.75 * ρ_gas / (2rₚ) * (1 - ε) / ε / P₀
eq = [
    # dU/dt + T/P * ∂/∂Z (P/T * v) = 0
    Dt(U(t, z)) ~ - T(t, z) / P(t, z) * Dz((P(t, z) / T(t, z) * v(t, z))),
    # dV/dt - Ω₁ * ∂²T/∂Z² + Ω₂ * ∂/∂Z (v * P) = 0
    Dt(V(t, z)) ~ Ω₁ * Dzz(T(t, z)) - Ω₂ * Dz(v(t, z) * P(t, z)),
    # U = log(P) - log(T)
    0 ~ U(t, z) - log(P(t,z)) + log(T(t,z)),
    # V = T + Ω₂P
    0 ~ V(t,z) - T(t,z) - Ω₂ * P(t,z),
    # Ergun
    Dz(P(t,z)) ~ a * v(t, z) + b * v(t, z)^2
]

domains = [
    z ∈ Interval(0.0, 1.0)
    t ∈ Interval(0.0, 10.0)
]

v0(z) = exp(-100 * z^2)
DP0 = a * 1.0 + b * 1.0^2
# IC and BC
bcs = [P(0, z) ~ 1.0,
       T(0, z) ~ 1.0,
       v(0, z) ~ v0(z),
       U(0, z) ~ 0.0,
       V(0, z) ~ 1.0 + Ω₂ * 1.0,
       
       P(t, 1) ~ 1.0,
       Dz(P(t, 0)) ~ a * v(t, 0) + b * v(t, 0)^2,
       Dz(T(t, 1)) ~ 0.0,
       Dz(T(t, 0)) ~ - Peₕ * (1.0 - T(t, 0)),
       v(t, 0) ~ 1.0]

@named pdesys = PDESystem(eq, bcs, domains, [t, z], [U(t, z), V(t, z), T(t, z), P(t, z), v(t, z)])

N = 20

order = 2 # This may be increased to improve accuracy of some schemes

# Integers for x and y are interpreted as number of points. Use a Float to directtly specify stepsizes dx and dy.
discretization = MOLFiniteDifference([z=>N], t, approx_order=order)

println("Discretization:")
@time prob = discretize(pdesys,discretization)