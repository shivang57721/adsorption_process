"""
Simulation of adsorption equations with ModelingToolkit.jl
dU/dt + T/P * ∂/∂Z (P/T * v) = 0
dV/dt - Ω₁ * ∂²T/∂Z² + Ω₂ ∂/∂Z (v * P) = 0
U = log(P) - log(T)
V = T + Ω₂P
∂P/∂Z + F(v) = 0 (Ergun's equation)
"""
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

function simple_upwind!(f_zf, f, f_left, f_right)
    N = length(f)
    f_zf[1] = f_left
    f_zf[N+1] = f_right
    @inbounds for j in 1:(N-1)
        f_zf[j+1] = f[j]
    end
end

@constants begin
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

    Dₗ       = 0.7Dₘ + 0.5v₀ * 2rₚ     # Axial dispersion
    Pe   = v₀ * L / Dₗ                 # Peclet number
    Peₕ  = (ε * v₀ * L) / (K_gas / (ρ_gas[1] * Cₚ_gas))      # Heat Peclet number

    Ω₂ = (Cₚ_gas / R) * (P₀ / T₀) / ((1 - ε) / ε * (ρₛ * Cₚₛ))
end

@parameters begin
    N = 10,
    L = 1
end

@variables begin
    P̅[1:N](t) = 1.0
    T̅[1:N](t) = 1.0
    U[1:N](t) = 0.0
    V[1:N](t) = 1.0 + Ω₂

    # Velocity
    v[1](t) = 1.0
    v[2:N+1](t) = 0.0
end

@equations begin
    P̅_zf = zeros(N+1)
    P̅_left = P̅[1] + (150μ[1] * (L * ΔZ/2)/(4rₚ²) * (1 - ε)^2 / ε^2 * 1.0 + 
                  + 1.75 * (L * ΔZ/2) * ρ_gas[1] / (2rₚ) * (1 - ε) / ε * 1.0 * abs(1.0)) / P₀
    P̅_right = 1
    simple_upwind!(P̅_zf, P̅, P̅_left, P̅_right)

    T̅_zf = zeros(N+1)
    T̅_left = (T̅[1] + 1.0 * Peₕ * ΔZ/2) / (1 + 1.0 * Peₕ * ΔZ/2)
    T̅_right = T̅[N]
    simple_upwind!(T̅_zf, T̅, T̅_left, T̅_right)

    # dU/dt = - T/P * ∂/∂Z (P/T * v)
    U_flux = P̅_zf ./ T̅_zf .* v
    for j in 1:N
        D(U[j]) ~ - T[j] / P[j] * (U_flux[j+1] - U_flux[j])
    end
end