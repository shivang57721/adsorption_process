module PressureVelocity

include("AdsorptionModel.jl")
using .AdsorptionModel

export pressure_velocity!, PressureParams, compute_velocity

Base.@kwdef struct PressureParams
    N::Int64
    L::Int64
    ε::Float64
    rₚ::Float64
    v_feed::Int64
    Pₕ::Float64
    T_feed::Float64
    μ_N₂::Float64
    μ_CO₂::Float64
    μ_H₂O::Float64
    R::Float64
end

function pressure_velocity!(dP̅, P̅, params, t)
    # Discretization parameters
    N      = params.N               # Number of finite volume elements
    ΔZ     = 1 / N                  # Space discretization  

    ε      = params.ε
    rₚ     = params.rₚ
    rₚ²    = rₚ^2
    P₀     = params.Pₕ
    T_feed = params.T_feed

    μ_N₂   = params.μ_N₂
    μ_CO₂  = params.μ_CO₂
    μ_H₂O  = params.μ_H₂O

    v₀     = params.v_feed
    L      = params.L
    R      = params.R

    y_CO₂ = fill(0.0, N)
    y_H₂O = fill(0.0, N)
    y_N₂  = fill(1.0, N)
    y = Dict("CO₂" => y_CO₂, "H₂O" => y_H₂O, "N₂" => y_N₂)

    # Compute P̅ and T̅ and y at cell faces using WENO
    μ = y_CO₂ * μ_CO₂ .+ y_H₂O * μ_H₂O .+ y_N₂ * μ_N₂
    ρ_gas = P₀ / (R * T_feed) * (y_CO₂ * 44.009 .+ y_H₂O * 18.01528 .+ y_N₂ * 28.0134) * 1e-3
    P̅_zf = similar(P̅, N + 1)
    P̅_left = P̅[1] + (150μ[1] * (ΔZ/2)/(4rₚ²) * (1 - ε)^2 / ε^2 * v₀ + 1.75 * (ΔZ/2) * ρ_gas[1] / (2rₚ) * (1 - ε) / ε * v₀ * abs(v₀)) / P₀
    P̅_right = 1
    WENO!(P̅_zf, P̅, P̅_left, P̅_right)

    # Compute local velcoity at cell faces
    v̅_zf = compute_velocity(P̅, y, params)

    P̅_flux = P̅_zf .* v̅_zf

    for j in 1:N
        dP̅[j] = - 1 / ΔZ * (P̅_flux[j+1] - P̅_flux[j])
    end
end


end