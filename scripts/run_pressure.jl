using Revise

include("../src/PressureVelocity.jl"); Revise.track("src/PressureVelocity.jl")
using.PressureVelocity

using Plots
using DifferentialEquations
using Sundials
using Random
Random.seed!(1) 

# Computes velocity from pressure drop using Ergun's equation
function compute_velocity(P̅, y, params; P̅_right=1)
    # Computes velocity from pressure drop using Ergun's equation
    N = params.N
    ΔZ = 1 / N

    ε = params.ε
    rₚ = params.rₚ
    P₀ = params.Pₕ
    v₀ = params.v_feed
    T_feed = params.T_feed
    R = params.R
    μ_N₂   = params.μ_N₂
    μ_CO₂  = params.μ_CO₂
    μ_H₂O  = params.μ_H₂O

    if y === nothing
        y_CO₂ = fill(0.0, N)
        y_H₂O = fill(0.0, N)
        y_N₂ = fill(1.0, N)
    else
        y_CO₂ = y["CO₂"]
        y_H₂O = y["H₂O"]
        y_N₂  = y["N₂"]
    end

    μ = y_CO₂ * μ_CO₂ .+ y_H₂O * μ_H₂O .+ y_N₂ * μ_N₂ # viscosity in each cell
    ρ_gas = P₀ / (R * T_feed) * (y_CO₂ * 44.009 .+ y_H₂O * 18.01528 .+ y_N₂ * 28.0134) * 1e-3 # gas density in each cell

    b = 150μ * (1 - ε) ./ (2rₚ * 1.75ρ_gas)
    c = ε^3 * 2rₚ ./ (1.75ρ_gas * (1 - ε))

    v̅_zf = similar(P̅, N + 1)

    # Calculate superficial velocity using Ergun equation
    v̅_zf[1] = 1 * ε
    for j in 1:N
        if j == N
            ΔPΔZ = P₀ * (P̅_right - P̅[j]) / (ΔZ/2)  
        else 
            ΔPΔZ = P₀ * (P̅[j+1] - P̅[j]) / ΔZ
        end

        v̅_zf[j+1] = - 1/2 * sign(ΔPΔZ) * (-b[j] + sqrt(b[j]^2 + sign(ΔPΔZ)* 4c[j] * ΔPΔZ)) / v₀
    end

    # Return interstitial velocity
    return v̅_zf / ε
end

function WENO!(f_zf, f, f_left, f_right)
    N = length(f)
    f_zf[1]   = f_left
    f_zf[N+1] = f_right

    δ = 1e-10
    for j in 1:N-1
        # Calculate f at face j + 1/2
        if j == 1
            # Use f₁ - f₀ = 2(f₁ - f₀₅)
            f₀ = 2 * f_left - f[1]
            α₀ⱼ = (2/3) / (f[j+1] - f[j] + δ)^4
            α₁ⱼ = (1/3) / (f[j] - f₀ + δ)^4
            value = α₀ⱼ / (α₀ⱼ + α₁ⱼ) * (1/2 * (f[j] + f[j+1])) +
                        α₁ⱼ / (α₀ⱼ + α₁ⱼ) * (3/2 * f[j] - 1/2 * f₀)
        else
            # interior points
            α₀ⱼ = (2/3) / (f[j+1] - f[j] + δ)^4 
            α₁ⱼ = (1/3) / (f[j] - f[j-1] + δ)^4
            value = α₀ⱼ / (α₀ⱼ + α₁ⱼ) * (1/2 * (f[j] + f[j+1])) +
                        α₁ⱼ / (α₀ⱼ + α₁ⱼ) * (3/2 * f[j] - 1/2 * f[j-1])
        end
        f_zf[j+1] = isnan(value) ? 1/2 * (f[j] + f[j+1]) : value
    end
end

params = PressureParams(
    N = 200,
    L = 1,
    ε = 0.37,
    rₚ = 1e-3,
    v_feed = 1,
    Pₕ = 1.0e5,
    T_feed = 298.15,
    μ_N₂ = 1.789e-5,
    μ_CO₂ = 1.5003e-5,
    μ_H₂O = 1.0e-5,
    R = 8.314
)

N  = params.N
P = fill(1.0, N)
# Δ = (1.1 - 1) / N
# P = 1.1 .- Δ * Array(1:N)

tspan = (0.0, 1000.0)
prob = ODEProblem(pressure_velocity!, P, tspan, params)

sol = solve(prob,  Rosenbrock23(); saveat=0.001, verbose = true)

plot(sol(1), title="pressure")

plot(compute_velocity(sol(300), nothing, params), title="velocity")

###############################################


P̅ = sol(1000)
v̅_zf = compute_velocity(sol(1000), nothing, params)

rₚ = params.rₚ
R = params.R
T = params.T_feed

ε      = params.ε
rₚ     = params.rₚ
rₚ²    = rₚ^2
P₀     = params.Pₕ
T_feed = params.T_feed
v₀     = params.v_feed
L      = params.L
μ_CO₂ = params.μ_CO₂
μ_H₂O = params.μ_H₂O
μ_N₂ = params.μ_N₂
ΔZ = 1/N

y_CO₂ = fill(0.0, N)
y_H₂O = fill(0.0, N)
y_N₂  = fill(1.0, N)

μ = y_CO₂ * μ_CO₂ .+ y_H₂O * μ_H₂O .+ y_N₂ * μ_N₂ # viscosity in each cell
ρ_gas = P₀ / (R * T_feed) * (y_CO₂ * 44.009 .+ y_H₂O * 18.01528 .+ y_N₂ * 28.0134) * 1e-3
P̅_zf = similar(P̅, N + 1)
P̅_left = P̅[1] - (150μ[1] * (ΔZ/2)/(4rₚ²) * (1 - ε)^2 / ε^3 * v₀ + 1.75 * (ΔZ/2) * ρ_gas[1] / (2rₚ) * (1 - ε) / ε * v₀ * abs(v₀)) / P₀
P̅_right = 1
WENO!(P̅_zf, P̅, P̅_left, P̅_right)

S = π * rₚ^2
n = v̅_zf .* S .* P̅_zf ./ (R * T)