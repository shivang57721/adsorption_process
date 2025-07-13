# Computes weno of u at the faces; result stored in place in u_zf
# Boundary Conditions are provided
function WENO!(f_zf, f, f_left, f_right; clamp_result=false)
    N = length(f)
    f_zf[1]   = f_left
    f_zf[N+1] = f_right

    δ = 1e-10
    @inbounds for j in 1:N-1
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
        if clamp_result
            f_zf[j+1] = clamp(f_zf[j+1], 0.0, 1.0)
        end
    end
end

# Computes velocity from pressure drop using Ergun's equation
function compute_velocity(P̅; y_CO₂, y_H₂O, y_N₂, params, t, P̅_right=1)
    # Computes velocity from pressure drop using Ergun's equation
    N = params.N
    ΔZ = 1 / N

    ε = params.ε
    L = params.L
    rₚ = params.rₚ
    P₀ = params.Pₕ
    v₀ = params.v_feed
    T_feed = params.T_feed
    R = params.R
    μ_N₂   = params.μ_N₂
    μ_CO₂  = params.μ_CO₂
    μ_H₂O  = params.μ_H₂O

    # Velocity ramping to reduce initial stiffness
    ramp_time = 10.0  # Time to reach full inlet velocity
    v_feed_ramp = 0.5 * (1 + tanh(5 * (t - ramp_time/2) / ramp_time))

    # Calculate viscosity and gas density in each cell
    μ = params.buffers["μ"]
    ρ_gas = params.buffers["ρ_gas"]
    @. μ = y_CO₂ * μ_CO₂ + y_H₂O * μ_H₂O + y_N₂ * μ_N₂
    @. ρ_gas = P₀ / (R * T_feed) * (y_CO₂ * 44.009 + y_H₂O * 18.01528 + y_N₂ * 28.0134) * 1e-3

    b = 150μ * (1 - ε) ./ (2rₚ * 1.75ρ_gas)
    c = ε^3 * 2rₚ ./ (1.75ρ_gas * (1 - ε))

    v̅_zf = params.buffers["v̅_zf"]

    # Calculate superficial velocity using Ergun equation
    v̅_zf[1] = v_feed_ramp * ε
    @inbounds for j in 1:N
        if j == N
            ΔPΔZ = P₀ * (P̅_right - P̅[j]) / (ΔZ/2)
        else 
            ΔPΔZ = P₀ * (P̅[j+1] - P̅[j]) / ΔZ
        end

        v̅_zf[j+1] = - 1/2 * sign(ΔPΔZ) * (-b[j] + sqrt(b[j]^2 + sign(ΔPΔZ) * 4c[j] * ΔPΔZ)) / (v₀ * L)
    end

    # Return interstitial velocity
    return v̅_zf / ε
end