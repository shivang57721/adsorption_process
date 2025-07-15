"""
DAE formulation of transport equations for adsorption modelling
This formulation includes algebraic constraints for y_N₂ normalization and velocity computation.
"""
module AdsorptionModelDAE
export adsorption_dae!, setup_dae_initial_conditions

include("util.jl")

function adsorption_dae!(res, du, u, params, t)
    #────────────────────────────────────────────
    # Unpack parameters
    #────────────────────────────────────────────
    
    # Discretization parameters
    N      = params.N               # Number of finite volume elements
    ΔZ     = 1 / N                  # Space discretization   
    
    # Column parameters
    L      = params.L               # Column length
    r_in   = params.r_in            # Inner column radius
    r_out  = params.r_out           # Outer column radius
    ε      = params.ε               # Bed porosity
    εₚ     = params.εₚ
    rₚ     = params.rₚ              # Particle radius
    τ      = params.τ               # Tortuosity
    rₚ²    = rₚ^2

    # Properties and constants
    ρₛ          = params.ρₛ              # Adsorbent density
    ρ_wall      = params.ρ_wall         # Column wall density
    Cₚ_gas      = params.Cₚ_gas          # Specific heat capacity of gas phase
    Cₚ_ads      = params.Cₚ_ads          # Specific heat capacity of adsorbed phase
    Cₚₛ          = params.Cₚₛ             # Specific heat capacity of adsorbent
    Cₚ_wall     = params.Cₚ_wall         # Specific heat capacity of column wall
    Dₘ          = params.Dₘ             # Molecular diffusivity
    γ           = params.γ              # Adiabatic constant
    K_wall      = params.K_wall         # Thermal conductivity of column wall
    K_gas       = params.K_gas          # Effective gas thermal conductivity
    h_in        = params.h_in           # Inside heat transfer coefficient
    h_out       = params.h_out          # Outside heat transfer coefficient
    R           = params.R              # Universal gas constant
    y_feed_CO₂  = params.y_feed_CO₂
    y_feed_N₂   = params.y_feed_N₂
    y_feed_H₂O  = params.y_feed_H₂O
    μ_N₂        = params.μ_N₂           # Fluid viscosity
    μ_CO₂       = params.μ_CO₂
    μ_H₂O       = params.μ_H₂O

    # Operating conditions for case studies
    v_feed  = params.v_feed             # Interstitial feed velocity
    T_feed  = params.T_feed             # Feed temperature
    Tₐ      = params.Tₐ                 # Ambient temperature
    Pₕ      = params.Pₕ                  # High pressure 
    Pᵢ      = params.Pᵢ                 # Intermediate pressure
    Pₗ       = params.Pₗ                 # Low pressure

    # Constants for dimensionless parameteres
    P₀     = Pₕ
    T₀     = T_feed
    v₀     = v_feed
    qₛ₀    = params.qₛ₀

    T̅_feed  = T_feed / T₀
    T̅ₐ      = Tₐ / T₀ 
    Dₚ      = Dₘ / τ                  # Macropore diffusivity
    Dₗ       = 0.7Dₘ + 0.5v₀ * 2rₚ     # Axial dispersion

    # Heat of adsorption
    ΔH_CO₂      = params.ΔH_CO₂
    ΔH_H₂O      = params.ΔH_H₂O

    # Isotherms
    q_star_CO₂ = params.q_star_CO₂ # function of temperature, partial pressure and q_H₂O
    q_star_H₂O = params.q_star_H₂O # function of temperature and partial pressure
    k_CO₂      = params.k_CO₂
    k_H₂O      = params.k_H₂O

    # Buffers
    μ = params.buffers["μ"]
    ρ_gas = params.buffers["ρ_gas"]
    Σxᵢ = params.buffers["Σxᵢ"]

    Ω₁ = params.buffers["Ω₁"]
    Ω₂ = params.buffers["Ω₂"]
    Ω₃ = params.buffers["Ω₃"]
    Ω₄ = params.buffers["Ω₄"]
    σ_CO₂ = params.buffers["σ_CO₂"]
    σ_H₂O = params.buffers["σ_H₂O"]

    #────────────────────────────────────────────
    # State vector layout for DAE:
    # Differential variables (du != 0): [y_CO₂, y_H₂O, x_CO₂, x_H₂O, P̅, T̅, T̅_wall]
    # Algebraic variables (du = 0):     [y_N₂, v̅]
    # Note: x_N₂ eliminated since dx_N₂ = 0 everywhere
    #────────────────────────────────────────────
    
    # Differential variable indices
    y_CO₂_idx = 1:N
    y_H₂O_idx = N+1:2N
    x_CO₂_idx = 2N+1:3N
    x_H₂O_idx = 3N+1:4N
    P̅_idx = 4N+1:5N
    T̅_idx = 5N+1:6N
    T̅_wall_idx = 6N+1:7N
    
    # Algebraic variable indices
    y_N₂_idx = 7N+1:8N         # y_N₂ (algebraic)
    v̅_idx = 8N+1:8N+(N+1)      # v̅ at faces (algebraic)

    # Extract state variables
    y_CO₂  = @view u[y_CO₂_idx]
    y_H₂O  = @view u[y_H₂O_idx]
    x_CO₂  = @view u[x_CO₂_idx]
    x_H₂O  = @view u[x_H₂O_idx]
    P̅      = @view u[P̅_idx]
    T̅      = @view u[T̅_idx]
    T̅_wall = @view u[T̅_wall_idx]
    y_N₂   = @view u[y_N₂_idx]      # Algebraic
    v̅_zf   = @view u[v̅_idx]        # Algebraic
    
    # x_N₂ is implicitly zero everywhere (not stored in state vector)

    # Extract residuals
    res_y_CO₂  = @view res[y_CO₂_idx]
    res_y_H₂O  = @view res[y_H₂O_idx]
    res_x_CO₂  = @view res[x_CO₂_idx]
    res_x_H₂O  = @view res[x_H₂O_idx]
    res_P̅      = @view res[P̅_idx]
    res_T̅      = @view res[T̅_idx]
    res_T̅_wall = @view res[T̅_wall_idx]
    res_y_N₂   = @view res[y_N₂_idx]    # Algebraic constraint
    res_v̅      = @view res[v̅_idx]      # Algebraic constraint

    # Extract derivatives (only for differential variables)
    dy_CO₂  = @view du[y_CO₂_idx]
    dy_H₂O  = @view du[y_H₂O_idx]
    dx_CO₂  = @view du[x_CO₂_idx]
    dx_H₂O  = @view du[x_H₂O_idx]
    dP̅      = @view du[P̅_idx]
    dT̅      = @view du[T̅_idx]
    dT̅_wall = @view du[T̅_wall_idx]

    #────────────────────────────────────────────
    # Compute auxiliary variables
    #────────────────────────────────────────────
    
    # Dimensionless groups
    Pe   = v₀ * L / Dₗ                                       # Peclet number
    Peₕ  = (ε * v₀ * L) / (K_gas / (ρ_gas[1] * Cₚ_gas))      # Heat Peclet number
    ψ    = (R * T₀ * qₛ₀) / Pₕ * (1 - ε) / ε * ρₛ

    @. Σxᵢ  = x_CO₂ + x_H₂O
    @. Ω₁   = (K_gas / (v₀ * ε * L)) / ((1 - ε) / ε * (ρₛ * Cₚₛ + ρ_gas * qₛ₀ * Cₚ_ads * Σxᵢ))
    @. Ω₂   = ((Cₚ_gas / R) * (P₀ / T₀)) / ((1 - ε) / ε * (ρₛ * Cₚₛ + ρ_gas * qₛ₀ * Cₚ_ads * Σxᵢ))
    @. Ω₃   = (ρ_gas * Cₚ_ads * qₛ₀) / (ρₛ * Cₚₛ + ρ_gas * qₛ₀ * Cₚ_ads * Σxᵢ)
    @. Ω₄   = (2h_in / r_in) * (L / v₀) / ((1 - ε) * (ρₛ * Cₚₛ + ρ_gas * qₛ₀ * Cₚ_ads * Σxᵢ))
    @. σ_CO₂ = (ρ_gas * qₛ₀ / T₀) * (- ΔH_CO₂) / (ρₛ * Cₚₛ + ρ_gas * qₛ₀ * Cₚ_ads * Σxᵢ)
    @. σ_H₂O = (ρ_gas * qₛ₀ / T₀) * (- ΔH_H₂O) / (ρₛ * Cₚₛ + ρ_gas * qₛ₀ * Cₚ_ads * Σxᵢ)

    Π₁   = K_wall / (ρ_wall * Cₚ_wall * v₀ * L)
    Π₂   = 2r_in * h_in / (r_out^2 - r_in^2) * L / (ρ_wall * Cₚ_wall * v₀)
    Π₃   = 2r_out * h_out / (r_out^2 - r_in^2) * L / (ρ_wall * Cₚ_wall * v₀)

    #────────────────────────────────────────────
    # Algebraic constraints
    #────────────────────────────────────────────
    
    # Constraint 1: y_CO₂ + y_H₂O + y_N₂ = 1
    @. res_y_N₂ = y_N₂ - (1.0 - y_CO₂ - y_H₂O)
    
    # μ and ρ_gas computation (needed for velocity)
    @. μ = y_CO₂ * μ_CO₂ + y_H₂O * μ_H₂O + y_N₂ * μ_N₂
    @. ρ_gas = P₀ / (R * T_feed) * (y_CO₂ * 44.009 + y_H₂O * 18.01528 + y_N₂ * 28.0134) * 1e-3

    res_v̅[1] = v̅_zf[1] - 1.0
    
    # Velocity constraints using Ergun equation
    @inbounds for j in 1:N
        if j == N
            ΔPΔZ = P₀ / L * (1.0 - P̅[j]) / (ΔZ/2)  # P̅_right = 1
        else 
            ΔPΔZ = P₀ / L * (P̅[j+1] - P̅[j]) / ΔZ
        end
        
        # Ergun equation: ΔP/ΔZ = 150μ/(4rₚ²) * (1-ε)²/ε² * v + 1.75ρ/(2rₚ) * (1-ε)/ε * v|v|
        res_v̅[j+1] = ΔPΔZ - (150 * μ[j]) / (4rₚ²) * (1 - ε)^2 / ε^2 * v̅_zf[j+1] - 
                          1.75 * ρ_gas[j] / (2rₚ) * (1 - ε) / ε * v̅_zf[j+1] * abs(v̅_zf[j+1])
    end

    # Compute P̅ and T̅ and y at cell faces using WENO
    P̅_zf = params.buffers["P̅_zf"]
    P̅_left = P̅[1] + (150μ[1] * (L * ΔZ/2)/(4rₚ²) * (1 - ε)^2 / ε^2 * v̅_zf[1] + 
                        1.75 * (L * ΔZ/2) * ρ_gas[1] / (2rₚ) * (1 - ε) / ε * v̅_zf[1] * abs(v̅_zf[1])) * (v₀ / P₀)
    # P̅_left = P̅[1] + (v̅_zf[1] * ΔZ/2) / ((4/150) * ε^2 / (1 - ε)^2 * rₚ² * (P₀ / (μ[1] * v₀ * L)))
    P̅_right = 1
    WENO!(P̅_zf, P̅, P̅_left, P̅_right)

    T̅_zf = params.buffers["T̅_zf"]
    T̅_left = (T̅[1] + v̅_zf[1] * Peₕ * ΔZ/2) / (1 + v̅_zf[1] * Peₕ * ΔZ/2)
    T̅_right = T̅[N]
    WENO!(T̅_zf, T̅, T̅_left, T̅_right)
    
    #────────────────────────────────────────────
    # Differential equations (residual form: res = du - f(u))
    #────────────────────────────────────────────

    #────────────────────────────────────────────
    # Solid phase mass balance equation
    #────────────────────────────────────────────
    # H₂O
    α_H₂O        = k_H₂O * L / v₀
    @. res_x_H₂O = dx_H₂O - α_H₂O * (q_star_H₂O.(T₀ * T̅, P₀ * P̅ .* y_H₂O) / qₛ₀ - x_H₂O)

    # CO₂
    α_CO₂        = k_CO₂ * L / v₀
    @. res_x_CO₂ = dx_CO₂ - α_CO₂ * (q_star_CO₂.(T₀ * T̅, P₀ * P̅ .* y_CO₂, x_H₂O * qₛ₀) / qₛ₀ - x_CO₂)

    #────────────────────────────────────────────
    # Column energy balance equation
    #────────────────────────────────────────────
    T̅_flux = params.buffers["T̅_flux"]
    @. T̅_flux = v̅_zf * P̅_zf

    @inbounds for j in 1:N
        # BC for diffusion term
        if j == 1
            diffusion = 1/ΔZ * ((T̅[j+1] - T̅[j])/ΔZ - (-Peₕ * (T̅_feed - T̅_left)))
        elseif j == N
            diffusion = 1/ΔZ * (0 - (T̅[j] - T̅[j-1])/ΔZ)
        else
            diffusion = 1/ΔZ * ((T̅[j+1] - T̅[j])/ΔZ - (T̅[j] - T̅[j-1])/ΔZ)
        end

        Σdxⱼ   = dx_CO₂[j] + dx_H₂O[j]
        Σσⱼdxⱼ = σ_CO₂[j] * dx_CO₂[j] + σ_H₂O[j] * dx_H₂O[j]

        dT̅_computed = Ω₁[j] * diffusion +
                      - Ω₂[j] * 1/ΔZ * (T̅_flux[j+1] - T̅_flux[j]) +
                      - Ω₃[j] * T̅[j] * Σdxⱼ + 
                      + Σσⱼdxⱼ +
                      - Ω₄[j] * (T̅[j] - T̅_wall[j]) +
                      - Ω₂[j] * dP̅[j]
        
        res_T̅[j] = dT̅[j] - dT̅_computed
    end

    #────────────────────────────────────────────
    # Total mass balance equation
    #────────────────────────────────────────────
    P̅_flux = params.buffers["P̅_flux"]
    @. P̅_flux = (P̅_zf / T̅_zf) * v̅_zf
    @inbounds for j in 1:N
        Σdxⱼ   = dx_CO₂[j] + dx_H₂O[j] 
        dP̅_computed = - T̅[j] / ΔZ * (P̅_flux[j+1] - P̅_flux[j]) - ψ * T̅[j] * Σdxⱼ + P̅[j]/T̅[j] * dT̅[j]
        res_P̅[j] = dP̅[j] - dP̅_computed
    end

    #────────────────────────────────────────────
    # Component mass balance equations (only CO₂ and H₂O as differential)
    #────────────────────────────────────────────
    y_vars  = (y_CO₂, y_H₂O)
    res_y_vars = (res_y_CO₂, res_y_H₂O)
    dy_vars = (dy_CO₂, dy_H₂O)
    dx_vars = (dx_CO₂, dx_H₂O)
    y_feed_vals = (y_feed_CO₂, y_feed_H₂O)

    yᵢ_flux = params.buffers["y_flux"]

    @inbounds for i in 1:2  # 1:CO₂, 2:H₂O (N₂ is algebraic)
        yᵢ       = y_vars[i]
        res_yᵢ   = res_y_vars[i]
        dyᵢ      = dy_vars[i]
        dxᵢ      = dx_vars[i]
        yᵢ_feed  = y_feed_vals[i]

        # Compute yᵢ at cell faces using WENO
        yᵢ_zf    = similar(yᵢ, N + 1)
        yᵢ_left  = (yᵢ[1] + yᵢ_feed * v̅_zf[1] * Pe * ΔZ/2) / (1 + v̅_zf[1] * Pe * ΔZ/2)
        yᵢ_right = yᵢ[N]
        WENO!(yᵢ_zf, yᵢ, yᵢ_left, yᵢ_right)

        @. yᵢ_flux = yᵢ_zf * (P̅_zf / T̅_zf) * v̅_zf

        @inbounds for j in 1:N
            # BC for diffusion term
            if j == 1
                diffusion = 1/ΔZ * (P̅_zf[j+1] / T̅_zf[j+1] * (yᵢ[j + 1] - yᵢ[j]) / ΔZ 
                                                - P̅_zf[j] / T̅_zf[j] * (-Pe * (yᵢ_feed - yᵢ_left)))
            elseif j == N
                diffusion = 1/ΔZ * (0 - P̅_zf[j] / T̅_zf[j] * (yᵢ[j] - yᵢ[j - 1]) / ΔZ)
            else
                diffusion = 1/ΔZ * (P̅_zf[j+1] / T̅_zf[j+1] * (yᵢ[j + 1] - yᵢ[j]) / ΔZ 
                                                - P̅_zf[j] / T̅_zf[j] * (yᵢ[j] - yᵢ[j - 1]) / ΔZ)
            end

            dyᵢ_computed = 1/Pe * T̅[j]/P̅[j] * diffusion +
                          - T̅[j]/P̅[j] * 1/ΔZ * (yᵢ_flux[j+1] - yᵢ_flux[j]) + 
                          - ψ * T̅[j]/P̅[j] * dxᵢ[j] - yᵢ[j]/P̅[j] * dP̅[j] + yᵢ[j]/T̅[j] * dT̅[j]
            
            res_yᵢ[j] = dyᵢ[j] - dyᵢ_computed
        end
    end
    #────────────────────────────────────────────
    # Wall energy balance equation
    #────────────────────────────────────────────
    @inbounds for j in 1:N
        # BC for diffusion term
        if j == 1
            diffusion = 1/ΔZ * ((T̅_wall[j+1] - T̅_wall[j])/ΔZ - (T̅_wall[j] - T̅ₐ)/(ΔZ/2))
        elseif j == N
            diffusion = 1/ΔZ * ((T̅ₐ - T̅_wall[j])/(ΔZ/2) - (T̅_wall[j] - T̅_wall[j-1])/ΔZ)
        else
            diffusion = 1/ΔZ * ((T̅_wall[j+1] - T̅_wall[j])/ΔZ - (T̅_wall[j] - T̅_wall[j-1])/ΔZ)
        end

        dT̅_wall_computed = Π₁ * diffusion + Π₂ * (T̅[j] - T̅_wall[j]) - Π₃ * (T̅_wall[j] - T̅ₐ)
        res_T̅_wall[j] = dT̅_wall[j] - dT̅_wall_computed
    end

    return nothing
end

function setup_dae_initial_conditions(N, y0_CO₂, y0_H₂O, x0_CO₂, x0_H₂O, P̅0, T̅0, T̅_wall0, params)
    """
    Setup initial conditions for DAE formulation
    Returns u0, du0, and differential_vars vector
    Note: x_N₂ eliminated since it's zero everywhere
    """
    
    # State vector size
    n_diff = 7*N        # y_CO₂, y_H₂O, x_CO₂, x_H₂O, P̅, T̅, T̅_wall (x_N₂ eliminated)
    n_alg = N + (N+1)   # y_N₂, v̅
    n_total = n_diff + n_alg
    
    # Initialize state vector
    u0 = zeros(n_total)
    du0 = zeros(n_total)    
    
    # Differential variables
    u0[1:N] = y0_CO₂
    u0[N+1:2N] = y0_H₂O
    u0[2N+1:3N] = x0_CO₂
    u0[3N+1:4N] = x0_H₂O
    u0[4N+1:5N] = P̅0
    u0[5N+1:6N] = T̅0
    u0[6N+1:7N] = T̅_wall0
    
    # Algebraic variables (initial guess)
    u0[7N+1:8N] = 1.0 .- y0_CO₂ .- y0_H₂O  # y_N₂ = 1 - y_CO₂ - y_H₂O
    u0[8N+1] = 1.0
    u0[8N+2:8N+(N+1)] = fill(0.0, N)  # Initial velocity guess

    # Fix du0
    res = similar(u0)
    adsorption_dae!(res, du0, u0, params, 0)
    du0[4N+1] = -res[4N+1]
    adsorption_dae!(res, du0, u0, params, 0)
    du0[1] = - res[1]
    du0[N+1] = -res[N+1]
    du0[5N+1] = -res[5N+1]
    du0[4N+1] += -res[4N+1]

    
    # Differential_vars: true for differential, false for algebraic
    differential_vars = vcat(
        fill(true, n_diff),   # Differential variables
        fill(false, n_alg)    # Algebraic variables
    )
    
    return u0, du0, differential_vars
end

end