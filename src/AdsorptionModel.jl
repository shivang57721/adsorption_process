"""
Semi-discretization of transport equations for adsorption modelling
"""

include("util.jl")

function adsorption_equations!(du, u, params, t)
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

    # Unpack buffers
    μ = params.buffers["μ"]
    ρ_gas = params.buffers["ρ_gas"]
    Σxᵢ = params.buffers["Σxᵢ"]

    Ω₁ = params.buffers["Ω₁"]
    Ω₂ = params.buffers["Ω₂"]
    Ω₃ = params.buffers["Ω₃"]
    Ω₄ = params.buffers["Ω₄"]
    σ_CO₂ = params.buffers["σ_CO₂"]
    σ_H₂O = params.buffers["σ_H₂O"]
    σ_N₂ = params.buffers["σ_N₂"]

    #────────────────────────────────────────────
    # Unpack u and du
    #────────────────────────────────────────────
    # Unpack u
    y_CO₂  = @view u[1:N]
    y_N₂   = @view u[N+1:2N]
    y_H₂O  = @view u[2N+1:3N]
    x_CO₂  = @view u[3N+1:4N]
    x_H₂O  = @view u[4N+1:5N]
    P̅      = @view u[5N+1:6N]
    T̅      = @view u[6N+1:7N]
    T̅_wall = @view u[7N+1:8N]

    # Unpack du
    dy_CO₂  = @view du[1:N]
    dy_N₂   = @view du[N+1:2N]
    dy_H₂O  = @view du[2N+1:3N]
    dx_CO₂  = @view du[3N+1:4N]
    dx_H₂O  = @view du[4N+1:5N]
    dP̅      = @view du[5N+1:6N]
    dT̅      = @view du[6N+1:7N]
    dT̅_wall = @view du[7N+1:8N]
    
    # μ and ρ_gas
    @. μ = y_CO₂ * μ_CO₂ + y_H₂O * μ_H₂O + y_N₂ * μ_N₂ # viscosity in each cell
    @. ρ_gas = P₀ / (R * T_feed) * (y_CO₂ * 44.009 + y_H₂O * 18.01528 + y_N₂ * 28.0134) * 1e-3 # gas density in each cell

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

    """ !info!
    _zf means value at the cell faces. zf[j] is the value at j - 1/2
    _left means boundary condition at Z=0
    _right means boundary condition at Z=1
    """

    # Compute local velcoity at cell boundaries
    v̅_zf = params.buffers["v̅_zf"]
    compute_velocity!(v̅_zf, P̅; y_CO₂, y_H₂O, y_N₂, params, t, P̅_right=1.0)

    # Compute P̅ and T̅ and y at cell faces using WENO
    P̅_zf = params.buffers["P̅_zf"]
    P̅_left = P̅[1] + (150μ[1] * (L * ΔZ/2)/(4rₚ²) * (1 - ε)^2 / ε^2 * 1.0 + 
                  + 1.75 * (L * ΔZ/2) * ρ_gas[1] / (2rₚ) * (1 - ε) / ε * 1.0 * abs(1.0)) / P₀
    P̅_right = 1
    WENO!(P̅_zf, P̅, P̅_left, P̅_right)

    T̅_zf = params.buffers["T̅_zf"]
    T̅_left = (T̅[1] + 1.0 * Peₕ * ΔZ/2) / (1 + 1.0 * Peₕ * ΔZ/2)
    T̅_right = T̅[N]
    WENO!(T̅_zf, T̅, T̅_left, T̅_right)

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

        dT̅[j]     = Ω₁[j] * diffusion +
                    - Ω₂[j] * 1/ΔZ * (T̅_flux[j+1] - T̅_flux[j]) +
                    - Ω₃[j] * T̅[j] * Σdxⱼ + 
                    + Σσⱼdxⱼ +
                    - Ω₄[j] * (T̅[j] - T̅_wall[j])
                    - Ω₂[j] * dP̅[j]
    end

    #────────────────────────────────────────────
    # Total mass balance equation
    #────────────────────────────────────────────
    P̅_flux = params.buffers["P̅_flux"]
    @. P̅_flux = (P̅_zf / T̅_zf) * v̅_zf
    @inbounds for j in 1:N
        Σdxⱼ   = dx_CO₂[j] + dx_H₂O[j]
        dP̅[j] = - T̅[j] / ΔZ * (P̅_flux[j+1] - P̅_flux[j]) - ψ * T̅[j] * Σdxⱼ + P̅[j]/T̅[j] * dT̅[j]
    end

    #────────────────────────────────────────────
    # Component mass balance equations
    #────────────────────────────────────────────
    y_vars  = (y_CO₂, y_N₂, y_H₂O)
    dy_vars = (dy_CO₂, dy_N₂, dy_H₂O)
    dx_vars = (dx_CO₂, nothing, dx_H₂O)  # dx_N₂ = 0, so use nothing as placeholder
    y_feed_vals = (y_feed_CO₂, y_feed_N₂, y_feed_H₂O)

    yᵢ_flux = params.buffers["y_flux"]

    @inbounds for i in 1:3  # 1:CO₂, 2:N₂, 3:H₂O
        yᵢ       = y_vars[i]
        dyᵢ      = dy_vars[i]
        dxᵢ      = dx_vars[i]
        yᵢ_feed  = y_feed_vals[i]

        # Compute yᵢ at cell faces using WENO
        yᵢ_zf    = params.buffers["y_zf"]
        yᵢ_left  = (yᵢ[1] + yᵢ_feed * Pe * ΔZ/2) / (1 + Pe * ΔZ/2)
        yᵢ_right = yᵢ[N]
        WENO!(yᵢ_zf, yᵢ, yᵢ_left, yᵢ_right; clamp_result=true)

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

            if i == 2  # N₂ component
                dyᵢ[j] = 1/Pe * T̅[j]/P̅[j] * diffusion +
                            - T̅[j]/P̅[j] * 1/ΔZ * (yᵢ_flux[j+1] - yᵢ_flux[j]) + 
                            - yᵢ[j]/P̅[j] * dP̅[j] + yᵢ[j]/T̅[j] * dT̅[j]
            else
                dyᵢ[j] = 1/Pe * T̅[j]/P̅[j] * diffusion +
                            - T̅[j]/P̅[j] * 1/ΔZ * (yᵢ_flux[j+1] - yᵢ_flux[j]) + 
                            - ψ * T̅[j]/P̅[j] * dxᵢ[j] - yᵢ[j]/P̅[j] * dP̅[j] + yᵢ[j]/T̅[j] * dT̅[j]
            end
        end
    end

    #────────────────────────────────────────────
    # Wall energy balance equation
    #────────────────────────────────────────────
    @inbounds for j in 1:N
        if j == 1
            diffusion = 1/ΔZ * ((T̅_wall[j+1] - T̅_wall[j])/ΔZ - (T̅_wall[j] - T̅ₐ)/(ΔZ/2))
        elseif j == N
            diffusion = 1/ΔZ * ((T̅ₐ - T̅_wall[j])/(ΔZ/2) - (T̅_wall[j] - T̅_wall[j-1])/ΔZ)
        else
            diffusion = 1/ΔZ * ((T̅_wall[j+1] - T̅_wall[j])/ΔZ - (T̅_wall[j] - T̅_wall[j-1])/ΔZ)
        end

        dT̅_wall[j] = Π₁ * diffusion + Π₂ * (T̅[j] - T̅_wall[j]) - Π₃ * (T̅_wall[j] - T̅ₐ)
    end
end