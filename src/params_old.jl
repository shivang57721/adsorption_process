Base.@kwdef struct AdsorptionParams
    # Discretization parameters
    N::Int

    # Column geometry
    L::Float64
    r_in::Float64
    r_out::Float64
    ε::Float64         # Bed porosity
    εₚ::Float64
    rₚ::Float64         # Particle radius
    τ::Float64         # Tortuosity

    # Properties and constants
    ρₛ::Float64           # Adsorbent density
    ρ_wall::Float64
    Cₚ_gas::Float64
    Cₚ_ads::Float64
    Cₚₛ::Float64          # Specific heat of adsorbent
    Cₚ_wall::Float64
    Dₘ::Float64               # Molecular diffusivity
    γ::Float64                # Adiabatic constant
    K_wall::Float64
    K_gas::Float64
    h_in::Float64
    h_out::Float64
    R::Float64                   # Universal gas constant
    y_feed_N₂::Float64
    y_feed_CO₂::Float64
    y_feed_H₂O::Float64
    μ_N₂::Float64           # Fluid viscosity
    μ_CO₂::Float64
    μ_H₂O::Float64

    # Operating conditions
    v_feed::Float64
    T_feed::Float64
    Tₐ::Float64
    Pₕ::Float64
    Pᵢ::Float64
    Pₗ::Float64
    ΔH_CO₂::Float64
    ΔH_H₂O::Float64

    # Reference parameters for nondimensionalization
    qₛ₀::Float64

    # Isotherm functions
    q_star_CO₂::Function  # (T, p_CO2, q_H2O) → q_CO2^*
    q_star_H₂O::Function  # (T, p_H2O) → q_H2O^*
    k_CO₂::Float64
    k_H₂O::Float64

    buffers::Dict{String, Any} = Dict()
end

Base.@kwdef struct IsothermParams
    T₀::Float64
    nₛ₀::Float64
    b₀::Float64
    t₀::Float64
    ΔH₀::Float64
    χ::Float64
    γ_iso::Float64
    β::Float64
    c_G::Float64
    K_ads::Float64
    cₘ::Float64
    R::Float64
    α::Float64
end