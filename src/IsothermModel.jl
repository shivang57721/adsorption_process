module IsothermModel

export IsothermParams

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

function q_star_H₂O(T, P_H₂O, params)
    c_G     = params.c_G
    K_ads   = params.K_ads
    cₘ      = params.cₘ
    
    # Calculate relative humidity from T and P_H₂O
    Pₛ(temp) = 611.21 * exp((18.678 - temp / 234.5) * temp / (temp + 273.15 - 16.01))
    x = P_H₂O / Pₛ(T - 273)

    return cₘ * (c_G * K_ads * x) / ((1 - K_ads * x) * (1 + (c_G - 1) * K_ads * x))
end

function q_star_CO₂(T, P_CO₂, q_H₂O, params)
    # Isotherm parameters
    T₀      = params.T₀
    nₛ₀     = params.nₛ₀
    b₀      = params.b₀
    t₀      = params.t₀
    ΔH₀     = params.ΔH₀
    χ       = params.χ
    γ_iso   = params.γ_iso
    β       = params.β
    R       = params.R
    α       = params.α

    nₛ(T) = nₛ₀ * exp(χ * (1 - T/T₀)) * (1 / (1 - γ_iso * q_H₂O))
    b(T) = b₀ * exp(ΔH₀/(R * T₀) * (T₀/T - 1)) * (1 + β * q_H₂O)
    t(T) = t₀ + α * (1 - T₀/T)

    base = max(b(T) * P_CO₂, 0)
    val = nₛ(T) * (b(T) * P_CO₂) / (1 + base^t(T))^(1/t(T))
    return isnan(val) ? 0 : val
end

end