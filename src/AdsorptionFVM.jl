"""
Simulation of adsorption equations with VoronoiFVM.jl
"""

const ε_total  = 0.742
const ε_bed    = 0.403
const Cₚ_gas   = 30.7
const Cₚ_solid = 1070
const Cₚ_ads   = 30.7

const R        = 8.314
const T_feed   = 298.15
const P_high   = 1e5
const ρ_gas    = P_high / (R * T_feed) * (0.0004 * 44.009 + 0.1 * 18.01528 + 0.8996 * 28.0134) * 1e-3

const μ = 1.789e-5
const rₚ = 1e-3
const dₚ = 2rₚ

darcy_velocity(u, data) = data.k * (u[data.ip, 1] - u[data.ip, 2])

Base.@kwdef struct AdsorptionData
    k = 1 / (150μ * (1 - ε_bed)^2 / (ε_bed^3 * dₚ^2))

    ip = 1
    ic = 2
    iT = 3
end

function flux!(f, u, edge, data)
    # Flux for c is just v * c at edge
    vh = darcy_velocity(u, data)
    f[data.ic] = vh * ()
    

end
