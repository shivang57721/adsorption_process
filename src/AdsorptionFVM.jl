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
const K_gas = 0.09
Base.@kwdef struct AdsorptionData
    k = - 1 / (150μ * (1 - ε_bed)^2 / (ε_bed^3 * dₚ^2))

    ic = 1
    iT = 2
    ip = 3

    Γ_in = 1
    Γ_out = 2
end

darcy_velocity(u, data) = data.k * R * (u[data.ip, 2] - u[data.ip, 1])

function flux(f, u, edge, data)
    # Flux for c is just v * c at edge
    vh = darcy_velocity(u, data)
    f[data.ic] = vh * u[data.ic]
    f[data.iT] = - ε_bed * K_gas * (u[data.iT, 2] - u[data.iT, 1]) + Cₚ_gas * vh * u[data.iT]
end

function storage(y, u, node, data)
    y[data.ic] = ε_total * u[data.ic]
    y[data.iT] = (ε_total * Cₚ_gas + ρ_gas * Cₚ_solid + ρ_gas * Cₚ_ads) * u[data.iT] - ε_total * u[data.ip]
end

# P = c R T
function reaction(y, u, node, data)
    y[data.ip] = u[data.ip] - u[data.ic] * R * u[data.iT]
end

function bcondition(y, u, node, data)
    boundary_dirichlet!(y, u, bnode; species = data.ic, region = data.Γ_in, value = P_high / (R * T_feed))
    boundary_dirichlet!(y, u, bnode; species = data.iT, region = data.Γ_in, value = T_feed)
    boundary_dirichlet!(y, u, bnode; species = data.ip, region = data.Γ_in, value = P_high)

    # Maybe change these conditions with boutflow conditions
    boundary_neumann!(y, u, bnode; species = data.ic, region = data.Γ_out, value = 0)
    boundary_neumann!(y, u, bnode; species = data.iT, region = data.Γ_out, value = 0)
    boundary_neumann!(y, u, bnode; species = data.ip, region = data.Γ_out, value = 0)
end

X = 0:0.1:1
data = AdsorptionData(;)
grid = VoronoiFVM.Grid(X)
sys = VoronoiFVM.System(
    grid;
    storage,
    flux,
    reaction,
    bcondition,
    data,
    species = [1,2,3]
)

sol = solve(sys)