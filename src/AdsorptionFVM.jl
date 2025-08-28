"""
Simulation of adsorption equations with VoronoiFVM.jl
"""
const ε_total  = 0.742
const ε_bed    = 0.403
const C_solid  = 1000
const Cₚ_CO2   = 37.1
const Cₚ_H2O   = 33.6
const Cₚ_N2    = 29.1

# Define feed conditions for clarity
const R        = 8.314
const T_feed   = 298.15
const P_out    = 1e5
const ρ_gas    = P_out / (R * T_feed) * (0.0004 * 44.009 + 0.1 * 18.01528 + 0.8996 * 28.0134) * 1e-3
const ρ_bed    = 507

const u_feed = 1.0
const c_total_feed = P_out / (R * T_feed)
const c_N2_feed = 0.86 * c_total_feed
const c_CO2_feed = 0.04 * c_total_feed
const c_H2O_feed = 0.10 * c_total_feed
const C_gas_feed = c_CO2_feed * Cₚ_CO2 + c_H2O_feed * Cₚ_H2O + c_N2_feed * Cₚ_N2

const Dₘ  = 4.3e-6
const μ = 1.789e-5
const dₚ = 0.003

using DifferentialEquations
using VoronoiFVM
using LinearAlgebra
using Plots
Base.@kwdef struct AdsorptionData
    k = 150μ * (1 - ε_bed)^2 / (ε_bed^3 * dₚ^2)

    iN2  = 1
    iCO2 = 2
    iH2O = 3
    iT   = 4
    ip   = 5

    Γ_in = 1
    Γ_out = 2
end

darcy_velocity(u, data) = 1/data.k * (u[data.ip, 1] - u[data.ip, 2])
D_L(u; γ₁=0.7, γ₂=0.5) = γ₁ * Dₘ + γ₂ * dₚ * u / ε_bed
K_L(u, C_gas) = D_L(u) * C_gas

function flux(f, u, edge, data)
    vh = darcy_velocity(u, data)

    # Concentrations at the edge
    c_N2_edge  = (u[data.iN2, 1]  + u[data.iN2, 2])  / 2
    c_CO2_edge = (u[data.iCO2, 1] + u[data.iCO2, 2]) / 2
    c_H2O_edge = (u[data.iH2O, 1] + u[data.iH2O, 2]) / 2
    T_edge     = (u[data.iT, 1]   + u[data.iT, 2])   / 2

    # Total concentrations at nodes 1 and 2 of the edge
    c_total_1 = u[data.iN2, 1] + u[data.iCO2, 1] + u[data.iH2O, 1]
    c_total_2 = u[data.iN2, 2] + u[data.iCO2, 2] + u[data.iH2O, 2]
    # Total concentration at the edge
    c_total_edge = (c_total_1 + c_total_2) / 2
    
    # Mole fractions (y = c_i / c_total) at nodes 1 and 2
    y_N2_1  = u[data.iN2, 1]  / c_total_1; y_N2_2  = u[data.iN2, 2]  / c_total_2
    y_CO2_1 = u[data.iCO2, 1] / c_total_1; y_CO2_2 = u[data.iCO2, 2] / c_total_2
    y_H2O_1 = u[data.iH2O, 1] / c_total_1; y_H2O_2 = u[data.iH2O, 2] / c_total_2

    # Calculate fluxes for each component
    # Advection part + Dispersion part
    f[data.iN2]  = vh * c_N2_edge  - ε_bed * D_L(vh) * c_total_edge * (y_N2_2  - y_N2_1)
    f[data.iCO2] = vh * c_CO2_edge - ε_bed * D_L(vh) * c_total_edge * (y_CO2_2 - y_CO2_1)
    f[data.iH2O] = vh * c_H2O_edge - ε_bed * D_L(vh) * c_total_edge * (y_H2O_2 - y_H2O_1)
    
    C_gas = c_CO2_edge * Cₚ_CO2 + c_H2O_edge * Cₚ_H2O + c_N2_edge * Cₚ_N2
    f[data.iT] = C_gas * vh * T_edge + ε_bed * K_L(vh, C_gas) * (u[data.iT, 1] - u[data.iT, 2])
end

function storage(y, u, node, data)
    y[data.iN2]  = ε_total * u[data.iN2]
    y[data.iCO2] = ε_total * u[data.iCO2]
    y[data.iH2O] = ε_total * u[data.iH2O]
    C_gas = u[data.iCO2] * Cₚ_CO2 + u[data.iH2O] * Cₚ_H2O + u[data.iN2] * Cₚ_N2
    y[data.iT]   = (ε_total * C_gas + ρ_bed * C_solid) * u[data.iT] - ε_total * u[data.ip]
end

function reaction(y, u, node, data)
    # P = (c_N2 + c_CO2 + c_H2O) * R * T
    c_total_node = u[data.iN2] + u[data.iCO2] + u[data.iH2O]
    y[data.ip] = u[data.ip] - c_total_node * R * u[data.iT]
end

function bcondition(y, u, bnode, data)
    # Boundary conditions at z=0
    boundary_neumann!(y, u, bnode; species=data.iN2, region=data.Γ_in, value = u_feed * c_N2_feed)
    boundary_neumann!(y, u, bnode; species=data.iCO2, region=data.Γ_in, value = u_feed * c_CO2_feed)
    boundary_neumann!(y, u, bnode; species=data.iH2O, region=data.Γ_in, value = u_feed * c_H2O_feed)
    boundary_neumann!(y, u, bnode; species=data.iT, region=data.Γ_in, value = u_feed * C_gas_feed * T_feed)

    # boundary conditions at z=1
    boundary_dirichlet!(y, u, bnode; species=data.ip, region=data.Γ_out, value = P_out)
end

# Outflow boundary conditions
function boutflow(y, u, edge, data)
    vh = darcy_velocity(u, data)
    y[data.iN2]   = -vh * u[data.iN2, outflownode(edge)]
    y[data.iCO2]  = -vh * u[data.iCO2, outflownode(edge)]
    y[data.iH2O]  = -vh * u[data.iH2O, outflownode(edge)]

    C_gas = u[data.iCO2, outflownode(edge)] * Cₚ_CO2 + 
            u[data.iH2O, outflownode(edge)] * Cₚ_H2O + 
            u[data.iN2, outflownode(edge)] * Cₚ_N2

    y[data.iT]   = -vh * C_gas * u[data.iT, outflownode(edge)]
end

X = 0:0.01:1
data = AdsorptionData()
grid = VoronoiFVM.Grid(X)
sys = VoronoiFVM.System(
    grid;
    storage,
    flux,
    reaction,
    bcondition,
    boutflow,
    data,
    outflowboundaries = [data.Γ_out],
    species = [1,2,3,4,5]
)

diffeqmethods = Dict(
    "Rosenbrock23 (Rosenbrock)"    => Rosenbrock23,
    "QNDF2 (Like matlab's ode15s)" => QNDF2,
    "FBDF"                         => FBDF,
    "Implicit Euler"               => ImplicitEuler
)

inival = unknowns(sys)

inival[data.ip, :]   .= P_out
inival[data.iT, :]   .= T_feed
inival[data.iN2, :]  .= c_total_feed # Initially 100% N2
inival[data.iCO2, :] .= 0.0
inival[data.iH2O, :] .= 0.0

state = VoronoiFVM.SystemState(sys)
problem = ODEProblem(state, inival, (0, 2))
odesol = solve(problem, FBDF())
sol = reshape(odesol, sys; state)

idx = 300
plot(sol.u[idx][data.iT, :], title=sol.t[idx])

plot(sol.u[idx][data.iCO2, :] .+ sol.u[idx][data.iH2O, :] .+ sol.u[idx][data.iN2, :])