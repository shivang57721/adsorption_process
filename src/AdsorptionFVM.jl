"""
Simulation of adsorption equations with VoronoiFVM.jl
"""
const ε_total  = 0.742
const ε_bed    = 0.403
const C_solid  = 1000
const Cₚ_CO2   = 37.1
const Cₚ_H2O   = 33.6
const Cₚ_N2    = 29.1
const C_wall   = 4.0e6

# Define feed conditions for clarity
const R        = 8.314
const T_amb    = 298.15
const P_out    = 1e5
const ρ_gas    = P_out / (R * T_amb) * (0.0004 * 44.009 + 0.1 * 18.01528 + 0.8996 * 28.0134) * 1e-3
const ρ_bed    = 507

const u_feed = 1.0
const c_total_feed = P_out / (R * T_amb)
const c_N2_feed = 0.86 * c_total_feed
const c_CO2_feed = 0.04 * c_total_feed
const c_H2O_feed = 0.10 * c_total_feed
const C_gas_feed = c_CO2_feed * Cₚ_CO2 + c_H2O_feed * Cₚ_H2O + c_N2_feed * Cₚ_N2

const Dₘ  = 4.3e-6
const μ = 1.789e-5
const dₚ = 0.003

# Experimental params
const Rₒ = 0.02
const Rᵢ = 0.0125
const h_L = 236
const h_W = 220
const a_wall = π * (Rₒ^2 - Rᵢ^2)

using DifferentialEquations
using VoronoiFVM
using LinearAlgebra
using Plots
Base.@kwdef struct AdsorptionData
    iN2     = 1
    iCO2    = 2
    iH2O    = 3
    iT      = 4
    ip      = 5
    iT_wall = 6

    Γ_in = 1
    Γ_out = 2
end

darcy_velocity(u, data) = begin
    k = 150μ * (1 - ε_bed)^2 / (ε_bed^3 * dₚ^2)
    1/k * (u[data.ip, 1] - u[data.ip, 2])
end

ergun_velocity(u, data) = begin
    b = 150μ * (1 - ε_bed)/(dₚ * 1.75 * ρ_gas)
    c = ε_bed^ 3 * dₚ / (1.75 * ρ_gas * (1 - ε_bed))
    dpdz = u[data.ip, 1] - u[data.ip, 2]
    sign(dpdz) * (-b + sqrt(b^2 + sign(dpdz)* 4c * dpdz))/2
end

D_L(u; γ₁=0.7, γ₂=0.5) = γ₁ * Dₘ + γ₂ * dₚ * u / ε_bed
K_L(u, C_gas) = D_L(u) * C_gas

function flux_upwind(f, u, edge, data)
    vh = ergun_velocity(u, data)

    # Concentrations at nodes 1 and 2
    c_N2_1  = u[data.iN2, 1];  c_N2_2  = u[data.iN2, 2]
    c_CO2_1 = u[data.iCO2, 1]; c_CO2_2 = u[data.iCO2, 2]
    c_H2O_1 = u[data.iH2O, 1]; c_H2O_2 = u[data.iH2O, 2]
    T_1     = u[data.iT, 1];   T_2     = u[data.iT, 2]

    # == UPWINDING SCHEME START ==
    # Determine the upwind concentrations based on flow direction (vh)
    c_N2_upwind, c_CO2_upwind, c_H2O_upwind, T_upwind = if vh > 0
        (c_N2_1, c_CO2_1, c_H2O_1, T_1) # Flow is from node 1 to 2
    else
        (c_N2_2, c_CO2_2, c_H2O_2, T_2) # Flow is from node 2 to 1
    end
    # == UPWINDING SCHEME END ==

    # Total concentrations for dispersion term
    c_total_1 = u[data.iN2, 1] + u[data.iCO2, 1] + u[data.iH2O, 1]
    c_total_2 = u[data.iN2, 2] + u[data.iCO2, 2] + u[data.iH2O, 2]
    c_total_edge = (c_total_1 + c_total_2) / 2
    
    # Mole fractions for dispersion term
    y_N2_1  = u[data.iN2, 1]  / c_total_1; y_N2_2  = u[data.iN2, 2]  / c_total_2
    y_CO2_1 = u[data.iCO2, 1] / c_total_1; y_CO2_2 = u[data.iCO2, 2] / c_total_2
    y_H2O_1 = u[data.iH2O, 1] / c_total_1; y_H2O_2 = u[data.iH2O, 2] / c_total_2

    # Calculate fluxes for each component using UPWINDED values for advection
    f[data.iN2]  = vh * c_N2_upwind  - ε_bed * D_L(vh) * c_total_edge * (y_N2_2  - y_N2_1)
    f[data.iCO2] = vh * c_CO2_upwind - ε_bed * D_L(vh) * c_total_edge * (y_CO2_2 - y_CO2_1)
    f[data.iH2O] = vh * c_H2O_upwind - ε_bed * D_L(vh) * c_total_edge * (y_H2O_2 - y_H2O_1)
    
    C_gas_upwind = c_CO2_upwind * Cₚ_CO2 + c_H2O_upwind * Cₚ_H2O + c_N2_upwind * Cₚ_N2
    C_gas_edge = (c_total_1 * (u[data.iCO2,1]*Cₚ_CO2 + u[data.iH2O,1]*Cₚ_H2O + u[data.iN2,1]*Cₚ_N2) + c_total_2 * (u[data.iCO2,2]*Cₚ_CO2 + u[data.iH2O,2]*Cₚ_H2O + u[data.iN2,2]*Cₚ_N2)) / (2*c_total_edge) # Weighted average
    
    f[data.iT] = vh * C_gas_upwind * T_upwind + ε_bed * K_L(vh, C_gas_edge) * (T_1 - T_2)
end

function storage(y, u, node, data)
    y[data.iN2]  = ε_total * u[data.iN2]
    y[data.iCO2] = ε_total * u[data.iCO2]
    y[data.iH2O] = ε_total * u[data.iH2O]
    C_gas = u[data.iCO2] * Cₚ_CO2 + u[data.iH2O] * Cₚ_H2O + u[data.iN2] * Cₚ_N2
    y[data.iT]   = (ε_total * C_gas + ρ_bed * C_solid) * u[data.iT] - ε_total * u[data.ip]
    y[data.iT_wall] = u[data.iT_wall]
end

function reaction(y, u, node, data)
    # P = (c_N2 + c_CO2 + c_H2O) * R * T
    c_total_node = u[data.iN2] + u[data.iCO2] + u[data.iH2O]
    y[data.ip] = u[data.ip] - c_total_node * R * u[data.iT]

    # T_wall term
    y[data.iT] = 2h_L/Rᵢ * (u[data.iT] - u[data.iT_wall])
    y[data.iT_wall] = - 2π/(C_wall * a_wall) * (h_L * Rᵢ * (u[data.iT] - u[data.iT_wall]) - h_W * Rₒ * (u[data.iT_wall] - T_amb))
end

function bcondition(y, u, bnode, data)
    # Boundary conditions at z=0
    boundary_neumann!(y, u, bnode; species=data.iN2, region=data.Γ_in, value = u_feed * c_N2_feed)
    boundary_neumann!(y, u, bnode; species=data.iCO2, region=data.Γ_in, value = u_feed * c_CO2_feed)
    boundary_neumann!(y, u, bnode; species=data.iH2O, region=data.Γ_in, value = u_feed * c_H2O_feed)
    boundary_neumann!(y, u, bnode; species=data.iT, region=data.Γ_in, value = u_feed * C_gas_feed * T_amb)

    # boundary conditions at z=1
    boundary_dirichlet!(y, u, bnode; species=data.ip, region=data.Γ_out, value = P_out)
end

# Outflow boundary conditions
function boutflow(y, u, edge, data)
    vh = ergun_velocity(u, data)
    y[data.iN2]   = -vh * u[data.iN2, outflownode(edge)]
    y[data.iCO2]  = -vh * u[data.iCO2, outflownode(edge)]
    y[data.iH2O]  = -vh * u[data.iH2O, outflownode(edge)]

    C_gas = u[data.iCO2, outflownode(edge)] * Cₚ_CO2 + 
            u[data.iH2O, outflownode(edge)] * Cₚ_H2O + 
            u[data.iN2, outflownode(edge)] * Cₚ_N2

    y[data.iT]   = -vh * C_gas * u[data.iT, outflownode(edge)]
end

X = 0:0.005:1
data = AdsorptionData(;)
grid = VoronoiFVM.Grid(X)
sys = VoronoiFVM.System(
    grid;
    storage,
    flux = flux_upwind,
    reaction,
    bcondition,
    boutflow,
    data,
    outflowboundaries = [data.Γ_out],
    species = [1,2,3,4,5,6]
)

diffeqmethods = Dict(
    "Rosenbrock23 (Rosenbrock)"    => Rosenbrock23,
    "QNDF2 (Like matlab's ode15s)" => QNDF2,
    "FBDF"                         => FBDF,
    "Implicit Euler"               => ImplicitEuler
)

inival = unknowns(sys)

inival[data.ip, :]      .= P_out
inival[data.iT, :]      .= T_amb
inival[data.iT_wall, :] .= T_amb
inival[data.iN2, :]     .= c_total_feed # Initially 100% N2
inival[data.iCO2, :]    .= 0.0
inival[data.iH2O, :]    .= 0.0

state = VoronoiFVM.SystemState(sys)
problem = ODEProblem(state, inival, (0, 2))
odesol = solve(problem, FBDF())
sol = reshape(odesol, sys; state)

idx = 150
plot(sol.u[idx][data.iCO2, :], title=sol.t[idx])

plot(sol.u[idx][data.iCO2, :] .+ sol.u[idx][data.iH2O, :] .+ sol.u[idx][data.iN2, :])