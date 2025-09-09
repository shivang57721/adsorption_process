"""
Simulation of adsorption equations with VoronoiFVM.jl
"""

using DifferentialEquations
using VoronoiFVM
using LinearAlgebra
using Plots
include("params.jl")

phys_params = PhysicalParams()
col_params = ColumnParams()
sorb_params = SorbentParams()
op_params = OperatingParams()

# Unpack phys_params
const R = phys_params.R
const γ₁ = phys_params.γ₁
const γ₂ = phys_params.γ₂
const Dₘ = phys_params.Dₘ
const μ = phys_params.μ
const C_solid = phys_params.C_solid
const Cₚ_CO2 = phys_params.Cₚ_CO2
const Cₚ_H2O = phys_params.Cₚ_H2O
const Cₚ_N2 = phys_params.Cₚ_N2
const Cₚ_wall = phys_params.Cₚ_wall
const MW_CO2 = phys_params.MW_CO2
const MW_H2O = phys_params.MW_H2O
const MW_N2  = phys_params.MW_N2

# Unpack col_params
const Rᵢ = col_params.Rᵢ
const Rₒ = col_params.Rₒ
const L = col_params.L
const h_L = col_params.h_L
const h_W = col_params.h_W

# Unpack sorb_params
const ε_bed = sorb_params.ε_bed
const ε_total = sorb_params.ε_total
const dₚ = sorb_params.dₚ
const ρ_bed = sorb_params.ρ_bed
const ρ_particle = sorb_params.ρ_particle
const k_CO2 = sorb_params.k_CO2
const k_H2O = sorb_params.k_H2O
const ΔH_CO2 = sorb_params.ΔH_CO2
const ΔH_H2O = sorb_params.ΔH_H2O

# CO2 isotherm params
const qinf0_dry_CO2 = sorb_params.qinf0_dry_CO2
const chi_dry_CO2 = sorb_params.chi_dry_CO2
const T0_dry_CO2 = sorb_params.T0_dry_CO2
const DH_dry_CO2 = sorb_params.DH_dry_CO2
const b0_dry_CO2 = sorb_params.b0_dry_CO2
const t0_dry_CO2 = sorb_params.t0_dry_CO2
const alfa_dry_CO2 = sorb_params.alfa_dry_CO2

const qinf0_wet_CO2 = sorb_params.qinf0_wet_CO2
const chi_wet_CO2 = sorb_params.chi_wet_CO2
const T0_wet_CO2 = sorb_params.T0_wet_CO2
const DH_wet_CO2 = sorb_params.DH_wet_CO2
const b0_wet_CO2 = sorb_params.b0_wet_CO2
const t0_wet_CO2 = sorb_params.t0_wet_CO2
const alfa_wet_CO2 = sorb_params.alfa_wet_CO2
const A = sorb_params.A
const gamma = sorb_params.gamma
const beta = sorb_params.beta

# GAB H2O isotherm params
const qm = sorb_params.qm
const C = sorb_params.C
const D = sorb_params.D
const F = sorb_params.F
const G = sorb_params.G

# Unpack op_params
const u_feed = op_params.u_feed
const T_amb = op_params.T_amb
const T_feed = op_params.T_feed
const P_out = op_params.P_out
const y_CO2_feed = op_params.y_CO2_feed
const y_H2O_feed = op_params.y_H2O_feed
const y_N2_feed = op_params.y_N2_feed

# Other parameters than depend on the previous ones
const ρ_gas  = P_out / (R * T_feed) * (y_CO2_feed * MW_CO2 + y_H2O_feed * MW_H2O + y_N2_feed * MW_N2) * 1e-3
const a_wall =  π * (Rₒ^2 - Rᵢ^2)
const D_L    = γ₁ * Dₘ + γ₂ * dₚ * u_feed / ε_bed

const c_total_feed = P_out / (R * T_feed)
const c_N2_feed  = y_N2_feed * c_total_feed
const c_CO2_feed = y_CO2_feed * c_total_feed
const c_H2O_feed = y_H2O_feed * c_total_feed
const C_gas_feed = c_CO2_feed * Cₚ_CO2 + c_H2O_feed * Cₚ_H2O + c_N2_feed * Cₚ_N2
const K_L = D_L * C_gas_feed

darcy_velocity(u, data) = begin
    k = - 150μ * (1 - ε_bed)^2 / (ε_bed^3 * dₚ^2)
    1/k * (u[data.ip, 2] - u[data.ip, 1])
end

function flux_exponential(f, u, edge, data)
    vh = darcy_velocity(u, data)

    # --- Calculate effective dispersion and Péclet number ---
    c_total_1 = u[data.iN2, 1] + u[data.iCO2, 1] + u[data.iH2O, 1]
    c_total_2 = u[data.iN2, 2] + u[data.iCO2, 2] + u[data.iH2O, 2]
    c_total_edge = (c_total_1 + c_total_2) / 2

    # --- Compute Mole fractions ---
    y_N2_1  = u[data.iN2, 1]  / c_total_1; y_N2_2  = u[data.iN2, 2]  / c_total_2
    y_CO2_1 = u[data.iCO2, 1] / c_total_1; y_CO2_2 = u[data.iCO2, 2] / c_total_2
    y_H2O_1 = u[data.iH2O, 1] / c_total_1; y_H2O_2 = u[data.iH2O, 2] / c_total_2

    bp, bm = fbernoulli_pm(vh / (ε_bed * D_L))
    f[data.iN2] = ε_bed * D_L * c_total_edge * (bm * y_N2_1 - bp * y_N2_2)
    f[data.iCO2] = ε_bed * D_L * c_total_edge * (bm * y_CO2_1 - bp * y_CO2_2)
    f[data.iH2O] = ε_bed * D_L * c_total_edge * (bm * y_H2O_1 - bp * y_H2O_2)

    # --- Thermal Flux (Energy Balance) ---
    C_gas_1 = u[data.iCO2,1]*Cₚ_CO2 + u[data.iH2O,1]*Cₚ_H2O + u[data.iN2,1]*Cₚ_N2
    C_gas_2 = u[data.iCO2,2]*Cₚ_CO2 + u[data.iH2O,2]*Cₚ_H2O + u[data.iN2,2]*Cₚ_N2
    C_gas_edge = (C_gas_1 + C_gas_2) / 2

    T_1     = u[data.iT, 1];   T_2     = u[data.iT, 2]
    D = ε_bed * K_L
    bp, bm = fbernoulli_pm(vh * C_gas_edge / D)
    f[data.iT] = D * (bm * T_1 - bp * T_2)
end

function storage(y, u, node, data)
    y[data.iN2]  = ε_total * u[data.iN2]
    y[data.iCO2] = ε_total * u[data.iCO2]
    y[data.iH2O] = ε_total * u[data.iH2O]

    C_gas = u[data.iCO2] * Cₚ_CO2 + u[data.iH2O] * Cₚ_H2O + u[data.iN2] * Cₚ_N2
    C_ads = u[data.iq_CO2] * Cₚ_CO2 + u[data.iq_H2O] * Cₚ_H2O
    y[data.iT]   = (ε_total * C_gas + ρ_bed * C_solid + ρ_bed * C_ads) * u[data.iT] - ε_total * u[data.ip]
    y[data.iT_wall] = u[data.iT_wall]

    y[data.iq_CO2] = u[data.iq_CO2]
    y[data.iq_H2O] = u[data.iq_H2O]
end

function reaction(y, u, node, data)
    # P = (c_N2 + c_CO2 + c_H2O) * R * T
    c_total_node = u[data.iN2] + u[data.iCO2] + u[data.iH2O]
    y[data.ip] = u[data.ip] - c_total_node * R * u[data.iT]

    # T_wall term
    y[data.iT] = 2h_L/Rᵢ * (u[data.iT] - u[data.iT_wall])
    y[data.iT_wall] = - 2π/(Cₚ_wall * a_wall) * (h_L * Rᵢ * (u[data.iT] - u[data.iT_wall]) - h_W * Rₒ * (u[data.iT_wall] - T_amb))

    # q term
    p_H2O = u[data.ip] * u[data.iH2O] / c_total_node
    p_CO2 = u[data.ip] * u[data.iCO2] / c_total_node
    qstar_H2O = q_star_H2O(qm, C, D, F, G, u[data.iT], p_H2O)
    qstar_CO2 = q_star_CO2(qinf0_dry_CO2, chi_dry_CO2, T0_dry_CO2, b0_dry_CO2, DH_dry_CO2, t0_dry_CO2, alfa_dry_CO2, u[data.iT], R, p_CO2, gamma, beta, qstar_H2O)

    y[data.iq_H2O] = - k_H2O * (qstar_H2O - u[data.iq_H2O])
    y[data.iq_CO2] = - k_CO2 * (qstar_CO2 - u[data.iq_CO2])

    y[data.iH2O] = ρ_bed * k_H2O * (qstar_H2O - u[data.iq_H2O])
    y[data.iCO2] = ρ_bed * k_CO2 * (qstar_CO2 - u[data.iq_CO2])

    y[data.iT] += ρ_bed * ΔH_CO2 * k_CO2 * (qstar_CO2 - u[data.iq_CO2]) + ρ_bed * ΔH_H2O * k_H2O * (qstar_H2O - u[data.iq_H2O])
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

Base.@kwdef struct AdsorptionData
    iN2     = 1
    iCO2    = 2
    iH2O    = 3
    iT      = 4
    ip      = 5
    iT_wall = 6
    iq_CO2  = 7
    iq_H2O  = 8

    Γ_in = 1
    Γ_out = 2
end

function simulation(;storage=storage, flux=flux_exponential, reaction=reaction, bcondition=bcondition, boutflow=boutflow, N=10, tads=3600)
    X = range(0, L, N)
    data = AdsorptionData(;)
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
        species = [1,2,3,4,5,6,7,8]
    )

    inival = unknowns(sys)

    inival[data.ip, :]      .= P_out
    inival[data.iT, :]      .= T_amb
    inival[data.iT_wall, :] .= T_amb
    inival[data.iN2, :]     .= c_total_feed # Initially 100% N2
    inival[data.iCO2, :]    .= 0.0
    inival[data.iH2O, :]    .= 0.0
    inival[data.iq_CO2, :]  .= 0.0
    inival[data.iq_H2O, :]  .= 0.0

    state = VoronoiFVM.SystemState(sys)
    problem = ODEProblem(state, inival, (0, tads))
    odesol = solve(problem, Rodas5P())
    sol(t) = reshape(odesol(t), sys)

    return sol, data
end

# t = 3600
# plot(sol(t)[data.iCO2, :], title="CO2 concentration in column at time=$(t)s")

# ts = 1:3600*2

# plot(ts, [sol(t)[data.iCO2, end] / (sol(t)[data.iCO2, end] + sol(t)[data.iN2, end] +sol(t)[data.iH2O, end])  for t in ts], title="Concentration of CO2 at the end of the column")
# xlabel!("Time (s)")