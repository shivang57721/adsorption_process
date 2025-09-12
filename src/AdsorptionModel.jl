"""
Simulation of adsorption equations with VoronoiFVM.jl
"""

using DifferentialEquations
using Sundials
using VoronoiFVM
using LinearAlgebra
using Plots
include("params.jl")

phys_consts = PhysicalConstants()
col_params = ColumnParams()
sorb_params = SorbentParams()

# Unpack phys_consts
const R = phys_consts.R
const γ₁ = phys_consts.γ₁
const γ₂ = phys_consts.γ₂
const μ = phys_consts.μ
const Cₚ_CO2 = phys_consts.Cₚ_CO2
const Cₚ_H2O = phys_consts.Cₚ_H2O
const Cₚ_N2 = phys_consts.Cₚ_N2
const Cₚ_wall = phys_consts.Cₚ_wall
const MW_CO2 = phys_consts.MW_CO2
const MW_H2O = phys_consts.MW_H2O
const MW_N2  = phys_consts.MW_N2

# Unpack col_params
const Rᵢ = col_params.Rᵢ
const Rₒ = col_params.Rₒ
const L = col_params.L
const h_L = col_params.h_L
const h_W = col_params.h_W
const a_wall =  π * (Rₒ^2 - Rᵢ^2)

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
const Dₘ = sorb_params.Dₘ
const C_solid = sorb_params.C_solid

q_star_H2O = sorb_params.q_star_H2O
q_star_CO2 = sorb_params.q_star_CO2

darcy_velocity(u, data) = begin
    k = - 150μ * (1 - ε_bed)^2 / (ε_bed^3 * dₚ^2)
    1/k * (u[data.ip, 2] - u[data.ip, 1])
end

ergun_velocity(u, data) = begin
    params = data.params
    ρ_gas = params.P_out / (R * params.T_amb) * (params.y_CO2_feed * MW_CO2 + params.y_H2O_feed * MW_H2O + params.y_N2_feed * MW_N2) * 1e-3
    b = 150μ * (1 - ε_bed)/(dₚ * 1.75 * ρ_gas)
    c = ε_bed^ 3 * dₚ / (1.75 * ρ_gas * (1 - ε_bed))
    dpdz = u[data.ip, 1] - u[data.ip, 2]
    sign(dpdz) * (-b + sqrt(b^2 + sign(dpdz)* 4c * dpdz))/2
end

function flux_exponential(f, u, edge, data)
    params = data.params
    vh = ergun_velocity(u, data)

    # --- Calculate effective dispersion and Péclet number ---
    c_total_1 = u[data.iN2, 1] + u[data.iCO2, 1] + u[data.iH2O, 1]
    c_total_2 = u[data.iN2, 2] + u[data.iCO2, 2] + u[data.iH2O, 2]
    c_total_edge = (c_total_1 + c_total_2) / 2

    # --- Compute Mole fractions ---
    y_N2_1  = u[data.iN2, 1]  / c_total_1; y_N2_2  = u[data.iN2, 2]  / c_total_2
    y_CO2_1 = u[data.iCO2, 1] / c_total_1; y_CO2_2 = u[data.iCO2, 2] / c_total_2
    y_H2O_1 = u[data.iH2O, 1] / c_total_1; y_H2O_2 = u[data.iH2O, 2] / c_total_2

    bp, bm = fbernoulli_pm(vh / (ε_bed * params.D_L))
    f[data.iN2] = ε_bed * params.D_L * c_total_edge * (bm * y_N2_1 - bp * y_N2_2)
    f[data.iCO2] = ε_bed * params.D_L * c_total_edge * (bm * y_CO2_1 - bp * y_CO2_2)
    f[data.iH2O] = ε_bed * params.D_L * c_total_edge * (bm * y_H2O_1 - bp * y_H2O_2)

    # --- Thermal Flux (Energy Balance) ---
    C_gas_1 = u[data.iCO2,1]*Cₚ_CO2 + u[data.iH2O,1]*Cₚ_H2O + u[data.iN2,1]*Cₚ_N2
    C_gas_2 = u[data.iCO2,2]*Cₚ_CO2 + u[data.iH2O,2]*Cₚ_H2O + u[data.iN2,2]*Cₚ_N2
    C_gas_edge = (C_gas_1 + C_gas_2) / 2

    T_1     = u[data.iT, 1];   T_2     = u[data.iT, 2]
    D = ε_bed * params.K_L
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
    params = data.params
    # P = (c_N2 + c_CO2 + c_H2O) * R * T
    c_total_node = u[data.iN2] + u[data.iCO2] + u[data.iH2O]
    y[data.ip] = u[data.ip] - c_total_node * R * u[data.iT]

    # T_wall term
    y[data.iT] = 2h_L/Rᵢ * (u[data.iT] - u[data.iT_wall])
    y[data.iT_wall] = - 2π/(Cₚ_wall * a_wall) * (h_L * Rᵢ * (u[data.iT] - u[data.iT_wall]) - h_W * Rₒ * (u[data.iT_wall] - params.T_amb))

    # q term
    p_H2O = u[data.ip] * u[data.iH2O] / c_total_node
    p_CO2 = u[data.ip] * u[data.iCO2] / c_total_node
    qstar_H2O = q_star_H2O(u[data.iT], p_H2O)
    qstar_CO2 = q_star_CO2(u[data.iT], R, p_CO2, qstar_H2O)

    y[data.iq_H2O] = - k_H2O * (qstar_H2O - u[data.iq_H2O])
    y[data.iq_CO2] = - k_CO2 * (qstar_CO2 - u[data.iq_CO2])

    y[data.iH2O] = ρ_bed * k_H2O * (qstar_H2O - u[data.iq_H2O])
    y[data.iCO2] = ρ_bed * k_CO2 * (qstar_CO2 - u[data.iq_CO2])

    y[data.iT] += ρ_bed * ΔH_CO2 * k_CO2 * (qstar_CO2 - u[data.iq_CO2]) + ρ_bed * ΔH_H2O * k_H2O * (qstar_H2O - u[data.iq_H2O])
end

function bcondition(y, u, bnode, data)
    params = data.params
    
    # Inlet feed flow boundary conditions
    boundary_neumann!(y, u, bnode; species=data.iN2, region=data.Γ_in, value = params.u_feed * params.c_N2_feed)
    boundary_neumann!(y, u, bnode; species=data.iCO2, region=data.Γ_in, value = params.u_feed * params.c_CO2_feed)
    boundary_neumann!(y, u, bnode; species=data.iH2O, region=data.Γ_in, value = params.u_feed * params.c_H2O_feed)
    boundary_neumann!(y, u, bnode; species=data.iT, region=data.Γ_in, value = params.u_feed * params.C_gas_feed * params.T_feed)

    # Dirichlet boundary conditions for Pressure
    if params.step_name == Pressurization
        boundary_dirichlet!(y, u, bnode; species=data.ip, region=data.Γ_in, value = params.P_out_func(bnode.time))
    elseif params.step_name != Cooling
        boundary_dirichlet!(y, u, bnode; species=data.ip, region=data.Γ_out, value = params.P_out_func(bnode.time))
    end
end

# Outflow boundary conditions
function boutflow(y, u, edge, data)
    params = data.params
    if params.step_name == Cooling
        # Both sides are closed
        return nothing
    end
    
    if params.step_name == Pressurization
        # Left boundary is open
        if outflownode(edge) == data.Γ_in
            vh = ergun_velocity(u, data)
            y[data.iN2] = -vh * params.c_N2_feed
            y[data.iCO2]  = -vh * params.c_CO2_feed
            y[data.iH2O]  = -vh * params.c_H2O_feed

            y[data.iT]   = -vh * params.C_gas_feed * params.T_feed
        end
    else
        # Right boundary is open
        if outflownode(edge) == data.Γ_out
            vh = ergun_velocity(u, data)
            y[data.iN2]   = -vh * u[data.iN2, outflownode(edge)]
            y[data.iCO2]  = -vh * u[data.iCO2, outflownode(edge)]
            y[data.iH2O]  = -vh * u[data.iH2O, outflownode(edge)]

            C_gas = u[data.iCO2, outflownode(edge)] * Cₚ_CO2 + 
                    u[data.iH2O, outflownode(edge)] * Cₚ_H2O + 
                    u[data.iN2, outflownode(edge)] * Cₚ_N2

            y[data.iT]   = -vh * C_gas * u[data.iT, outflownode(edge)]
        end
    end
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

    params::OperatingParameters
end

function run_simulation(; N=10, cycle_steps, num_cycles=1)
    for step in cycle_steps
        update_derived_params!(step, phys_consts, sorb_params)
    end

    # --- System Initialization (only done once) ---
    data = AdsorptionData(; params=deepcopy(cycle_steps[1]))
    X = range(0, L, length=N)
    grid = VoronoiFVM.Grid(X)
    sys = VoronoiFVM.System(
        grid;
        storage,
        flux=flux_exponential,
        reaction,
        bcondition,
        boutflow,
        data,
        outflowboundaries=[data.Γ_in, data.Γ_out],
        species=[1, 2, 3, 4, 5, 6, 7, 8]
    )

    # --- Initial Conditions for the very first cycle ---
    inival = unknowns(sys)
    inival[data.ip, :]      .= cycle_steps[1].P_out
    inival[data.iT, :]      .= cycle_steps[1].T_feed
    inival[data.iT_wall, :] .= cycle_steps[1].T_amb
    inival[data.iN2, :]     .= cycle_steps[1].c_total_feed
    inival[data.iCO2, :]    .= 0.0
    inival[data.iH2O, :]    .= 0.0
    inival[data.iq_CO2, :]  .= 0.0
    inival[data.iq_H2O, :]  .= 0.0

    # --- Simulation Loop for Multiple Cycles ---
    solutions = []
    u_current = inival

    all_steps = repeat(cycle_steps, num_cycles)

    for step in all_steps
        println("Running step $(step.step_name)")
        if step.step_name == PressurizationReset # Doesn't work currently
            u_current[data.ip, :]      .= step.P_out
            u_current[data.iCO2, :]    .= step.c_CO2_feed
            u_current[data.iH2O, :]    .= step.c_H2O_feed
            u_current[data.iN2, :]     .= step.c_N2_feed
            u_current[data.iT, :]      .= step.T_feed
            u_current[data.iT_wall, :] .= step.T_feed
            continue
        end

        step.P_out_func = t -> step.P_out + (u_current[data.ip, end] - step.P_out) * exp(- 0.11 * t)
        copy_params!(data.params, step)

        problem = ODEProblem(VoronoiFVM.SystemState(sys), u_current, (0, step.duration))
        # odesol = solve(problem, Rodas5P(), reltol=1e-5, abstol=1e-5)
        odesol = solve(problem, FBDF())#, reltol = 1e-8, abstol = 1e-8)

        push!(solutions, odesol)
        
        # Update state for the next step
        u_current = odesol.u[end]
    end

    # --- Generalized Solution Interpolation Function ---
    # Pre-calculate the end time of each step
    end_times = cumsum([s.duration for s in all_steps])
    start_times = start_times = [0; end_times[1:end-1]]
    
    function sol(t)
        # Find the index of the first step that ends after time t
        idx = findfirst(x -> t <= x, end_times)
        return reshape(solutions[idx](t - start_times[idx]), sys)
    end

    return sol, data
end

# t = 3600
# plot(sol(t)[data.iCO2, :], title="CO2 concentration in column at time=$(t)s")

# ts = 1:3600*2

# plot(ts, [sol(t)[data.iCO2, end] / (sol(t)[data.iCO2, end] + sol(t)[data.iN2, end] +sol(t)[data.iH2O, end])  for t in ts], title="Concentration of CO2 at the end of the column")
# xlabel!("Time (s)")