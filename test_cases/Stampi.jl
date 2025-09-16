include("../src/AdsorptionModel.jl")

# 1. Column Parameters
col_params = ColumnParams(
    Rᵢ = 0.04,
    Rₒ = 0.041,
    L = 0.01,
    h_L = 3.0,    # from heat_transfer_coeff_hL
    h_W = 26.0    # from heat_transfer_coeff_wall
)

isotherm_params = (;
    # --- Toth for CO2 ---
    ns0 = 2.38,      # [mol/Kg]
    chi = 0,
    T_ref = 296,
    b0 = 7.07e3,    # [1/bar]
    t0 = 0.4148,    # [-]
    DH = -5.7047e4, # [J/mol]
    alfa = -1.606,

    # --- Toth Water Interaction ---
    gamma = 0.0016, # [kg/mol]
    beta = 59.1,    # [kg/mol]

    # --- GAB for H2O ---
    Cg = 0.1489,        # [-]
    K_ads = 0.5751,     # [-]
    Cm = 36.48         # [mol/kg]
)

# 2. Sorbent Parameters
sorb_params = SorbentParams(
    ε_bed = 0.092,
    ε_total = 0.9,
    dₚ = 0.0075,
    ρ_bed = 55.4,
    ρ_particle = 61.0,
    k_CO2 = 2.0e-4,
    k_H2O = 2.0e-3,
    ΔH_CO2 = -57000.0,
    ΔH_H2O = -49000.0,
    Dₘ = 4.3e-6,
    C_solid = 2070.0,

    q_star_H2O = (T, p_H2O) -> GAB_isotherm_H2O(T, p_H2O, isotherm_params),
    q_star_CO2 = (T, R, p_CO2, q_H2O) -> Toth_isotherm_CO2_modified_H2O_Stampi(T, R, p_CO2, q_H2O, isotherm_params)
)

cross_section = π * col_params.Rᵢ^2

adsorption = OperatingParameters(;
                    step_name = Adsorption,
                    u_feed = 50e-6/cross_section,
                    T_feed = 20+273,
                    T_amb = 20+273,
                    P_out = 1e5,
                    y_CO2_feed = 0.0004,
                    y_H2O_feed = 0.0115,
                    duration = 13772)

heating = OperatingParameters(;
                step_name = Heating,
                u_feed = 0.0,
                T_feed = 288.15,
                T_amb = 95+273,
                P_out = 50e2,
                y_CO2_feed = 0.0,
                y_H2O_feed = 1.0,
                duration = 3600,
                ΔT = 5)

desorption = OperatingParameters(;
                step_name = Desorption,
                u_feed = 25e-6/cross_section,
                T_amb = 95+273,
                T_feed = 95+273,
                P_out = 50e2,
                y_CO2_feed = 0.0,
                y_H2O_feed = 1.0,
                duration = 30000)

pressurization = OperatingParameters(;
                    step_name = Pressurization,
                    u_feed = 0,
                    T_feed = 20+273,
                    T_amb = 20+273,
                    P_out = 1e5,
                    y_CO2_feed = 0.0004,
                    y_H2O_feed = 0.0115,
                    duration = 30)

cycle_steps = [adsorption, heating, desorption, pressurization]
num_cycles = 4
@time sol, index_data = run_simulation(;N=10, cycle_steps, col_params, sorb_params, num_cycles)

# Compare plots with digitized data
using CSV
using DataFrames

duration_per_cycle = sum([s.duration for s in cycle_steps])
last_cycle = duration_per_cycle*(num_cycles - 1) : duration_per_cycle*num_cycles
df = CSV.File("test_cases/Stampi_Bombelli_2020/pressure.csv") |> DataFrame
p = plot(df[:, :x], df[:, " y"] .* 1e5, label="Stampi")
plot!([sol(t)[index_data.ip, end] for t in last_cycle], label="Julia Simulation")
xlabel!("Time (s)")
ylabel!("Pressure (Pa)")
savefig(p, "test_cases/Stampi_Bombelli_2020/figures/pressure.png")

df = CSV.File("test_cases/Stampi_Bombelli_2020/qco2.csv") |> DataFrame
p = plot(df[:, :x], df[:, " y"], label="Stampi")
plot!([sol(t)[index_data.iq_CO2, end] for t in last_cycle], label="Julia Simulation")
xlabel!("Time (s)")
ylabel!("q CO2")
savefig(p, "test_cases/Stampi_Bombelli_2020/figures/qco2.png")

df = CSV.File("test_cases/Stampi_Bombelli_2020/qh2o.csv") |> DataFrame
p = plot(df[:, :x], df[:, " y"], label="Stampi")
plot!([sol(t)[index_data.iq_H2O, end] for t in last_cycle], label="Julia Simulation")
xlabel!("Time (s)")
ylabel!("q H2O")
savefig(p, "test_cases/Stampi_Bombelli_2020/figures/qh2o.png")

df = CSV.File("test_cases/Stampi_Bombelli_2020/T.csv") |> DataFrame
p = plot(df[:, :x], df[:, " y"], label="Stampi")
plot!([sol(t)[index_data.iT, end] for t in last_cycle], label="Julia Simulation")
xlabel!("Time (s)")
ylabel!("Temperature")
savefig(p, "test_cases/Stampi_Bombelli_2020/figures/T.png")
