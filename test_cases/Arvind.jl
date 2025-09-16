include("../src/AdsorptionModel.jl")

# 1. Column Parameters
col_params = ColumnParams(
    Rᵢ = 0.04,
    Rₒ = 0.041,
    L = 0.01,
    h_L = 3.0, # from heat_transfer_coeff_hL
    h_W = 26.0 # from heat_transfer_coeff_wall
)

isotherm_params = (;
    ns0 = 2.2,
    chi = 0,
    T_ref = 296,
    DH =  -6.0e4,
    b0 = 3.73e4,
    t0 = 0.4247,
    alfa = -0.4921,
    gamma = 0.00958,
    beta = 3.448,

    Cg = 0.1489,
    K_ads = 0.5751,
    Cm = 36.48
)

# 2. Sorbent Parameters
sorb_params = SorbentParams(
    ε_bed = 0.092,
    ε_total = 0.9616,
    dₚ = 0.0075,
    ρ_bed = 55.0,
    ρ_particle = 61.0,
    k_CO2 = 2.0e-4,
    k_H2O = 2.0e-1,
    ΔH_CO2 = -60000.0,
    ΔH_H2O = -43800.0,
    Dₘ = 2.05e-5,
    C_solid = 2070.0,
    # The functions are already set by the struct's default definition
    q_star_CO2 = (T, R, p_CO2, q_H2O) -> Toth_isotherm_CO2_modified_H2O_Stampi(T,R,p_CO2,q_H2O,isotherm_params),
    q_star_H2O = (T, p_H2O) -> GAB_isotherm_H2O(T, p_H2O, isotherm_params)
)

adsorption = OperatingParameters(;
                    step_name = Adsorption,
                    u_feed = 0.01,
                    T_feed = 20+273,
                    T_amb = 20+273,
                    P_out = 1e5,
                    y_CO2_feed = 0.0004,
                    y_H2O_feed = 0.0115,
                    duration = 13772)

heating = OperatingParameters(;
                step_name = Heating,
                u_feed = 0.0,
                T_feed = 120+273.15,
                T_amb = 120+273.15,
                P_out = 0.05e5,
                y_CO2_feed = 0.0,
                y_H2O_feed = 1.0,
                duration = 348)

desorption = OperatingParameters(;
                step_name = Desorption,
                u_feed = 0.005,
                T_amb = 120+273.15,
                T_feed = 120+273.15,
                P_out = 0.05e5,
                y_CO2_feed = 0.0,
                y_H2O_feed = 1.0,
                duration = 30000)

cooling = OperatingParameters(;
            step_name = Cooling,
            u_feed = 0.0,
            T_feed = 20+273,
            T_amb = 20+273,
            P_out = 0.05e5,
            y_CO2_feed = 0.0,
            y_H2O_feed = 1.0,
            duration = 800)

pressurization = OperatingParameters(;
                    step_name = Pressurization,
                    u_feed = 0,
                    T_feed = 20+273,
                    T_amb = 20+273,
                    P_out = 1e5,
                    y_CO2_feed = 0.0004,
                    y_H2O_feed = 0.0115,
                    duration = 60)

cycle_steps = [adsorption, heating, desorption, cooling, pressurization]
num_cycles = 4
@time sol, index_data = run_simulation(;N=10, cycle_steps, col_params, sorb_params, num_cycles)

# Compare plots with digitized data
using CSV
using DataFrames

duration_per_cycle = sum([s.duration for s in cycle_steps])
last_cycle = duration_per_cycle*(num_cycles - 1) : duration_per_cycle*num_cycles

df = CSV.File("test_cases/Arvind_2023/qco2.csv") |> DataFrame
p = plot(df[:, :x], df[:, " y"], label="Arvind")
plot!([sol(t)[index_data.iq_CO2, end] for t in last_cycle], label="Julia Simulation")
xlabel!("Time (s)")
ylabel!("q CO2")
savefig(p, "test_cases/Arvind_2023/figures/qco2.png")

df = CSV.File("test_cases/Arvind_2023/qh2o.csv") |> DataFrame
p = plot(df[:, :x], df[:, " y"], label="Arvind")
plot!([sol(t)[index_data.iq_H2O, end] for t in last_cycle], label="Julia Simulation")
xlabel!("Time (s)")
ylabel!("q H2O")
savefig(p, "test_cases/Arvind_2023/figures/qh2o.png")

df = CSV.File("test_cases/Arvind_2023/T.csv") |> DataFrame
p = plot(df[:, :x], df[:, " y"], label="Arvind")
plot!([sol(t)[index_data.iT, end] for t in last_cycle], label="Julia Simulation")
xlabel!("Time (s)")
ylabel!("Temperature")
savefig(p, "test_cases/Arvind_2023/figures/T.png")
