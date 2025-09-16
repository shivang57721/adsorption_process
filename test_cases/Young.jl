include("../src/AdsorptionModel.jl")

col_params = ColumnParams(
    Rᵢ = 0.04,
    Rₒ = 0.041,
    L = 0.01,
    h_L = 14,
    h_W = 22000
)

isotherm_params = (;
    ns0 = 4.86,
    chi = 0.0,      # Temperature dependence parameter for saturation capacity [dimensionless]
    T_ref = 298.15,    # Reference temperature [K]
    DH = -117798.0, # Adsorption enthalpy [J mol⁻¹]
    b0 = 2.85e-16,  # Affinity constant at infinite temperature [bar⁻¹]
    t0 = 0.209,     # Heterogeneity parameter at T0 [dimensionless]
    alfa =  0.523,  # Temperature dependence parameter for t [dimensionless]
    gamma = -0.137,         # CO₂-H₂O interaction parameter for saturation capacity [kg mol⁻¹]
    beta = 5.612           # CO₂-H₂O interaction parameter for affinity [kg mol⁻¹]
)

sorb_params = SorbentParams(
    ε_bed = 0.4,
    ε_total = 0.54,
    dₚ = 0.00052,
    ρ_bed = 528,
    ρ_particle = 880,
    k_CO2 = 0.003,
    k_H2O = 0.0086,
    ΔH_CO2 = -70000,
    ΔH_H2O = -46000,
    Dₘ = 1.3e-5,
    C_solid = 1580,
    q_star_CO2 = (T, R, p_CO2, q_H2O) -> Toth_isotherm_CO2_modified_H2O(T, R, p_CO2, q_H2O, isotherm_params),
    q_star_H2O = GAB_isotherm_H2O_Tfunction_Resins
)

RH=0.55
adsorption = OperatingParameters(;
                    step_name = Adsorption,
                    u_feed = 7.06*1e-2,
                    T_feed = 288.15,
                    T_amb = 288.15,
                    P_out = 1.01325e5,
                    y_CO2_feed = 0.0004,
                    y_H2O_feed = RH*Psat_H2O(288.15-273)/(1.01325e5),
                    duration = 8000)

heating = OperatingParameters(;
                step_name = Heating,
                u_feed = 0.0,
                T_feed = 288.15,
                T_amb = 373.15,
                P_out = 0.2e5,
                y_CO2_feed = 0.0,
                y_H2O_feed = 1.0,
                duration = 2400)

desorption = OperatingParameters(;
                step_name = Desorption,
                u_feed = 0.0,
                T_amb = 373.15,
                T_feed = 393.15,
                P_out = 0.2e5,
                y_CO2_feed = 0.0,
                y_H2O_feed = 1.0,
                duration = 20000)

cooling = OperatingParameters(;
            step_name = Cooling,
            u_feed = 0.0,
            T_amb = 288.15,
            T_feed = 120+273,
            P_out = 0.2e5,
            y_CO2_feed = 0.0,
            y_H2O_feed = 1.0,
            duration = 400)

pressurization = OperatingParameters(;
                    step_name = Pressurization,
                    u_feed = 0,
                    T_amb = 288.15,
                    T_feed = 288.15,
                    P_out = 1.01325e5,
                    y_CO2_feed = 0.0004,
                    y_H2O_feed = 0.0095,
                    duration = 30)

cycle_steps = [adsorption, heating, desorption, cooling, pressurization]
num_cycles = 4
@time sol, index_data = run_simulation(;N=10, cycle_steps, col_params, sorb_params, num_cycles)

# Compare plots with digitized data
using CSV
using DataFrames

duration_per_cycle = sum([s.duration for s in cycle_steps])
last_cycle = duration_per_cycle*(num_cycles - 1) : duration_per_cycle*num_cycles
df = CSV.File("test_cases/Young_2021/pressure.csv") |> DataFrame
p = plot(df[:, :x], df[:, " y"] .* 1e5, label="Young")
plot!([sol(t)[index_data.ip, end] for t in last_cycle], label="Julia Simulation")
xlabel!("Time (s)")
ylabel!("Pressure (Pa)")
savefig(p, "test_cases/Young_2021/figures/pressure.png")

df = CSV.File("test_cases/Young_2021/qco2.csv") |> DataFrame
p = plot(df[:, :x], df[:, " y"], label="Young")
plot!([sol(t)[index_data.iq_CO2, end] for t in last_cycle], label="Julia Simulation")
xlabel!("Time (s)")
ylabel!("q CO2")
savefig(p, "test_cases/Young_2021/figures/qco2.png")

df = CSV.File("test_cases/Young_2021/qh2o.csv") |> DataFrame
p = plot(df[:, :x], df[:, " y"], label="Young")
plot!([sol(t)[index_data.iq_H2O, end] for t in last_cycle], label="Julia Simulation")
xlabel!("Time (s)")
ylabel!("q H2O")
savefig(p, "test_cases/Young_2021/figures/qh2o.png")

df = CSV.File("test_cases/Young_2021/yco2.csv") |> DataFrame
p = plot(df[:, :x], df[:, " y"], label="Young")
plot!([sol(t)[index_data.iCO2, end] / (sol(t)[index_data.iN2, end] + sol(t)[index_data.iH2O, end] + sol(t)[index_data.iCO2, end]) for t in last_cycle], label="Julia Simulation")
xlabel!("Time (s)")
ylabel!("y CO2")
savefig(p, "test_cases/Young_2021/figures/yco2.png")

df = CSV.File("test_cases/Young_2021/T.csv") |> DataFrame
p = plot(df[:, :x], df[:, " y"], label="Young")
plot!([sol(t)[index_data.iT, end] for t in last_cycle], label="Julia Simulation")
xlabel!("Time (s)")
ylabel!("Temperature")
savefig(p, "test_cases/Young_2021/figures/T.png")

using XLSX
df = XLSX.readtable("test_cases/Python_sim_Young/x_sol_ads_cycle3.xlsx", 1) |> DataFrame
xs = [parse(Float64, t) for t in names(df)]
p = plot(xs, Vector(df[10, :]), label="Python simulation")
plot!([sol(t)[index_data.iCO2, end] for t in last_cycle[1]:(last_cycle[1] + adsorption.duration)], label="Julia simulation")
ylabel!("Concentration of CO2")
title!("Adsorption")
savefig(p, "test_cases/Python_sim_Young/figures/ads_co2.png")

df = XLSX.readtable("test_cases/Python_sim_Young/x_sol_heating_cycle3.xlsx", 1) |> DataFrame
xs = [parse(Float64, t) for t in names(df)]
p = plot(xs, Vector(df[10, :]), label="Python simulation")
ylabel!("Concentration of CO2")
title!("Heating")
plot!([sol(t)[index_data.iCO2, end] for t in (last_cycle[1] + adsorption.duration):(last_cycle[1] + adsorption.duration + heating.duration)], label="Julia simulation")
savefig(p, "test_cases/Python_sim_Young/figures/heating_co2.png")

df = XLSX.readtable("test_cases/Python_sim_Young/x_sol_des_cycle3.xlsx", 1) |> DataFrame
xs = [parse(Float64, t) for t in names(df)]
p = plot(xs, Vector(df[10, :]), label="Python simulation")
ylabel!("Concentration of CO2")
title!("Desorption")
plot!([sol(t)[index_data.iCO2, end] for t in (last_cycle[1] + adsorption.duration+heating.duration):(last_cycle[1] + adsorption.duration + heating.duration+desorption.duration)], label="Julia simulation")
savefig(p, "test_cases/Python_sim_Young/figures/des_co2.png")

df = XLSX.readtable("test_cases/Python_sim_Young/x_sol_cooling_cycle3.xlsx", 1) |> DataFrame
xs = [parse(Float64, t) for t in names(df)]
p = plot(xs, Vector(df[10, :]), label="Python simulation")
ylabel!("Concentration of CO2")
title!("Cooling")
plot!([sol(t)[index_data.iCO2, end] for t in (last_cycle[1] + adsorption.duration+heating.duration+desorption.duration):(last_cycle[1] + adsorption.duration + heating.duration+desorption.duration+cooling.duration)], label="Julia simulation")
savefig(p, "test_cases/Python_sim_Young/figures/cooling_co2.png")