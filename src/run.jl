include("../src/AdsorptionModel.jl")

using DifferentialEquations
@show 0.55 * Psat_H2O(288.15 - 273) / 1.01325e5
adsorption = OperatingParameters(;
                    step_name = "Adsorption",
                    u_feed = 0.1,
                    T_feed = 288.15,
                    T_amb = 288.15,
                    y_CO2_feed = 0.0004,
                    y_H2O_feed = 0.01,
                    P_out = 1e5,
                    duration = 3600*2)

heating = OperatingParameters(;
                    step_name = "Heating",
                    u_feed = 0.0,
                    T_feed = 378.17,
                    T_amb = 99.35+273.15,
                    P_out = 1e5,
                    y_CO2_feed = 0.0,
                    y_H2O_feed = 0.0,
                    duration = 3600*2)

desorption = OperatingParameters(;
                step_name = "Desorption",
                u_feed = 0.1,
                T_feed = 378.17,
                T_amb = 99.35+273.15,
                P_out = 1e5,
                y_CO2_feed = 0.0,
                y_H2O_feed = 0.0,
                duration = 3600*2)

cooling = OperatingParameters(;
                    step_name = "Cooling",
                    u_feed = 0.0,
                    T_feed = 378.17,
                    T_amb = 288.15,
                    P_out = 1e5,
                    y_CO2_feed = 0.0,
                    y_H2O_feed = 0.0,
                    duration = 3600*2)

cycle_steps = [adsorption, heating, desorption, cooling]
@time sol, index_data = run_simulation(;N=10, cycle_steps, num_cycles=1)

ts = 1:3600*8*1
p = plot(ts./3600, [sol(t)[index_data.iT, end] for t in ts])