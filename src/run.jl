include("../src/AdsorptionModel.jl")

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
                    T_amb = 372.5,
                    P_out = 1e5,
                    y_CO2_feed = 0.0,
                    y_H2O_feed = 0.0,
                    duration = 3600*2)

desorption = OperatingParameters(;
                step_name = "Desorption",
                u_feed = 0.1,
                T_feed = 372.5,
                T_amb = 372.5,
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

cycle_steps = [adsorption, heating, desorption]
@time sol, index_data = run_simulation(;N=10, cycle_steps, num_cycles=1)

ts = 1:3600*6
p = plot(ts./3600, [sol(t)[index_data.ip, 1] for t in ts])

# sol(3600*6)