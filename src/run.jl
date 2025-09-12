include("../src/AdsorptionModel.jl")

adsorption = OperatingParameters(;
                    step_name = Adsorption,
                    u_feed = 7.06*1e-2,
                    T_feed = 288.15,
                    T_amb = 288.15,
                    P_out = 1.01325e5,
                    y_CO2_feed = 0.0004,
                    y_H2O_feed = 0.0095,
                    duration = 8000)

heating = OperatingParameters(;
                step_name = Heating,
                u_feed = 0.0,
                T_feed = 288.15,
                T_amb = 99.35+273.15,
                P_out = 0.12e5,
                y_CO2_feed = 0.0,
                y_H2O_feed = 1.0,
                duration = 1147)

desorption = OperatingParameters(;
                step_name = Desorption,
                u_feed = 0.0,
                T_amb = 373.15,
                T_feed = 120+273,
                P_out = 0.12e5,
                y_CO2_feed = 0.0,
                y_H2O_feed = 1.0,
                duration = 21377)

cooling = OperatingParameters(;
            step_name = Cooling,
            u_feed = 0.0,
            T_amb = 288.15,
            T_feed = 120+273,
            P_out = 0.12e5,
            y_CO2_feed = 0.0,
            y_H2O_feed = 1.0,
            duration = 400)

pressurization = OperatingParameters(;
                    step_name = Pressurization,
                    u_feed = 0.0,
                    T_amb = 288.15,
                    T_feed = 288.15,
                    P_out = 1.01325e5,
                    y_CO2_feed = 0.0004,
                    y_H2O_feed = 0.0095,
                    duration = 60)

pressurization_reset = OperatingParameters(;
                    step_name = PressurizationReset,
                    u_feed = 0.0,
                    T_amb = 288.15,
                    T_feed = 288.15,
                    P_out = 1.01325e5,
                    y_CO2_feed = 0.0004,
                    y_H2O_feed = 0.0095,
                    duration = 0.0)

cycle_steps = [adsorption, heating, desorption, cooling, pressurization]
num_cycles = 4
@time sol, index_data = run_simulation(;N=10, cycle_steps, num_cycles)


total_time = sum([s.duration for s in cycle_steps]) * num_cycles
ts = 1:total_time
# plot(ts, [(sol(t)[index_data.iCO2, end] + sol(t)[index_data.iH2O, end] + sol(t)[index_data.iN2, end]) for t in ts])
plot(ts, [sol(t)[index_data.iq_CO2, end] for t in ts])

# plot(sol(8000+1147+21377+400)[index_data.ip, :])