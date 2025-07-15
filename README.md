# adsorption_process
Numerical simulation of adsorption.
Run `scripts/run_adsorption.jl` to run the simulation, where the equations are written as an ODE system.
Run `scripts/run_dae.jl` to run the simulation, where the equations are written as a DAE system.
Currently adsorption is "turned off", by setting the constants k_CO₂ = 0 and k_H₂O = 0 in `simulations/default_params.jl`.