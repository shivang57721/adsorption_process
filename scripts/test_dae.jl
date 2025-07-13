using Revise
include("../src/AdsorptionModelDAE.jl"); Revise.track("src/AdsorptionModelDAE.jl")
using .AdsorptionModelDAE

using DifferentialEquations
using Sundials
using Plots

# Load params file
include("../simulations/default_params.jl")

function setup_dae_problem(params)
    """
    Create a DAE problem from the existing ODE setup with consistent initial conditions
    """
    N = params.N

    # Set up initial conditions
    y0_CO₂  = fill(0.0, N)
    y0_H₂O  = fill(0.0, N)
    x0_CO₂  = fill(0.0, N)
    x0_H₂O  = fill(0.0, N)
    P̅0      = fill(1.0, N)
    T̅0      = fill(1.0, N)
    T̅_wall0 = fill(1.0, N)
    
    u0, du0, differential_vars = setup_dae_initial_conditions(
        N, y0_CO₂, y0_H₂O, x0_CO₂, x0_H₂O, P̅0, T̅0, T̅_wall0, params
    )
    
    # Time span
    tspan = (0.0, 200.0)
    
    # Create DAE problem
    prob = DAEProblem(adsorption_dae!, du0, u0, tspan, params; 
                      differential_vars=differential_vars)
    
    return prob
end

function extract_solution(sol, N)
    """
    Extract solution components from DAE solution
    Note: x_N₂ eliminated (implicitly zero everywhere)
    """
    # Differential variables
    y_CO₂  = t -> sol(t)[1:N]
    y_H₂O  = t -> sol(t)[N+1:2N]
    x_CO₂  = t -> sol(t)[2N+1:3N]
    x_H₂O  = t -> sol(t)[3N+1:4N]
    P̅      = t -> sol(t)[4N+1:5N]
    T̅      = t -> sol(t)[5N+1:6N]
    T̅_wall = t -> sol(t)[6N+1:7N]
    
    # Algebraic variables
    y_N₂   = t -> sol(t)[7N+1:8N]
    v̅_zf   = t -> sol(t)[8N+1:8N+(N+1)]
    
    # x_N₂ is implicitly zero everywhere
    x_N₂   = t -> zeros(N)
    
    return (y_CO₂=y_CO₂, y_H₂O=y_H₂O, y_N₂=y_N₂, 
            x_CO₂=x_CO₂, x_N₂=x_N₂, x_H₂O=x_H₂O,
            P̅=P̅, T̅=T̅, T̅_wall=T̅_wall, v̅_zf=v̅_zf)
end

function test_dae_residuals(params)
    """
    Test the DAE residuals at initial conditions to debug issues
    """
    println("Testing DAE residuals at initial conditions...")
    
    N = params.N
    
    # Set up initial conditions
    y0_CO₂  = fill(0.0, N)
    y0_H₂O  = fill(0.0, N)
    x0_CO₂  = fill(0.0, N)
    x0_H₂O  = fill(0.0, N)
    P̅0      = fill(1.0, N)
    T̅0      = fill(1.0, N)
    T̅_wall0 = fill(1.0, N)
    
    u0, du0, differential_vars = setup_dae_initial_conditions(
        N, y0_CO₂, y0_H₂O, x0_CO₂, x0_H₂O, P̅0, T̅0, T̅_wall0, params
    )
    
    # Test residual calculation
    res = similar(u0)
    try
        adsorption_dae!(res, du0, u0, params, 0.0)
        
        # Check residuals for different components
        y_CO₂_res = res[1:N]
        y_H₂O_res = res[N+1:2N]
        x_CO₂_res = res[2N+1:3N]
        x_H₂O_res = res[3N+1:4N]
        P̅_res = res[4N+1:5N]
        T̅_res = res[5N+1:6N]
        T̅_wall_res = res[6N+1:7N]
        y_N₂_res = res[7N+1:8N]
        v̅_res = res[8N+1:8N+(N+1)]
        
        println("Max residuals:")
        println("  y_CO₂: ", maximum(abs.(y_CO₂_res)))
        println("  y_H₂O: ", maximum(abs.(y_H₂O_res)))
        println("  x_CO₂: ", maximum(abs.(x_CO₂_res)))
        println("  x_H₂O: ", maximum(abs.(x_H₂O_res)))
        println("  P̅: ", maximum(abs.(P̅_res)))
        println("  T̅: ", maximum(abs.(T̅_res)))
        println("  T̅_wall: ", maximum(abs.(T̅_wall_res)))
        println("  y_N₂ (algebraic): ", maximum(abs.(y_N₂_res)))
        println("  v̅ (algebraic): ", maximum(abs.(v̅_res)))
        
        overall_max = maximum(abs.(res))
        println("Overall max residual: $overall_max")
        
        if overall_max < 1e-10
            println("✓ Initial conditions appear consistent!")
        else
            println("⚠ Initial conditions may need adjustment")
        end
        
    catch e
        println("❌ Error evaluating DAE residuals: $e")
        rethrow(e)
    end
end

# Example usage:
if !@isdefined(TESTING) || !TESTING
    # First test the residuals at initial conditions
    test_dae_residuals(params)
    
    println("\nSetting up and solving DAE problem...")
    prob = setup_dae_problem(params)
    
    @time sol = solve(prob,
    DABDF2(autodiff=AutoFiniteDiff()),
        abstol = 1e-6,
        reltol = 1e-5,
        dtmin = 1e-20,
        saveat = 0.1,
        verbose = true
    )
    
    println("Extracting results...")
    N = params.N
    results = extract_solution(sol, N)
    
    # Plot results
    t_final = sol.t[end]
    
    plot(results.P̅(t_final), title="Pressure", xlabel="Cell", ylabel="P̅")
    plot!(results.T̅(t_final), title="Temperature", xlabel="Cell", ylabel="T̅")
    plot!(results.y_CO₂(t_final), title="CO₂ concentration", xlabel="Cell", ylabel="y_CO₂")
    plot!(results.y_H₂O(t_final), title="H₂O concentration", xlabel="Cell", ylabel="y_H₂O")
    plot!(results.y_N₂(t_final), title="N₂ concentration", xlabel="Cell", ylabel="y_N₂")
    
    # Check normalization constraint
    y_sum = results.y_CO₂(t_final) + results.y_H₂O(t_final) + results.y_N₂(t_final)
    println("Max deviation from y_sum = 1: ", maximum(abs.(y_sum .- 1.0)))
end
