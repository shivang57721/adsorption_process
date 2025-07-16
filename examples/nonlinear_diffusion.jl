"""
Solving the nonlinear diffusion equation ∂ₜu - Δuᵐ = 0
Also called the Porous Medium Problem
"""
using DifferentialEquations
using VoronoiFVM
using LinearAlgebra

# Exact solution to the problem is the Barenblatt solution (which has finite support for m>1)
function barenblatt(x, t, m)
    tx = t^(-1.0 / (m + 1.0))
    xx = x * tx
    xx = xx * xx
    xx = 1 - xx * (m - 1) / (2.0 * m * (m + 1))
    if xx < 0.0
        xx = 0.0
    end
    return tx * xx^(1.0 / (m - 1.0))
end

function create_porous_medium_problem(n, m)
    h = 1.0 / convert(Float64, n / 2)
    X = collect(-1:h:1)
    grid = VoronoiFVM.Grid(X)

    function flux!(f, u, edge, data)
        f[1] = u[1, 1]^m - u[1, 2]^m
    end

    storage!(f, u, node, data) = f[1] = u[1]

    sys = VoronoiFVM.System(grid, flux = flux!, storage = storage!, species = 1)
    return sys, X
end

diffeqmethods = Dict(
    "Rosenbrock23 (Rosenbrock)"    => Rosenbrock23,
    "QNDF2 (Like matlab's ode15s)" => QNDF2,
    "FBDF"                         => FBDF,
    "Implicit Euler"               => ImplicitEuler
)

function run_diffeq(; n = 20, m = 2, t0 = 0.001, tend = 0.01, solver = nothing)
    sys, X = create_porous_medium_problem(n, m)
    inival = unknowns(sys)
    inival[1, :] .= map(x -> barenblatt(x, t0, m), X)
    state = VoronoiFVM.SystemState(sys)
    problem = ODEProblem(state, inival, (t0, tend))
    odesol = solve(problem, solver)
    sol = reshape(odesol, sys; state)
    err = norm(sol[1, :, end] - map(x -> barenblatt(x, tend, m), X))
    return sol, sys, err
end

for method in diffeqmethods
    run_diffeq(m = 2, n = 10, solver = method.second()) # "Precompile"
end

method = "Rosenbrock23 (Rosenbrock)"
t2 = @elapsed sol2, sys2, err2 = run_diffeq(m = 2, n = 10, solver = diffeqmethods[method]())
history_summary(sol2)


"""
Example with changing the mass matrix
Consider the same equation, but now ∂u^(1/m) - Δu = 0
Expressed as a DAE 
    ∂ₜw - Δu = 0
    wᵐ - u = 0
"""
# Now u[1] = u, u[2] = w

# First equation has a storage term (but no reaction)
function dae_storage!(y, u, node, data)
    y[1] = u[2]
    return nothing
end

# Second equation has a reaction term (but no storage)
function dae_reaction!(y, u, node, data)
    y[2] = u[2]^m - u[1]
    return nothing
end

# First equation has a flux term
function flux!(f, u, edge, data)
    f[1] = u[1, 1]^m - u[1, 2]^m
end

dae_physics = VoronoiFVM.Physics(
    flux = flux!,
    storage = dae_storage!,
    reaction = dae_reaction!
)

n, m = 10, 2
h = 1.0 / convert(Float64, n / 2)
X = collect(-1:h:1)
grid = VoronoiFVM.Grid(X)
dae_sys = VoronoiFVM.System(grid, dae_physics, species = [1, 2])
dae_inival = unknowns(dae_sys)


