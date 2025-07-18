using Revise
using Test
using VoronoiFVM
using Plots

Base.@kwdef struct FlowTransportData
    k = 1.0
    v_in = 1.0
    c_in = 0.5
    D = 1.0
    Γ_in = 1
    Γ_out = 2
    ip = 1
    ic = 2
end

X = 0:0.1:1

darcyvelo(u, data) = data.k * (u[data.ip, 1] - u[data.ip, 2])

function flux(y, u, edge, data)
    vh = darcyvelo(u, data)
    y[data.ip] = vh

    bp, bm = fbernoulli_pm(vh / data.D)
    y[data.ic] = data.D * (bm * u[data.ic, 1] - bp * u[data.ic, 2])

    return nothing
end

function bcondition(y, u, bnode, data)
    boundary_dirichlet!(y, u, bnode; species = data.ip, region = data.Γ_in, value = 1.0)
    boundary_dirichlet!(y, u, bnode; species = data.ip, region = data.Γ_out, value = 0)

    boundary_dirichlet!(y, u, bnode; species = data.ic, region = data.Γ_in, value = data.c_in)
    boundary_neumann!(y, u, bnode; species = data.ic, region = data.Γ_out, value = 0)
    return nothing
end

# function boutflow(y, u, edge, data)
#     y[data.ic] = -darcyvelo(u, data) * u[data.ic, outflownode(edge)]
#     return nothing
# end

function flowtransportsystem(grid; kwargs...)
    data = FlowTransportData(; kwargs...)
    return VoronoiFVM.System(
        grid;
        flux,
        bcondition,
        data,
        species = [1, 2, 3],
    )
end

function checkinout(sys, sol)
    data = sys.physics.data
    tfact = TestFunctionFactory(sys)
    tf_in = testfunction(tfact, [data.Γ_out], [data.Γ_in])
    tf_out = testfunction(tfact, [data.Γ_in], [data.Γ_out])
    return (; in = integrate(sys, tf_in, sol), out = integrate(sys, tf_out, sol))
end

# 1D case
grid = VoronoiFVM.Grid(X)
sys1 = flowtransportsystem(grid);
sol1 = solve(sys1; verbose = "n");

plot(sol1[1, :])

t1 = checkinout(sys1, sol1)