using Printf
using VoronoiFVM
using ExtendableGrids
using GridVisualize
using GLMakie

# Mutable struct to hold problem parameters
# This encapsulates all physical and numerical parameters
mutable struct ProblemData
    D::Float64      ## Diffusion coefficient
    v::Vector{Float64}  ## Velocity vector
end

# Bernoulli function used in the exponential fitting discretization
function bernoulli(x)
    if abs(x) < nextfloat(eps(typeof(x)))
        return 1
    end
    return x / (exp(x) - 1)
end

function exponential_flux!(f, u, edge, data)
    vh = project(edge, data.v)
    Bplus = data.D * bernoulli(vh / data.D)
    Bminus = data.D * bernoulli(-vh / data.D)
    f[1] = Bminus * u[1, 1] - Bplus * u[1, 2]
    return nothing
end

function central_flux!(f, u, edge, data)
    vh = project(edge, data.v)
    f[1] = data.D * (u[1, 1] - u[1, 2]) + vh * (u[1,1] + u[1,2]) / 2
    return nothing
end

function outflow!(f, u, node, data)
    if node.region == 2
        f[1] = data.v[1] * u[1]
    end
    return nothing
end

function main(; n = 10, D = 0.01, v = 1.0, tend = 100)

    # Create a one-dimensional discretization
    h = 1.0 / n
    grid = simplexgrid(0:h:1)

    # Initialize problem parameters in data structure
    problem_data = ProblemData(D, [v])

    sys = VoronoiFVM.System(
        grid,
        VoronoiFVM.Physics(;
            flux = central_flux!,
            breaction = outflow!,
            data = problem_data
        )
    )

    # Add species 1 to region 1
    enable_species!(sys, 1, [1])

    # Set boundary conditions
    boundary_neumann!(sys, 1, 1, 0.0)

    # Create a solution array
    inival = unknowns(sys)
    inival[1, :] .= map(x -> 1 - 2x, grid)

    # Transient solution of the problem
    control = VoronoiFVM.SolverControl()
    control.Δt = 0.01 * h
    control.Δt_min = 0.01 * h
    control.Δt_max = 0.1 * tend
    tsol = solve(sys; inival, times = [0, tend], control)

    vis = GridVisualizer(; Plotter = GLMakie, resolution=(600,400))

    # Use the record function to loop and save frames
    fig = reveal(vis)
    GLMakie.record(fig, "animation_.mp4", 1:length(tsol.t); framerate = 30) do i
        # Clear the previous plot to ensure a fresh frame
        # (for Makie, it's often better to create a new scene or clear axes)
        # For simplicity here, we'll just over-plot.
        
        scalarplot!(
            vis[1, 1], grid, tsol[1, :, i];
            flimits = (0, 1),
            colormap = :viridis,
            title = "t=$(round(tsol.t[i], digits=2))"
        )
    end
    return tsol
end

using Test
function runtests()
    tsol = main()
    @test maximum(tsol) <= 1.0 && maximum(tsol.u[end]) < 1.0e-20
    return nothing
end

runtests()