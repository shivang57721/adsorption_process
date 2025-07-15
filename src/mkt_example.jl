using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

@mtkmodel FOL begin
    @parameters begin
        τ = 3.0 # parameters
    end
    @variables begin
        x(t) = 0.0 # dependent variables
    end
    @equations begin
        D(x) ~ (1 - x) / τ
    end
end

using OrdinaryDiffEq
@mtkcompile fol = FOL()
prob = ODEProblem(fol, [], (0.0, 10.0), [])
sol = solve(prob)

using Plots
plot(sol)