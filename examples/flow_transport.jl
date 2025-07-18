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
    iv = 3
end

