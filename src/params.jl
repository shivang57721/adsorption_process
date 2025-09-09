Base.@kwdef struct PhysicalParams
    R::Float64 = 8.314

    γ₁::Float64 = 0.7
    γ₂::Float64 = 0.5
    Dₘ::Float64 = 4.3e-6 # Molecular diffusion
    μ::Float64  = 1.789e-5 # dynamic viscosity 

    C_solid::Float64 = 2070
    Cₚ_CO2::Float64 = 42.46
    Cₚ_H2O::Float64 = 73.1
    Cₚ_N2::Float64  = 29.1
    Cₚ_wall::Float64 = 4.0e6

    MW_CO2::Float64 = 44
    MW_N2::Float64  = 28
    MW_H2O::Float64 = 18
end

Base.@kwdef struct ColumnParams # Default params are Arvind
    Rᵢ::Float64 = 0.04          # Inner radius
    Rₒ::Float64 = 0.041         # Outer radius
    L::Float64  = 0.01          # Column length

    h_L::Float64 = 3        # Heat transfer coefficient from inside the column to the column wall
    h_W::Float64 = 26       # Heat transfer coefficient from column wall to environment
end

Base.@kwdef struct SorbentParams # Default params are Lewatit
    ε_bed::Float64   = 0.4
    ε_total::Float64 = 0.54
    dₚ::Float64      = 0.00052 # Particle size
    ρ_bed::Float64   = 528
    ρ_particle::Float64 = 880 # Particle density
    k_CO2::Float64 = 0.003
    k_H2O::Float64 = 0.0086
    ΔH_CO2::Float64 = -70000
    ΔH_H2O::Float64 = -46000

    # Improve these names later

    # CO2 Isotherm params
    qinf0_dry_CO2::Float64 = 4.86  # mol/kg
    chi_dry_CO2::Float64 = 0.0
    T0_dry_CO2::Float64 = 298.15
    DH_dry_CO2::Float64 = -117798.0 # J/mol
    b0_dry_CO2::Float64 = 2.85e-16 # 1/bar
    t0_dry_CO2::Float64 = 0.209  # 0.422
    alfa_dry_CO2::Float64 =  0.523 # 0.949

    qinf0_wet_CO2::Float64 = 9.035  # mol/kg
    chi_wet_CO2::Float64 = 0.0
    T0_wet_CO2::Float64 = 298.15
    DH_wet_CO2::Float64 = -203687.0 # J/mol
    b0_wet_CO2::Float64 = 1.23e-13 # 1/bar
    t0_wet_CO2::Float64 = 0.053  # 0.422
    alfa_wet_CO2::Float64 =  0.053 # 0.949
    A::Float64 = 1.532  # mol/kg

    gamma::Float64 = -0.137
    beta::Float64 = 5.612

    # GAB H2O isotherm params
    qm::Float64 = 3.63  # mol/kg
    C::Float64 = 47110  # J/mol
    D::Float64 = 0.023744  # 1/K
    F::Float64 = 57706  # J/mol
    G::Float64 = -47.814  # J/molK
end

Base.@kwdef struct OperatingParams
    u_feed::Float64 = 0.1
    T_amb::Float64  = 298.15
    T_feed::Float64 = 298.15
    P_out::Float64  = 1e5

    y_CO2_feed::Float64 = 0.0004
    y_H2O_feed::Float64 = 0.0115
    y_N2_feed::Float64  = 1 - 0.0004 - 0.0115
end

function q_star_H2O(qm, C, D, F, G, T, p_H2O)
    Psat_H2O(T) = 611.21 * exp((18.678 - T / 234.5) * T / (T + 273.15 - 16.01))

    E1 = C - exp(D * T)
    E2_9 = F + G * T
    E10 = -44.38 * T + 57220
    c = exp((E1-E10)/(R*T))
    k = exp((E2_9 - E10)/(R*T))

    x = p_H2O / Psat_H2O(T - 273)

    q_star = qm * k * c * x / ((1 - k * x) * (1 + (c - 1) * k * x))

    q_star
end

function q_star_CO2(ns0, chi, T_ref, b0, DH, t0, alfa, T, R, p_CO2, gamma, beta, q_H2O)
    ns = ns0 * exp(chi*(1 - T / T_ref)) * (1 / (1 - gamma * q_H2O))
    b = b0 * exp(-DH / R / T) * (1 + beta * q_H2O)
    t = t0 + alfa * (1 - T_ref / T)

    p_CO2_safe = max(eps(eltype(p_CO2)), p_CO2)

    ns * b * p_CO2 * 1e-5 / (1 + (b * p_CO2_safe * 1e-5)^t)^(1/t)
end