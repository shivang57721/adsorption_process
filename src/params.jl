Base.@kwdef struct PhysicalConstants
    R::Float64 = 8.314         # Ideal gas constant [J mol⁻¹ K⁻¹]

    γ₁::Float64 = 0.7            # Axial dispersion coefficient correlation parameter [dimensionless]
    γ₂::Float64 = 0.5            # Axial dispersion coefficient correlation parameter [dimensionless]
    μ::Float64  = 1.789e-5      # Dynamic viscosity of the gas mixture [Pa s] or [kg m⁻¹ s⁻¹]

    Cₚ_CO2::Float64 = 37.520      # Molar heat capacity of CO₂ gas [J mol⁻¹ K⁻¹]
    Cₚ_H2O::Float64 = 34.2       # Molar heat capacity of H₂O gas [J mol⁻¹ K⁻¹]
    Cₚ_N2::Float64  = 29.171       # Molar heat capacity of N₂ gas [J mol⁻¹ K⁻¹]
    Cₚ_wall::Float64 = 4.0e6     # Volumetric heat capacity of the column wall (ρ_wall * C_p_wall) [J m⁻³ K⁻¹]

    MW_CO2::Float64 = 44         # Molar weight of CO₂ [g mol⁻¹] or [kg kmol⁻¹]
    MW_N2::Float64  = 28         # Molar weight of N₂ [g mol⁻¹] or [kg kmol⁻¹]
    MW_H2O::Float64 = 18         # Molar weight of H₂O [g mol⁻¹] or [kg kmol⁻¹]
end

Base.@kwdef struct ColumnParams # Default params are Young
    Rᵢ::Float64 = 0.04           # Inner radius of the column [m]
    Rₒ::Float64 = 0.041          # Outer radius of the column [m]
    L::Float64  = 0.01           # Column length [m]

    h_L::Float64 = 14            # Heat transfer coefficient from gas to column wall [W m⁻² K⁻¹]
    h_W::Float64 = 22000         # Heat transfer coefficient from column wall to environment [W m⁻² K⁻¹]
end

Base.@kwdef struct SorbentParams # Default params are Lewatit
    ε_bed::Float64   = 0.4       # Bed voidage (interparticle porosity) [dimensionless]
    ε_total::Float64 = 0.54      # Total porosity (interparticle + intraparticle) [dimensionless]
    dₚ::Float64      = 0.00052   # Particle diameter [m]
    ρ_bed::Float64   = 528       # Bed density [kg m⁻³]
    ρ_particle::Float64 = 880    # Particle density [kg m⁻³]
    k_CO2::Float64 = 0.003       # Lumped mass transfer coefficient for CO₂ [s⁻¹]
    k_H2O::Float64 = 0.0086      # Lumped mass transfer coefficient for H₂O [s⁻¹]
    ΔH_CO2::Float64 = -70000     # Heat of adsorption for CO₂ [J mol⁻¹]
    ΔH_H2O::Float64 = -46000     # Heat of adsorption for H₂O [J mol⁻¹]
    Dₘ::Float64 = 1.3e-5       # Molecular diffusivity of the gas mixture [m² s⁻¹]
    C_solid::Float64 = 1580     # Heat capacity of the solid sorbent [J kg⁻¹ K⁻¹]

    q_star_CO2::Function = Toth_isotherm_CO2_modified_H2O
    q_star_H2O::Function = GAB_isotherm_H2O_Tfunction_Resins

    pho_bulk_packing::Float64 =  528
end

# A struct to hold all operating and derived parameters
@enum StepType Adsorption Heating Desorption Cooling Pressurization PressurizationReset
mutable struct OperatingParameters
    # Independent parameters (that you change between steps)
    u_feed::Float64
    T_feed::Float64
    T_amb::Float64
    y_CO2_feed::Float64
    y_H2O_feed::Float64
    duration::Float64
    step_name::StepType
    ΔT::Float64

    P_out::Float64
    P_out_func::Function
    
    # Dependent parameters (calculated from the above)
    y_N2_feed::Float64
    c_total_feed::Float64
    c_N2_feed::Float64
    c_CO2_feed::Float64
    c_H2O_feed::Float64

    D_L::Float64
    C_gas_feed::Float64
    K_L::Float64
    
    # Constructor to ensure derived parameters are always consistent
    function OperatingParameters(; u_feed, T_feed, T_amb, y_CO2_feed, y_H2O_feed, P_out, duration, step_name, ΔT=NaN)
        # Create a new instance
        p = new()

        # Assign independent parameters
        p.u_feed = u_feed
        p.T_feed = T_feed
        p.T_amb = T_amb
        p.P_out = P_out
        p.y_CO2_feed = y_CO2_feed
        p.y_H2O_feed = y_H2O_feed
        p.duration = duration
        p.step_name = step_name
        p.ΔT = ΔT

        return p
    end
end

# A helper function to recalculate dependent parameters
function update_derived_params!(p::OperatingParameters, phys_consts::PhysicalConstants, sorb_params::SorbentParams)
    p.y_N2_feed = 1.0 - p.y_CO2_feed - p.y_H2O_feed
    p.c_total_feed = p.P_out / (phys_consts.R * p.T_feed)
    p.c_N2_feed  = p.y_N2_feed * p.c_total_feed
    p.c_CO2_feed = p.y_CO2_feed * p.c_total_feed
    p.c_H2O_feed = p.y_H2O_feed * p.c_total_feed
    
    p.D_L = phys_consts.γ₁ * sorb_params.Dₘ + phys_consts.γ₂ * sorb_params.dₚ * p.u_feed / sorb_params.ε_bed
    p.C_gas_feed = p.c_CO2_feed * phys_consts.Cₚ_CO2 + p.c_H2O_feed * phys_consts.Cₚ_H2O + p.c_N2_feed * phys_consts.Cₚ_N2
    p.K_L = p.D_L * p.C_gas_feed

    return nothing
end

function copy_params!(dest::OperatingParameters, src::OperatingParameters)
    # This copies the values of each field from src to dest
    for name in fieldnames(OperatingParameters)
        setfield!(dest, name, getfield(src, name))
    end
    return nothing
end

Psat_H2O(T) = 611.21 * exp((18.678 - T / 234.5) * T / (T + 273.15 - 16.01))

function GAB_isotherm_H2O_Tfunction_Resins(T, p_H2O)
    # --- GAB H2O isotherm params ---
    qm::Float64 = 3.63            # Monolayer water capacity [mol kg⁻¹]
    C::Float64 = 47110            # GAB energy parameter [J mol⁻¹]
    D::Float64 = 0.023744          # GAB temperature dependence parameter [K⁻¹]
    F::Float64 = 57706            # GAB energy parameter [J mol⁻¹]
    G::Float64 = -47.814          # GAB energy parameter [J mol⁻¹ K⁻¹]

    E1 = C - exp(D * T)
    E2_9 = F + G * T
    E10 = -44.38 * T + 57220
    c = exp((E1-E10)/(R*T))
    k = exp((E2_9 - E10)/(R*T))

    x = p_H2O / Psat_H2O(T - 273)

    q_star = qm * k * c * x / ((1 - k * x) * (1 + (c - 1) * k * x))

    q_star
end

function Toth_isotherm_CO2_modified_H2O_Stampi(T, R, p_CO2, q_H2O, params)
    ns0 = params.ns0
    chi = params.chi
    T_ref = params.T_ref
    DH = params.DH
    b0 = params.b0
    t0 = params.t0
    alfa = params.alfa
    gamma = params.gamma
    beta = params.beta

    ns = ns0 * exp(chi*(1 - T / T_ref)) * (1 / (1 - gamma * q_H2O))
    b = b0 * exp(-DH / R / T_ref * (T_ref / T - 1)) * (1 + beta * q_H2O)
    t = t0 + alfa * (1 - T_ref / T)

    p_CO2_safe = max(eps(eltype(p_CO2)), p_CO2)

    ns * b * p_CO2 * 1e-5 / (1 + (b * p_CO2_safe * 1e-5) ^ t) ^ (1/t)
end

function Toth_isotherm_CO2_modified_H2O(T, R, p_CO2, q_H2O, params)
    ns0 = params.ns0
    chi = params.chi
    T_ref = params.T_ref
    DH = params.DH
    b0 = params.b0
    t0 = params.t0
    alfa = params.alfa
    gamma = params.gamma
    beta = params.beta

    ns = ns0 * exp(chi*(1 - T / T_ref)) * (1 / (1 - gamma * q_H2O))
    b = b0 * exp(-DH / R / T) * (1 + beta * q_H2O)
    t = t0 + alfa * (1 - T_ref / T)

    p_CO2_safe = max(eps(eltype(p_CO2)), p_CO2)
    
    if t ≥ 0
        if (b * p_CO2_safe * 1e-5) < 0
            @show b.value
            @show q_H2O.value
            @show p_CO2_safe.value
        end
        q_star = ns * b * p_CO2 * 1e-5 / (1 + (b * p_CO2_safe * 1e-5) ^ t) ^ (1/t)
    else
        q_star = ns * (1 + (b * p_CO2_safe * 1e-5) ^ (-t)) ^ (-1/t)
    end

    return q_star
end

function GAB_isotherm_H2O(T, p_H2O, params)
    Cm = params.Cm
    Cg = params.Cg
    K_ads = params.K_ads

    x = p_H2O / Psat_H2O(T-273)

    q_H2O = Cm * Cg * K_ads * x /((1 - K_ads * x) * (1 + (Cg - 1) * K_ads * x))

    return q_H2O
end