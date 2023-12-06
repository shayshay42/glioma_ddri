drug = "rg"

include("../../model/$(drug)_params.jl")

println("this is the drug: ", drug)
println("this is the IC50 for TMZ: ", IC50_1)
println("this the IC50 for $drug: ", IC50_2)

using DifferentialEquations, ModelingToolkit

function rg_pk(du, u, p, t)
    Cl2, ka2, Vpla, Q, Vtis = p
    AbsRG, PlaRG, TisRG = u
    dAbsRG = -ka2 * AbsRG
    dPlaRG = ka2 * AbsRG - (Cl2 / Vpla) * PlaRG + (Q / Vtis) * TisRG - (Q / Vpla) * PlaRG
    dTisRG = -(Q / Vtis) * TisRG + (Q / Vpla) * PlaRG
    du .= [dAbsRG, dPlaRG, dTisRG]
end
p = [Cl2, ka2, Vpla, Q, Vtis]
u0 = [0.0, 0.0, 0.0]
tspan = (0.0, 72.0)
prob = ODEProblem(drug_pk, u0, tspan, p)
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)
# Function to check if the field name matches 'PlaRG'
property_name_symbol = Symbol("Pla" * uppercase(drug))
drug_plasma_index = findfirst(isequal(property_name), fieldnames(typeof(ics)))
u0[drug_plasma_index] = IC50_2*(Vpla*1000)

tspan = (10,1).*hours
p = [ode_params; 1.0; 1.0;doses]
prob = ODEProblem(pk_pd!,u0,tspan,p)
sys = modelingtoolkitize(prob)
sys = structural_simplify(sys)
prob_jac = ODEProblem(sys, u0, tspan, p)

# Define the saving function
function save_func(u, t, integrator)
    gamma_1, psi, C0, D0, r, K, BW, IC50_1, Imax_1, IC50_2, gamma_2, Imax_2, xi, VD1, Cl1, k23, ka1, k32, Cl2, ka2, Vpla, Q, Vtis =  p[1:length(ode_params)]
    C, D, AbsTMZ, PlaTMZ, CSFTMZ, AbsRG, PlaRG, TisRG, cAUC = u

    # Calculations for E and delta as in pk_pd! function
    cCSFTMZ = CSFTMZ / 140
    cPlaRG = PlaRG / (1000 * Vpla)

    #domain error
    cCSFTMZ = erelu(cCSFTMZ)
    cPlaRG = erelu(cPlaRG)
    C = erelu(C)

    pi1 = psi * IC50_1
    exp1 = (cCSFTMZ/pi1)^gamma_1
    exp2 = ((xi * cPlaRG)/pi1)^gamma_2
    exp12 = exp1 * exp2

    Imax_sum = Imax_1 + Imax_2
    E_num = Imax_1 * exp1 + Imax_2 * exp2 + Imax_sum * exp12 - Imax_1 * Imax_2 * exp12
    E_den = exp1 + exp2 + exp12 + 1
    E = E_num / E_den

    log_C0K = log(C0/K)
    t1 = -log(log(C/K) / log_C0K) / r
    t2 = t1 + 72
    fun = K * exp(log_C0K * exp(-r * t2))
    delta = (E * fun) / (72 * C)

    # Return the values to save
    return (E, delta, PlaRG, cPlaRG, cCSFTMZ)
end

# Create a SavingCallback
saveat = 1.0
saving_cb = SavingCallback(save_func, SavedValues(Float64, Tuple{Float64, Float64, Float64, Float64, Float64}), saveat=saveat)

# Add the callback to the problem
cb = CallbackSet(saving_cb)

# Solve the problem with the callback
sol = solve(prob_jac, callback=cb)

# Access the saved values
saved_values = saving_cb.affect!.saved_values.saveval
timepoints = saving_cb.affect!.saved_values.t
