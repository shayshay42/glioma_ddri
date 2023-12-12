# want to reverse the ic50 parameter for each patient back to the amount in the absroption compartment
# or instead try many different amounts and see where it ends up with realtion to ic50
# plot the effect amount as the constant dopes as a fraction of the ic50 is given on the schedule

#change the order the equations are computed int he pkpd

drug = "rg"

include("../../../model/$(drug)_pkpd2.jl")
include("../../../model/$(drug)_dosing2.jl")
include("../../../model/$(drug)_params.jl")

println("this is the drug: ", drug)
println("this is the IC50 for TMZ: ", IC50_1, " converted to amount: ", IC50_1 * 140)
println("this the IC50 for $drug: ", IC50_2, " converted to amount: ", IC50_2 * 1000 * Vpla)

using DifferentialEquations, ModelingToolkit

@register dose_heaviside()

u0 = zeros(9)
u0[1] = 17.7
tspan = (0,end_time+7).*hours
p = [ode_params; doses]
prob = ODEProblem(pk_pd_alt!,u0,tspan,p)
# sys = modelingtoolkitize(prob)
# sys = structural_simplify(sys)
# func = ODEFunction(sys, jac=true)
# prob_jac = ODEProblem(func, u0, tspan, p)
sol = solve(prob, tstops=inject_times, d_discontinuities=inject_times)







# Define the saving function
function save_func(u, t, integrator)
    gamma_1, psi, C0, D0, r, K, BW, IC50_1, Imax_1, IC50_2, gamma_2, Imax_2, xi, VD1, Cl1, k23, ka1, k32, Cl2, ka2, Vpla, Q, Vtis =  p[1:length(ode_params)]
    C, D, AbsTMZ, PlaTMZ, CSFTMZ, AbsRG, PlaRG, TisRG, cAUC = u

    # Calculations for E and delta as in pk_pd! function
    cCSFTMZ = CSFTMZ / 140
    cPlaRG = PlaRG / (1000 * Vpla)

    #domain error
    # cCSFTMZ = erelu(cCSFTMZ)
    # cPlaRG = erelu(cPlaRG)
    # C = erelu(C)

    exp1 = (cCSFTMZ /(psi*IC50_1))^gamma_1
    exp2 = (cPlaRG  /(psi*IC50_2))^gamma_2
    exp12 = exp1 * exp2

    Imax_sum = Imax_1 + Imax_2
    E_num = Imax_1 * exp1 + Imax_2 * exp2 + Imax_sum * exp12 - Imax_1 * Imax_2 * exp12
    E_den = exp1 + exp2 + exp12 + 1
    E = E_num / E_den

    t = (-log(log(C/K) / log(C0/K)) / r) + 72
    fun = K * (C0/K) ^ exp(-r * t)
    delta = (E * fun) / (72 * C)

    # Return the values to save
    return (E, delta, PlaRG, cPlaRG, cCSFTMZ)
end

# Create a SavingCallback
saveat = 0.1
saving_cb = SavingCallback(save_func, SavedValues(Float64, Tuple{Float64, Float64, Float64, Float64, Float64}), saveat=saveat)

# Add the callback to the problem
cb = CallbackSet(saving_cb, hit)

#can i affect the dosing callback amoutn from avg_human_surface_area
tmz_treat_dose = 0.0#75.0*avg_human_surface_area
tmz_adjuv_dose = 0.0#150.0*avg_human_surface_area
# Solve the problem with the callback
edit_prob = remake(prob_jac, tspan=(0,600))
sol = solve(edit_prob, callback=hit)#, isoutofdomain=(y,p,t)->any(x->x<0,y))
