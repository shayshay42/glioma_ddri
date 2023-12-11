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

include("../../../scripts/setup/init_integrate.jl")

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
saveat = 0.1
saving_cb = SavingCallback(save_func, SavedValues(Float64, Tuple{Float64, Float64, Float64, Float64, Float64}), saveat=saveat)

# Add the callback to the problem
cb = CallbackSet(saving_cb, hit)

#can i affect the dosing callback amoutn from avg_human_surface_area
tmz_treat_dose = 0.0#75.0*avg_human_surface_area
tmz_adjuv_dose = 0.0#150.0*avg_human_surface_area
# Solve the problem with the callback
edit_prob = remake(prob_jac, tspan=(0,600))
sol = solve(edit_prob,callback=cb, saveat=saveat)

# Access the saved values
saved_values = saving_cb.affect!.saved_values.saveval
timepoints = saving_cb.affect!.saved_values.t

using Plots

# Extract E and delta values
E_values = [val[1] for val in saved_values] # Extracting first element of each tuple
delta_values = [val[2] for val in saved_values] # Extracting second element of each tuple
plarrg_values = [val[3] for val in saved_values] # Extracting second element of each tuple
cplarrg_values = [val[4] for val in saved_values] # Extracting second element of each tuple
ccsftmz_values = [val[5] for val in saved_values] # Extracting second element of each tuple

# Plot E and delta against time

plot(sol.t, sol[7,:]./(1000*Vpla), label="[Plasma RG]", xlabel="Time", ylabel="Value")
plot!(sol.t, ones(length(sol.t)).*(IC50_2), label="[IC50]")

p1 = plot(sol.t, sol[7,:], label="Plasma RG", xlabel="Time", ylabel="Value", legend=:topleft)
plot!(sol.t, ones(length(sol.t)).*(IC50_2 * 1000 * Vpla), label="IC50")
p2 = twinx(p1) 
plot!(p2, timepoints, E_values, label="Effect", color="red", ylabel="Effect", ylim=[0.0,1.0],ls=:dash, legend=:topright, alpha=0.5)

plot(sol.t, E_values, label="E", xlabel="Time", ylabel="Value", title="E over Time")
plot!(E_values, label="E_only")




plot(timepoints, delta_values, label="Delta", xlabel="Time", ylabel="Value", title="Delta over Time")
plot(timepoints, plarrg_values, label="PlaRG", xlabel="Time", ylabel="Value", title="PlaRG over Time")
plot(timepoints, cplarrg_values, label="cPlaRG", xlabel="Time", ylabel="Value", title="cPlaRG over Time")
plot(timepoints, ccsftmz_values, label="cCSFTMZ", xlabel="Time", ylabel="Value", title="cCSFTMZ over Time")