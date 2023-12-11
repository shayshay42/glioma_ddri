drug = "rg"

include("../../../model/$(drug)_pkpd2.jl")
include("../../../model/$(drug)_dosing2.jl")
include("../../../model/$(drug)_params.jl")

println("this is the drug: ", drug)
println("this is the IC50 for TMZ: ", IC50_1, " converted to amount: ", IC50_1 * 140)
println("this the IC50 for $drug: ", IC50_2, " converted to amount: ", IC50_2 * 1000 * Vpla)

# include("../../../scripts/setup/init_integrate.jl")

function dose_affect!(integrator)
    integrator.u[2] += integrator.u[end]
end
# cb = PeriodicCallback(dose_affect!, 24.0, initial_affect=true)

event_times = collect(1:18).*24.0
cb = PresetTimeCallback(event_times, dose_affect!)
tspan = (0.0, event_times[end]+(10.0*24.0))
    
p = ode_params
u0 = [17.7, 0.0, 0.0, 0.0, 5000.0]
prob = ODEProblem(simple_pkpd!, u0, tspan, p)
sol = solve(prob, callback=cb)

mults = [0.01, 0.1, 0.5, 1.0, 1.5, 2.0]
avg_human_surface_area=1.7
# dosage guess for each multiplier of IC50
initial_guess = [20.0, 100.0, 500.0, 1000.0, 2000.0, 3000.0].*avg_human_surface_area # Adjust this as needed
# initialize array to hold patient doses for each multiplier

# pk_param_values= Cl2, ka2, Vpla, Q, Vtis
inputs_per_mult = zeros(length(mults))
for (j,mult) in enumerate(mults)
    function bc_drug_pk!(residual, sol, p, t)
        residual[1] = maximum(sol[2,:]) - ((IC50_2 * (1000 * Vpla))*mult) # Ensure the second state matches ic50 amount conversion at end time
        residual[2] = sol[1][2] # Ensure the seocnd state is 0 at the beginning
        residual[3] = sol[1][3] # Ensure the third state is 0 at the beginning
        residual[4] = sol[1][1]
    end
    u0_temp = [0.0, 0.0, 0.0, initial_guess[j]]
    p_temp = pk_param_values
    prob_temp = BVProblem(drug_pk_dose!, bc_drug_pk!, u0_temp, tspan, p_temp)
    sol_temp = solve(prob_temp, Shooting(KenCarp4()), callback=cb)
    input_dose = sol_temp[4, 1]
end

# Define the saving function
function save_func(u, t, integrator)
    AbsRG, PlaRG, TisRG, dose = u

    cPlaRG = PlaRG / (1000 * Vpla)
    # cPlaRG = erelu(cPlaRG)
    # xi = IC50_1/IC50_2
    # pi1 = psi * IC50_1
    # exp1 = (cCSFTMZ/(psi*IC50_1))^gamma_1
    # exp2 = (cPlaRG /(psi*IC50_2))^gamma_2
    # exp12 = exp1 * exp2

    # Imax_sum = Imax_1 + Imax_2
    # E_num = Imax_1 * exp1 + Imax_2 * exp2 + Imax_sum * exp12 - Imax_1 * Imax_2 * exp12
    # E_num = Imax_2 * exp2
    # E_den = exp1 + exp2 + exp12 + 1
    # E_den = exp2 + 1
    # E = E_num / E_den
    exp2 = (cPlaRG /(psi*IC50_2))^gamma_2
    E = Imax_2 * exp2 / (exp2 + 1)

    # log_C0K = log(C0/K)
    # t1 = -log(log(C/K) / log_C0K) / r
    # t2 = t1 + 72
    # fun = K * exp(log_C0K * exp(-r * t2))
    # delta = (E * fun) / (72 * C)

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






