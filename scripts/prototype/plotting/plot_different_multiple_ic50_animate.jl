drug = "rg"

include("../../../model/$(drug)_pkpd2.jl")
include("../../../model/$(drug)_dosing2.jl")
include("../../../model/$(drug)_params.jl")

println("this is the drug: ", drug)
println("this is the IC50 for TMZ: ", IC50_1, " converted to amount: ", IC50_1 * 140)
# Determine IC50 values
IC50_value = IC50_2 * 1000 * Vpla
println("this the IC50 for $drug: ", IC50_2, " converted to amount: ",IC50_value)

# include("../../../scripts/setup/init_integrate.jl")

function dose_affect!(integrator)
    SciMLBase.set_proposed_dt!(integrator, 1)
    integrator.u[2] += integrator.u[end]
end
# cb = PeriodicCallback(dose_affect!, 24.0, initial_affect=true)

event_times = collect(1:18).*24.0
cb = PresetTimeCallback(event_times, dose_affect!)
tspan = (0.0, event_times[end]+(10.0*24.0))
    
p = ode_params
u0 = [17.7, 0.0, 0.0, 0.0, 5000.0]
prob = ODEProblem(simple_pkpd!, u0, tspan, p)
using ModelingToolkit
sys = modelingtoolkitize(prob)
sys = structural_simplify(sys)
func = ODEFunction(sys, jac=true)
prob_jac = ODEProblem(func, u0, tspan, p)
sol = solve(prob_jac, Rodas4P2(; linsolve = nothing), callback=cb)

# auxillary functions computing for plotting
cPlaRG_sol = sol[3,:] ./ (1000 * Vpla)
exp2_sol = (cPlaRG_sol ./(psi*IC50_2)).^gamma_2
E_sol = (exp2_sol.*Imax_2) ./ (exp2_sol .+ 1)
t_sol = (-log.(log.(sol[1,:]./K) ./ log(C0/K)) ./ r) .+ 72
fun_sol = K .* (C0/K).^exp.(-r .* t_sol)
delta_sol = (E_sol .* fun_sol) ./ (72 .* sol[1,:])

# Get max values for each y-axis
max_y1 = 40 # Replace with actual max value of data associated with primary y-axis
# Set the same relative position for IC50 on both y-axes
relative_position_ic50 = IC50_value / max_y1
# Apply the relative position to secondary y-axis
ic50_secondary_axis = 0.5/relative_position_ic50
max_y2 = ic50_secondary_axis
p = plot(sol.t, sol[3,:], label="PlaRG", xlabel="Time", ylabel="Value", legend=:topleft, ylim=[0.0, max_y1])
plot!(sol.t, ones(length(sol.t)).*(IC50_2 * 1000 * Vpla), label="IC50", ls=:dash, alpha=0.5, color="blue")
p2 = twinx(p)
plot!(p2, sol.t, E_sol, label="Effect", color="red", ylabel="Effect", ylim=[0.0, max_y2],ls=:dash, legend=:topright, alpha=0.5)
plot!(p2, sol.t, ones(length(sol.t)).*(0.5), label="IC50", color="red", ls=:dash, alpha=0.5)


mults = [0.01, 0.1, 0.5, 1.0, 1.5, 2.0]
avg_human_surface_area=1.7
# dosage guess for each multiplier of IC50
initial_guess = [0.0, 1.0, 2.0, 5.0, 10.0, 20.0, 100.0, 500.0, 1000.0, 2000.0, 3000.0].*avg_human_surface_area # Adjust this as needed
# initialize array to hold patient doses for each multiplier
gradation = [collect(0.01:0.01:0.09)...,collect(0.1:0.1:0.9)...,collect(0.91:0.01:0.99)...]
gradation = [0.000001, 0.00001, 0.0001, 0.001, collect(0.01:0.01:0.99)..., 0.999, 0.9999]
gradation = [0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 0.25, 0.5, 0.6, 0.75, 0.85, 0.99, 0.999, 0.9999]
#              11                                                               3000                56000
# pk_param_values= Cl2, ka2, Vpla, Q, Vtis
inputs_per_mult = zeros(length(gradation))
for (j,mult) in enumerate(gradation)
    function bc_drug_pk!(residual, sol, p, t)
        cPlaRG_sol = sol[3,:] ./ (1000 * Vpla)
        cPlaRG_sol = erelu.(cPlaRG_sol)
        exp2_sol = (cPlaRG_sol ./(psi*IC50_2)).^gamma_2
        E_sol = (exp2_sol.*Imax_2) ./ (exp2_sol .+ 1)
        residual[1] = maximum(E_sol) - mult # Ensure the plasma compartment matches ic50 amount conversion at end time
        residual[2] = sol[1][3] # Ensure the second state is 0 at the beginning
        residual[3] = sol[1][4] # Ensure the third state is 0 at the beginning
        residual[4] = sol[1][2]
        residual[5] = sol[1][1]
    end
    u0_temp = [17.7, 0.0, 0.0, 0.0, 1000.0]
    prob_temp = BVProblem(simple_pkpd!, bc_drug_pk!, u0_temp, tspan, ode_params)
    sol_temp = solve(prob_temp, Shooting(Rodas4P2()), callback=cb, alg_hints=[:stiff])
    inputs_per_mult[j] = sol_temp[1][5]
end
#plot bar plot of the doses for effect grdations
# p = bar(gradation, inputs_per_mult, label="Dose", xlabel="Effect", ylabel="Dose Amount (mg)", title="Dose vs. Effect", color="gray")
# i=1
# mask = zeros(length(inputs_per_mult))
# mask[i] = inputs_per_mult[i]
# bar!(gradation, mask, color="red")
# savefig(p, "results/$(drug)_dose_vs_effect.png")
# savefig(p, "results/$(drug)_dose_vs_effect.svg")

for (i, init_dose) in enumerate(inputs_per_mult)

    p = bar(gradation, inputs_per_mult, label="Dose", xlabel="Effect", ylabel="Dose Amount (mg)", title="Dose vs. Effect", color="gray")
    mask = zeros(length(inputs_per_mult))
    mask[i] = init_dose
    bar!(gradation, mask, color="red")

    prob_temp = remake(prob_jac, u0=[17.7, 0.0, 0.0, 0.0, init_dose])
    sol = solve(prob_temp, Rodas4P2(; linsolve = nothing), callback=cb)

    # auxillary functions computing for plotting
    cPlaRG_sol = sol[3,:] ./ (1000 * Vpla)
    exp2_sol = (cPlaRG_sol ./(psi*IC50_2)).^gamma_2
    E_sol = (exp2_sol.*Imax_2) ./ (exp2_sol .+ 1)
    t_sol = (-log.(log.(sol[1,:]./K) ./ log(C0/K)) ./ r) .+ 72
    fun_sol = K .* (C0/K).^exp.(-r .* t_sol)
    delta_sol = (E_sol .* fun_sol) ./ (72 .* sol[1,:])

    drug_auc = trapezoidal_rule(sol.t, sol[3,:])
    cell_auc = trapezoidal_rule(sol.t, sol[1,:])

    # Get max values for each y-axis
    max_y1 = 35#maximum(sol[3,:])*1.1 # Replace with actual max value of data associated with primary y-axis
    # Set the same relative position for IC50 on both y-axes
    relative_position_ic50 = IC50_value / max_y1
    # Apply the relative position to secondary y-axis
    ic50_secondary_axis = 0.5/relative_position_ic50
    max_y2 = ic50_secondary_axis

    p = plot(sol.t, sol[3,:], label="PlaRG", xlabel="Time", ylabel="Value", legend=:topleft, ylim=[0.0, max_y1], title="effect_$(gradation[i])_$(init_dose)")
    plot!(sol.t, ones(length(sol.t)).*(IC50_2 * 1000 * Vpla), label="IC50", ls=:dash, alpha=0.5, color="blue")
    annotate!(sol.t[1],max_y1*0.6, text("Drug AUC: $drug_auc", :left))
    annotate!(sol.t[1],max_y1*0.5, text("Cell AUC: $cell_auc", :left))
    p2 = twinx(p)
    plot!(p2, sol.t, E_sol, label="Effect", color="red", ylabel="Effect", ylim=[0.0, max_y2],ls=:dash, legend=:topright, alpha=0.5)
    plot!(p2, sol.t, ones(length(sol.t)).*(0.5), label="IC50", color="red", ls=:dash, alpha=0.5)
    savefig(p, "results/trial$(i)_$(drug)_effect_$(gradation[i])_$(init_dose).png")
end






using Plots
using Printf
gr()

# Set up the animation
anim = @animate for (i, init_dose) in enumerate(inputs_per_mult)
    # Bar plot
    gradation_labels = [string(x) for x in gradation]
    bar_plot = bar(gradation_labels, inputs_per_mult, label="Dose", xlabel="Effect", ylabel="Dose Amount (mg)", color="gray", ylim=[0.0, 10000])
    mask = zeros(length(inputs_per_mult))
    mask[i] = init_dose
    bar!(bar_plot, gradation_labels, mask, label="selected dose", color="red")
    formatted_init_dose = @sprintf("%.2f", init_dose)
    title!(bar_plot, "Effect: $(gradation[i]), Dose: $(formatted_init_dose)")

    # Solving the ODE problem (as per your existing code)
    prob_temp = remake(prob_jac, u0=[17.7, 0.0, 0.0, 0.0, init_dose])
    sol = solve(prob_temp, Rodas4P2(; linsolve = nothing), callback=cb)

    # auxillary functions computing for plotting
    cPlaRG_sol = sol[3,:] ./ (1000 * Vpla)
    exp2_sol = (cPlaRG_sol ./(psi*IC50_2)).^gamma_2
    E_sol = (exp2_sol.*Imax_2) ./ (exp2_sol .+ 1)
    t_sol = (-log.(log.(sol[1,:]./K) ./ log(C0/K)) ./ r) .+ 72
    fun_sol = K .* (C0/K).^exp.(-r .* t_sol)
    delta_sol = (E_sol .* fun_sol) ./ (72 .* sol[1,:])

    drug_auc = @sprintf("%.2f", trapezoidal_rule(sol.t, sol[3,:]))
    cell_auc = @sprintf("%.2f", trapezoidal_rule(sol.t, sol[1,:]))  

    # Get max values for each y-axis
    max_y1 = 35#maximum(sol[3,:])*1.1 # Replace with actual max value of data associated with primary y-axis
    # Set the same relative position for IC50 on both y-axes
    relative_position_ic50 = IC50_value / max_y1
    # Apply the relative position to secondary y-axis
    ic50_secondary_axis = 0.5/relative_position_ic50
    max_y2 = ic50_secondary_axis
    # Trajectory subplot
    times = sol.t./24.0
    trajectory_plot = plot(times, sol[3,:], label="PlaRG", xlabel="Time (days)", ylim=[0.0, max_y1], legend=:topleft)
    plot!(trajectory_plot, times, ones(length(sol.t)).*(IC50_2 * 1000 * Vpla), label="IC50", ls=:dash, alpha=0.5, color="blue")
    annotate!(trajectory_plot, times[1], max_y1 * 0.6, text("Drug AUC: $drug_auc", :left))
    annotate!(trajectory_plot, times[1], max_y1 * 0.5, text("Cell AUC: $cell_auc", :left))
    p2 = twinx(trajectory_plot)
    plot!(p2, times, E_sol, label="Effect", color="red", ylabel="Effect", ylim=[0.0, max_y2], ls=:dash, legend=:topright, alpha=0.5)
    plot!(p2, times, ones(length(sol.t)).*(0.5), label="IC50", color="red", ls=:dash, alpha=0.5)

    # Combine subplots into a single figure
    combined_plot = plot(bar_plot, trajectory_plot, layout = (1, 2), size=(1440, 720))
    
    # Optional: save individual frames if needed
    # savefig(combined_plot, "frame_$i.png")
end

# Save the animation
gif(anim, "RG_dose_effect_outputs_animation.gif", fps = 10)
# video(anim, "dose_trajectory_animation.mp4", fps = 10)


using Printf

anim = @animate for (i, init_dose) in enumerate(inputs_per_mult)
    formatted_init_dose = @sprintf("%.2f", init_dose)

    # Bar plot with thicker bars
    bar_width = 0.7  # Adjust this value as needed
    bar_plot = bar(gradation, inputs_per_mult, label="Dose", xlabel="Effect", ylabel="Dose Amount (mg)", color="gray", bar_width=bar_width, yaxis=:log)
    
    # Add arrow to indicate selected dose
    arrow_x = gradation[i]
    arrow_y = inputs_per_mult[i]
    annotate!(bar_plot, [(arrow_x, arrow_y, text("Selected", 8, :center, :above)), (arrow_x, arrow_y, arrow(0.5, 30, :red))])

    title!(bar_plot, "Effect: $(gradation[i]), Dose: $(formatted_init_dose)")

    # Solving the ODE problem (as per your existing code)
    prob_temp = remake(prob_jac, u0=[17.7, 0.0, 0.0, 0.0, init_dose])
    sol = solve(prob_temp, Rodas4P2(; linsolve = nothing), callback=cb)

    # auxillary functions computing for plotting
    cPlaRG_sol = sol[3,:] ./ (1000 * Vpla)
    exp2_sol = (cPlaRG_sol ./(psi*IC50_2)).^gamma_2
    E_sol = (exp2_sol.*Imax_2) ./ (exp2_sol .+ 1)
    t_sol = (-log.(log.(sol[1,:]./K) ./ log(C0/K)) ./ r) .+ 72
    fun_sol = K .* (C0/K).^exp.(-r .* t_sol)
    delta_sol = (E_sol .* fun_sol) ./ (72 .* sol[1,:])

    drug_auc = @sprintf("%.2f", trapezoidal_rule(sol.t, sol[3,:]))
    cell_auc = @sprintf("%.2f", trapezoidal_rule(sol.t, sol[1,:]))  

    # Get max values for each y-axis
    max_y1 = 35#maximum(sol[3,:])*1.1 # Replace with actual max value of data associated with primary y-axis
    # Set the same relative position for IC50 on both y-axes
    relative_position_ic50 = IC50_value / max_y1
    # Apply the relative position to secondary y-axis
    ic50_secondary_axis = 0.5/relative_position_ic50
    max_y2 = ic50_secondary_axis
    # Trajectory subplot
    trajectory_plot = plot(sol.t, sol[3,:], label="PlaRG", xlabel="Time", ylim=[0.0, max_y1], legend=:topleft)
    plot!(trajectory_plot, sol.t, ones(length(sol.t)).*(IC50_2 * 1000 * Vpla), label="IC50", ls=:dash, alpha=0.5, color="blue")
    annotate!(trajectory_plot, sol.t[1], max_y1 * 0.6, text("Drug AUC: $drug_auc", :left))
    annotate!(trajectory_plot, sol.t[1], max_y1 * 0.5, text("Cell AUC: $cell_auc", :left))
    p2 = twinx(trajectory_plot)
    plot!(p2, sol.t, E_sol, label="Effect", color="red", ylabel="Effect", ylim=[0.0, max_y2], ls=:dash, legend=:topright, alpha=0.5)
    plot!(p2, sol.t, ones(length(sol.t)).*(0.5), label="IC50", color="red", ls=:dash, alpha=0.5)

    # Combine subplots into a single figure
    combined_plot = plot(bar_plot, trajectory_plot, layout = (1, 2), size=(1080, 720))
    
    # Optional: save individual frames if needed
    # savefig(combined_plot, "frame_$i.png")
end


gif(anim, "dose_trajectory_animation3.gif", fps = 10)




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






