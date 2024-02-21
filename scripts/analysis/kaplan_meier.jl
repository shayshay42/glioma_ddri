using Plots, DifferentialEquations, MAT, Statistics, Serialization, ModelingToolkit
gr()

include("../../utilities/utils.jl")

drug = "gdc"
include("../../models/$(drug)_pkpd2.jl")
include("../../models/$(drug)_dosing2.jl")
include("../../models/$(drug)_params.jl")

drug_dosetimes = eval(Symbol("$(drug)_dosetimes"))

if drug == "rg"
   drug_params = ["Cl2", "ka2", "Vpla", "Q", "Vtis"]
else
   drug_params = ["ka2", "V2", "kel", "k12", "k21"]
end
TMZ_params = ["Cl1", "k23", "ka1"]

population = deserialize("./assets/$(drug)_vp5.jls")
optimum = deserialize("./results/optimum/optima_$(drug)_lfbgs_auc_maxinit_2023-11-17_dms.jls")
opt = Array{Float64, 2}(optimum[1:end-2,:])
loss_bit = optimum[end,:]

function death_time(v_init)
   v_init = max(v_init, 0.0)
   log((log(v_init/K)/log(113/K)))/r
end

# find the volume at a time where drug doages are zero in the blood
# use a callback to stop the simulation when the amount of u4 adn 7 are both below 1e-9
function find_volume_at_zero_dose(u, t, integrator)
   tmz_plasma = integrator.u[4]
   drug_plasma = integrator.u[7]
   (tmz_plasma < 1e-8 && drug_plasma < 1e-8) && integrator.t > (end_time+7.0)*hours#drug_dosetimes[end]#*hours
       # terminate!(integrator)
   # end
end
function affect!(integrator)
   integrator.terminate!
end

zero_dose_ = ContinuousCallback(find_volume_at_zero_dose, affect!)
cbset = CallbackSet(hit, zero_dose_)

u0 = zeros(Float64, 9)
u0[1] = 17.7
tspan = (0.0, end_time + 7.0) .* hours
ode_params = [ode_params; [0.1,1.0]]
p = [ode_params; doses]

prob = ODEProblem(pk_pd!, u0, tspan, p)
sys = modelingtoolkitize(prob)
sys = structural_simplify(sys)
prob_jac = ODEProblem(sys, u0, tspan, p, jac=true)
sol = solve(prob_jac, callback=cbset, saveat=0.5)


function get_dose(specification, patient_idx)
   if specification == "Min"
       return fill(min_drug_dosage/90, 90)  # Assuming min_dose is defined
   elseif specification == "Quarter"
       return doses ./ 4  # Assuming max_dose is defined
   elseif specification == "Half"
       return doses ./ 2
   elseif specification == "Max"
       return doses
   elseif specification == "Random"
       return rand(min_drug_dosage/90:doses[1], 90)
   elseif specification == "Optim"
       return opt[:,patient_idx]
   end
end

using Survival
# Initialize an empty plot
plt = plot(xlabel="Time (Days)", ylabel="Survival Probability", title="Kaplan-Meier Curve for $(uppercase(drug)) Dosing", grid=false)
nb_patients = size(population, 2)
for dose_spec in ["Optim", "Min", "Quarter", "Half", "Max", "Random"]
   println("\rDose Spec: $dose_spec")
   death_times = zeros(nb_patients)
   Threads.@threads for i in 1:nb_patients
       println("\rPatient: $i")
       dosing = get_dose(dose_spec, i)
       ode_p = [population[:, i];0.1;1.0]
       p = [ode_p;dosing]
       p_prob = remake(prob_jac, p=p)
       p_sol = solve(p_prob, callback=cbset, saveat=1)#, abstol=1e-10, reltol=1e-10,dtmax=1)#, alg_hints=[:stiff])
       sols = Array(p_sol)
       # println("Patient $i terminated at $(sols[1,end])")
       # println("Patient $i stopped at time $(p_sol.t[end]) vs $((end_time+7)*24) (hours)")
       #assuming it terminated when the drug was zero_dose_
       death_times[i] = death_time(sols[1,end]) + (end_time+7)*24
   end

   # Create an event indicator vector, all ones (1) since all are events
# Create an event indicator vector (all ones if all are death events)
#     event_indicators = ones(Int, length(death_times))

#     # Create a DataFrame with the required structure
#     data = DataFrame(time = death_times, event = event_indicators)

#     # Fit the Kaplan-Meier model
#     km = fit(KaplanMeier, data.time, data.event)

# # Extract survival probabilities and time points
#     survival_probabilities = km.survival
#     # Get unique event times from the original death_times data
#     unique_times = unique(sort(death_times))
#     # Ensure that the lengths of unique_times and survival_probabilities are the same
#     # This is necessary because the Kaplan-Meier method may have fewer points than unique_times
#     if length(unique_times) != length(survival_probabilities)
#         unique_times = unique_times[1:length(survival_probabilities)]
#     end


#  # Plot the Kaplan-Meier curve
#     plot!(plt,unique_times, survival_probabilities, label="$dose_spec")
   sorted_death_times = sort(death_times)

   # Initialize variables
   n = length(sorted_death_times)
   survival_probabilities = []
   times = []

   # Calculate survival probabilities
   for (i, time) in enumerate(sorted_death_times)
       if i == 1 || time != sorted_death_times[i-1]
           survival_probability = (n - i + 1) / n
           push!(survival_probabilities, survival_probability)
           push!(times, time)
       end
   end

   # Plot the Kaplan-Meier-type curve
   plot!(plt, times./24, survival_probabilities, label="$(dose_spec) Doses")

end

display(plt)

savefig(plt, "./results/$(drug)_kaplan_meier.svg")
savefig(plt, "./results/$(drug)_kaplan_meier.png")




























using Plots, DifferentialEquations, MAT, Statistics, Serialization, ModelingToolkit, Survival
gr()

# Assuming 'patients' is your vector of Patient structs

# Extract unique treatment conditions from all patients
conditions = [ "1.0e-5effect"
              ,"0.1effect"
              ,"0.25effect"
              ,"0.5effect"
              ,"0.75effect"
              ,"0.9effect"
              ,"0.99999effect"
              ,"random"
              ,"optimal"]

# Initialize an empty plot for the Kaplan-Meier curve
plt = plot(xlabel="Time (Days)", ylabel="Survival Probability", title="Kaplan-Meier Curve for Different Doses", grid=false)

# Iterate over each condition to perform survival analysis
for cond in conditions
    println("\rCondition: $cond")
    death_times = Float64[]

    for patient in patients
        if haskey(patient.output_measures, cond)
            # Assuming you can calculate or extract death time from the patient's data for the specific condition
            # For example, using a field from `patient.output_measures[cond]` or a custom function
            # Here, we're using a placeholder function `calculate_death_time` which you'll need to define
            death_time = calculate_death_time(patient.output_measures[cond])
            push!(death_times, death_time)
        end
    end

    # Sort death times and perform Kaplan-Meier analysis or equivalent
    # Here we simply plot sorted death times for illustration; adjust according to your survival analysis method
    sorted_death_times = sort(death_times)

    # Plot survival probabilities for this condition
    # Replace this with your actual method for calculating and plotting survival probabilities
    # This is a placeholder for illustration
    plot!(plt, sorted_death_times, label="$cond Doses")
end

display(plt)

# Save the plot
savefig(plt, "./results/kaplan_meier.svg")
savefig(plt, "./results/kaplan_meier.png")
