using Pkg
Pkg.activate(".")
Pkg.instantiate()

#holds the NLMEM parameters for the PK model
include("../assets/pk_params.jl")
include("../utilities/utils.jl")

#inlcude the functions that are used in setting up the struct
include("../scripts/setup/generate_vp.jl")
include("../scripts/setup/compute_dose_bvp.jl")
include("../scripts/setup/precompute_scale.jl")

function generate_population(drug, num_patients, seed)
    # Load the model and parameters
    include("../model/$(drug)_pkpd2.jl")
    include("../model/$(drug)_dosing2.jl")
    include("../model/$(drug)_params.jl")

    drug_dosetimes =  eval(Symbol(drug * "_dosetimes"))
    drug_vparams = eval(Symbol(uppercase(drug) * "_params"))

    # Generate and return the population
    population = generate_virtual_population(TMZ_params, drug_vparams, ode_params, param_order, num_patients, seed)
    return population, drug_dosetimes
end

function compute_doses(population, gradation)
    @info "Computing dose matrix via bvp..."
    patient_doses = compute_dose_effect(population, gradation)
    gradation = patient_doses.grad
    doses_matrix = patient_doses.doses_matrix
    retcodes = patient_doses.retcodes
    return gradation, doses_matrix, retcodes
end

num_patients = 1000
seed = 123
drug = "rg"

# Generate population
population, drug_dosetimes = generate_population(drug, num_patients, seed)

# Compute doses
gradation = [0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999]
gradation, doses_matrix, retcodes = compute_doses(population, gradation)
# gotta 
max_doses = doses_matrix[end,:]

# println(retcodes)

include("../scripts/setup/init_integrate.jl")
# Compute loss scaling
max_dose_min_tumor, min_dose_max_tumor = compute_loss_scaling_effect_based(population, max_doses, drug_dosetimes)


function simulate_to_plot(amount_per_mult_per_patient, population, state=2, mults=gradation, drug_dose_times=drug_dosetimes)
    sol = solve(prob, Rodas4P2(), callback=hit, saveat=0.1)
    timepoints = sol.t
    sols = zeros(length(mults), size(amount_per_mult_per_patient, 2), length(timepoints))
    for (j,mult) in enumerate(mults)
        Threads.@threads for i in 1:size(amount_per_mult_per_patient, 2)
            dose = amount_per_mult_per_patient[j,i]
            tmp_prob = remake(prob, p=[population[:,i]..., default_scaling..., ones(length(drug_dose_times)).*dose...])
            tmp_sol = solve(tmp_prob, Rodas4P2(), callback=hit, saveat=0.1)
            sols[j,i,:] = tmp_sol[state,:]
        end
    end
    return timepoints, sols
end

state = 7
states = OrderedDict(zip(["C", "D", "AbsTMZ", "PlaTMZ", "CSFTMZ", "AbsRG", "PlaRG", "TisRG", "cAUC"], 1:9))
state_name = states.keys[state]
times, sols = simulate_to_plot(doses_matrix, population, states[state_name], gradation, drug_dosetimes)
for (i,mult) in enumerate(gradation)
    data_matrix = sols[i,:,:]
    timepoints = times ./24

    title = "VP Trajectories - $drug model $state_name - $(mult) effect"
    means = vec(mean(data_matrix, dims=1))

    if state_name == "C"
        limit_y = [0, 40]
    else
        limit_y = nothing
    end
    # Plot the trajectories with a label for the legend
    plot(timepoints, data_matrix'[:,1]     , alpha=0.3, linecolor=:lightgray, grid=false, label="$state_name", legend=:topleft, ylim=limit_y)
    plot!(timepoints, data_matrix'[:,2:end], alpha=0.3, linecolor=:lightgray, grid=false, label=false)
    # Plot the mean and standard deviation with labels for the legend
    plot!(timepoints, means, 
        ribbon=(vec(std(data_matrix, dims=1)), vec(std(data_matrix, dims=1))), 
        fillalpha=0.5, linewidth=3, linecolor=:pink, fillcolor=:lightblue, 
        label="Mean ± SD")

    # Set the title and axis labels
    title!(title)
    xlabel!("Time (days)")
    if state_name == "C"
        ylabel!("Volume (mL)")
    else
        ylabel!("Amount (mg)")
    end
    savefig("./results/$(num_patients)_$(drug)_model_$(mult)_effect_spaghetti_$state_name.png")
end

