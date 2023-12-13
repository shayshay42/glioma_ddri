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
    return population, drug_dosetimes, drug_vparams
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
population, drug_dosetimes, drug_params = generate_population(drug, num_patients, seed)

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

# make heatmap of pearson correlation_heatmap
function compute_aucs(population, doses_matrix, gradation, drug_dose_times)
    @info "Computing AUCs..."
    drug_aucs = zeros(length(gradation), size(population, 2))
    cell_aucs = zeros(length(gradation), size(population, 2))
    for (j,mult) in enumerate(gradation)
        Threads.@threads for i in 1:size(population, 2)
            dose = doses_matrix[j,i]
            tmp_prob = remake(prob, p=[population[:,i]..., default_scaling..., ones(length(drug_dose_times)).*dose...])
            tmp_sol = solve(tmp_prob, Rodas4P2(), callback=hit, saveat=0.1)
            drug_aucs[j,i] = trapezoidal_rule(tmp_sol.t, tmp_sol[7,:])
            cell_aucs[j,i] = tmp_sol[end,end]
        end
    end
    return drug_aucs, cell_aucs
end

drug_aucs, cell_aucs = compute_aucs(population, doses_matrix, gradation, drug_dosetimes)

# compute the correlation bet ween the auc and the parameters across the population
pk_params = [TMZ_params.keys..., drug_params.keys...]
function compute_corr(drug_aucs, cell_aucs, gradation, population, doses_matrix, drug_dose_times)
    @info "Computing correlations..."
    # pk_param_names = collect(keys(pk_params))
    pk_indices = indexin(pk_params, param_order)
    pk_param_values= population[pk_indices, :]
    correlation_matrices = OrderedDict()
    for (j, mult) in enumerate(gradation)
        corr_mat = hcat(
            cor(pk_param_values', drug_aucs[j,:]),
            cor(pk_param_values', cell_aucs[j,:]),
        )
        correlation_matrices[mult] = corr_mat
    end
    return correlation_matrices
end
correlation_matrices = compute_corr(drug_aucs, cell_aucs, gradation, population, doses_matrix, drug_dosetimes)
# Before you start plotting, define the consistent color limits for your heatmaps
# You may need to adjust `min_color_val` and `max_color_val` according to your data range.
min_color_val = -1.0 # Assuming the correlation can be between -1 and 1
max_color_val = 1.0

# Make sure the results directory exists
results_dir = "./results/"

outs = ["Drug Exposure", "Tumor Burden"]
function create_sexy_heatmap(matrix, labels, title)
    # Define the size of the plot and the font size
    default(size=(800, 600), fontfamily="DejaVu Sans", legendfontsize=10, guidefontsize=12, tickfontsize=10, titlefontsize=14)

    # Create the heatmap with a modern color palette and additional styling
    heatmap = Plots.heatmap(
        outs, labels, matrix, 
        clims=(min_color_val, max_color_val), 
        color=cgrad(:RdBu, rev=true), # This is a modern, perceptually-uniform color map
        aspect_ratio=:equal, # Keeps cells square
        xrotation=45, # Rotate x-axis labels for better legibility
        title=title, # Add a title to the plot
        titlefontsize=16, # Adjust title font size
        xlabel="Outputs",
        ylabel="Parameters",
        size=(600, 500),
        grid=false,
        framestyle=:origin
    )
   
    # Adding colorbar with label
    # colorbar!(title="Correlation")

    return heatmap
end
# first_condition = first(conditions)
# first_matrix = correlation_matrices[first_condition]
# first_title = "$(drug) - $(first_condition)"
# sexy_heatmap = create_sexy_heatmap(first_matrix, labels, first_title)

# # To show the heatmap in a REPL or Jupyter Notebook
# display(sexy_heatmap)


# Function to save a heatmap
function save_heatmap(matrix, labels, filename, cond)
    heatmap = create_sexy_heatmap(matrix, labels, cond)

    # To show the heatmap in a REPL or Jupyter Notebook
    Plots.savefig(heatmap, filename * ".svg")
    Plots.savefig(heatmap, filename * ".png")
end

# Loop through the conditions and create + save heatmaps
for (condition, matrix) in correlation_matrices
    # Define the filenames for the heatmap
    svg_filename = joinpath(results_dir, "$(drug)_$(condition)_heatmap")
    # Create and save the heatmap
    save_heatmap(matrix, pk_params, svg_filename, condition)
end