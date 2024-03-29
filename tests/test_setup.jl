using Pkg
Pkg.activate(".")
Pkg.instantiate()

using Plots
using StatsPlots
gr()
# pyplot()


#add the utils function mainly logit and erelu used by functions in this file
include("../utilities/utils.jl")
#holds the NLMEM parameters for the PK model
include("../assets/pk_params.jl")

#specify the struct that holds the patient data
struct ConditionalOutput
    loss::Float64
    ftv::Float64
    drug_auc::Float64
    tumor_auc::Float64
    trajectory::Array{Float64,2}
end

struct Patient
    idx::Int
    ode_parameters::Vector{Float64}
    # minimal_tumor::Float64
    # maximal_tumor::Float64
    scaling::Vector{Float64}
    optimal_doses::Vector{Float64}
    output_measures::Dict{String, ConditionalOutput}
end

#inlcude the functions that are used in setting up the struct
include("../scripts/setup/generate_vp.jl")
include("../scripts/setup/compute_dose_bvp.jl")
include("../scripts/setup/precompute_scale.jl")
include("../scripts/setup/compute_outputs.jl")

function generate_patients_struct(num_patients, seed, drug; gradation=[1e-5, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0-1e-5])
    patients = Vector{Patient}(undef, num_patients)

    drug_params = eval(Symbol(uppercase(drug) * "_params"))
    drug_dose_times = eval(Symbol(drug * "_dosetimes"))
    # Generate population
    population = generate_virtual_population(TMZ_params, drug_params, ode_params, param_order, num_patients, seed)
    # gradation = [1e-5, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0 -1e-5]
    patient_doses = compute_dose_effect(population, gradation)
    # gradation = patient_doses.grad
    doses_matrix = patient_doses.doses_matrix
    # retcodes = patient_doses.retcodes
    # max_dose_min_tumor, min_dose_max_tumor = compute_loss_scaling(population, doses, fill(min_drug_dosage/length(doses), length(doses)))
    max_dose_min_tumor, min_dose_max_tumor = compute_loss_scaling_effect_based(population, doses_matrix[end,:], doses_matrix[1,:], drug_dose_times)
    scaling = [collect(s) for s in zip(max_dose_min_tumor, min_dose_max_tumor, doses_matrix[end,:], doses_matrix[1,:])]

    #the loss here is meaningless unless
    avg_max_dose_min_tumor, avg_min_dose_max_tumor = compute_loss_scaling_effect_based(population, doses_matrix[end,:], doses_matrix[1,:], drug_dose_times)
    avg_scale = [collect(s) for s in zip(avg_max_dose_min_tumor, avg_min_dose_max_tumor, avg_dose_per_gradation[end], avg_dose_per_gradation[1])]
                
    #want the average dose per patient
    avg_dose_per_gradation = mean(doses_matrix, dims=2)
    #store output in the struct
    effect_keys = [string(i)*"effect" for i in gradation]
    # get_outputs(population[:,k], ones(length(drug_dosetimes)).*avg_dose_per_gradation[1], scaling[k], drug)
    # Threads.@threads for i in 1:num_patients
    #     @info "Generating Patient $i"
    #     patients[i].idx = i
    #     patients[i].ode_parameters = population[:, i]
    #     patients[i].scaling = scaling[i]
    #     # patients[i].optimal_doses = doses_matrix[end,:]
    #     for (j, effect) in enumerate(effect_keys)
    #         patients[i].output_measures[effect] = get_outputs(population[:,i], ones(length(drug_dose_times)).*doses_matrix[j,i], scaling[i], drug)
    #     end
    #     for (j, dose) in enumerate(avg_dose_per_gradation)
    #         patients[i].output_measures[string(gradation[j])*"avg_effect"] = get_outputs(population[:,i], ones(length(drug_dose_times)).*dose, scaling[i], drug)
    #     end
    #     min_dose = doses_matrix[1,i]
    #     max_dose = doses_matrix[end,i]
    #     random_doses = min_dose .+ (max_dose - min_dose) .* rand(length(drug_dose_times))
    #     patients[i].output_measures["random"] = get_outputs(population[:,i], random_doses, scaling[i], drug)
    #     # patients[i] = Patient(i, population[:, i], max_dose_min_tumor[i], min_dose_max_tumor[i], OPTIMAL DOSE, DICTIONARY OF OUTPUTS)
    # end

    Threads.@threads for i in 1:num_patients
        @info "Generating Patient $i"
        output_measures = Dict()

        for (j, effect) in enumerate(effect_keys)
            output_measures[effect] = get_outputs(population[:,i], ones(length(drug_dose_times)).*doses_matrix[j,i], scaling[i], drug)
        end

        for (j, dose) in enumerate(avg_dose_per_gradation)
            output_measures[string(gradation[j])*"avg_effect"] = get_outputs(population[:,i], ones(length(drug_dose_times)).*dose, avg_scale, drug)
        end

        min_dose = doses_matrix[1,i]
        max_dose = doses_matrix[end,i]
        random_doses = min_dose .+ (max_dose - min_dose) .* rand(length(drug_dose_times))
        output_measures["random"] = get_outputs(population[:,i], random_doses, scaling[i], drug)

        # Create a new Patient instance
        patients[i] = Patient(i, population[:, i], scaling[i], zeros(length(drug_dose_times)), output_measures)
    end

    
    # OPTIMIZE

    # # put ODE parameters in the struct
    # for i in 1:num_patients
    #     patients[i].optimal_doses = optimal_doses[:,i]
    #     patients[i].output_measures["optimal"] = get_outputs(population[:,i], optimal_doses[:,i], scaling[i], drug)
    # end
    return patients
end

num_patients = 200
seed = 123
drug = "rg"

#loads the model equations, dosing event functions and the parameter values
include("../model/$(drug)_pkpd2.jl")
include("../model/$(drug)_dosing2.jl")
include("../model/$(drug)_params.jl")
include("../scripts/setup/init_integrate.jl")

patients = generate_patients_struct(num_patients, seed, drug)

#plot viiolins

conditions = keys(patients[1].output_measures)
conditions = [
        "1.0e-5effect"
        ,"0.1effect"
        ,"0.25effect"
        ,"0.5effect"
        ,"0.75effect"
        ,"0.9effect"
        ,"0.99999effect"

        ,"1.0e-5avg_effect"
        ,"0.1avg_effect"
        ,"0.25avg_effect"
        ,"0.5avg_effect"
        ,"0.75avg_effect"
        ,"0.9avg_effect"
        ,"0.99999avg_effect"

        ,"random"
]
colors = fill(:gray, length(conditions)) # Define a color for each condition

ftv_data = OrderedDict(cond => [] for cond in conditions)
drug_auc_data = OrderedDict(cond => [] for cond in conditions)
tumor_auc_data = OrderedDict(cond => [] for cond in conditions)

for patient in patients
    for cond in conditions
        push!(ftv_data[cond], patient.output_measures[cond].ftv)
        push!(drug_auc_data[cond], patient.output_measures[cond].drug_auc)
        push!(tumor_auc_data[cond], patient.output_measures[cond].tumor_auc)
    end
end

function create_combined_plot(metric_data, title)
    p = plot(title=title, legend=false, grid=false)

    for (i, cond) in enumerate(conditions)
        # Create box plot with reduced width and darker fill
        # Adjust the width parameter here to control the box width
        box_width = 0.05 # Adjust this value as needed
        # boxplot!(p, [cond], metric_data[cond], width=box_width, color=darker_colors[i], linecolor=:black, fillalpha=0.3, outliers_marker=:asterisk, outliers_color=:red, label=false)

        # Overlay with violin plot
        violin!(p, [cond], metric_data[cond], color=colors[i], alpha=0.7, label=false, xrotation=45)
    end
    # Rotate x-tick labels
    # plot!(p, xticks=(1:length(conditions), conditions), xrotation=45)
    return p
end

p2 = create_combined_plot(ftv_data, "FTV Distribution")
p3 = create_combined_plot(drug_auc_data, "Drug AUC Distribution")
p4 = create_combined_plot(tumor_auc_data, "Tumor AUC Distribution")

# Function to save plot in multiple formats
function save_plot(plot, base_filename)
    savefig(plot, "./results/outputs_violin/$(drug)" * base_filename * ".png")
    # savefig(plot, "./results/violin/$(drug)" * base_filename * ".svg")
end

# Save the plots
save_plot(p2, "FTV_Distribution")
save_plot(p3, "Drug_AUC_Distribution")
save_plot(p4, "Tumor_AUC_Distribution")