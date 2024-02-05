using Pkg
Pkg.activate(".")
Pkg.instantiate()

#add the utils function mainly logit and erelu used by functions in this file
include("../utilities/utils.jl")
#holds the NLMEM parameters for the PK model
include("../assets/pk_params.jl")

#specify the struct that holds the patient data
struct ConditionalOutput
    doses::Vector{Float64}
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
    # optimal_doses::Vector{Float64}
    output_measures::Dict{String, ConditionalOutput}
end

#inlcude the functions that are used in setting up the struct
include("../scripts/setup/generate_vp.jl")
include("../scripts/setup/compute_dose_bvp.jl")
include("../scripts/setup/precompute_scale.jl")
include("../scripts/setup/compute_outputs.jl")

function generate_patients_struct(num_patients, seed, drug; min_dose=true, gradation=[1e-5, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0-1e-5])
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
    maxi_doses = doses_matrix[end,:]
    if min_dose
        mini_doses = zeros(num_patients)
    else
        mini_doses = doses_matrix[1,:]
    end
    max_dose_min_tumor, min_dose_max_tumor = compute_loss_scaling_effect_based(population, maxi_doses, mini_doses, drug_dose_times)
    scaling = [collect(s) for s in zip(max_dose_min_tumor, min_dose_max_tumor, maxi_doses, mini_doses)]

    #want the average dose per patient
    avg_dose_per_gradation = mean(doses_matrix, dims=2)
        #the loss here is meaningless unless
    # print("avg_dose_per_gradation: ", avg_dose_per_gradation)
    # print("min avg doses", ones(length(num_patients)).*avg_dose_per_gradation[end])
    maxi_avg_doses = ones(num_patients).*avg_dose_per_gradation[end]
    if min_dose
        mini_avg_doses = zeros(num_patients)
    else
        mini_avg_doses = ones(num_patients).*avg_dose_per_gradation[1]
    end
    avg_max_dose_min_tumor, avg_min_dose_max_tumor = compute_loss_scaling_effect_based(population, maxi_avg_doses, mini_avg_doses, drug_dose_times)
    avg_scale = [collect(s) for s in zip(avg_max_dose_min_tumor, avg_min_dose_max_tumor, maxi_avg_doses, mini_avg_doses)]
               
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
    Random.seed!(seed)
    Threads.@threads for i in 1:num_patients
        @info "Generating Patient $i"
        output_measures = Dict()

        for (j, effect) in enumerate(effect_keys)
            output_measures[effect] = get_outputs(population[:,i], ones(length(drug_dose_times)).*doses_matrix[j,i], scaling[i], drug)
        end

        for (j, dose) in enumerate(avg_dose_per_gradation)
            output_measures[string(gradation[j])*"avg_effect"] = get_outputs(population[:,i], ones(length(drug_dose_times)).*dose, avg_scale[i], drug)
        end

        min_dose = doses_matrix[1,i]
        max_dose = doses_matrix[end,i]
        random_doses = min_dose .+ (max_dose - min_dose) .* rand(length(drug_dose_times))
        output_measures["random"] = get_outputs(population[:,i], random_doses, scaling[i], drug)

        # Create a new Patient instance
        # patients[i] = Patient(i, population[:, i], scaling[i], zeros(length(drug_dose_times)), output_measures)
        patients[i] = Patient(i, population[:, i], scaling[i], output_measures)
    end
    @info "Finished generating patients."
    
    # OPTIMIZE

    # # put ODE parameters in the struct
    # for i in 1:num_patients
    #     patients[i].optimal_doses = optimal_doses[:,i]
    #     patients[i].output_measures["optimal"] = get_outputs(population[:,i], optimal_doses[:,i], scaling[i], drug)
    # end
    return patients
end
