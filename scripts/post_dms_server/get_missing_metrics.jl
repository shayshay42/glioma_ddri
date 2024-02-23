# filename = "results/optim/rg_probmask_test_ADAM_finitediff_temp2_lr0.1_2024-02-08_400patients.jls"
filename = "results/optim/rg_probmask_2AUCloss_ADAM_finitediff_temp2_lr0.1_2024-02-21_100patients.jls"

using Dates, StatsPlots, Plots, Serialization, OrderedCollections
today = Dates.today()

drug = "rg"
#
include("../../model/$(drug)_pkpd2.jl")
include("../../model/$(drug)_dosing2.jl")
include("../../model/$(drug)_params.jl")
include("../../scripts/setup/init_integrate.jl")
include("../../assets/pk_params.jl")
include("../../scripts/setup/generate_vp_lhs.jl")
include("../../scripts/setup/compute_effect_dose_optim.jl")
include("../../scripts/setup/precompute_scale.jl")
include("../../scripts/setup/compute_outputs.jl")
include("../../src/setup.jl")
include("../../utilities/utils.jl")

# filename = "results/optim/rg_probmask_2AUCloss_test_ADAM_finitediff_temp2_lr0.1_2024-02-17_400patients.jls"
patients = deserialize(open(filename, "r"))

num_patients = length(patients)
vp_param_matrix = hcat([patients[i].ode_parameters for i in 1:num_patients]...)
gradation = [1e-5, 0.001, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0-1e-5]
patient_doses = compute_dose_optim(vp_param_matrix, gradation)
#show the unstable patients
# println(patient_doses.retcodes)
# check that non are above 1e-6 loss
for i in 1:num_patients
    for j in size(patient_doses.mins, 1)
        if patient_doses.mins[j,i] > 1e-6
            println("Patient $(i) has a loss of $(patient_doses.mins[j,i])")
        end
    end
end





dose_matrix = patient_doses.doses_matrix

maxi_doses = dose_matrix[end,:]
num_patients = length(patients)
mini_doses = zeros(num_patients)

maxi_doses = ones(num_dose_times).*(single_max)



# if effect
#     # mini_doses = doses_matrix[1,:]
#     max_dose_min_tumor, min_dose_max_tumor = compute_loss_scaling_effect_based(population, maxi_doses, mini_doses, drug_dose_times)
#     scaling = [collect(s) for s in zip(max_dose_min_tumor, min_dose_max_tumor, maxi_doses, mini_doses)]
# else
include("../../scripts/setup/init_integrate.jl")
vp_param_matrix = hcat([patients[i].ode_parameters for i in 1:num_patients]...)
max_dose_min_tumor, min_dose_max_tumor, max_dose_max_drug, min_dose_min_drug = compute_loss_scaling_AUC(vp_param_matrix, maxi_doses, mini_doses, drug)
scaling = [collect(s) for s in zip(max_dose_min_tumor, min_dose_max_tumor, max_dose_max_drug, min_dose_min_drug)]
# end
effect_keys = [string(i)*"effect" for i in gradation]
avg_dose_per_gradation = mean(dose_matrix, dims=2)
#turns out max effect dose is higher dose than single max
maxi_avg_doses = ones(num_patients).*avg_dose_per_gradation[end]
avg_max_dose_min_tumor, avg_min_dose_max_tumor, avg_max_dose_max_drug, avg_min_dose_min_drug = compute_loss_scaling_AUC(vp_param_matrix, maxi_avg_doses, mini_doses, drug)
avg_scale = [collect(s) for s in zip(avg_max_dose_min_tumor, avg_min_dose_max_tumor, avg_max_dose_max_drug, avg_min_dose_min_drug)]

#compute the metrics from the doses computed
drug_dose_times = eval(Symbol(drug * "_dosetimes"))
repaired_patients = Vector{Patient2}(undef, num_patients)
Threads.@threads for i in 1:num_patients
    @info "Completing Patient $i"
    output_measures = OrderedDict()

    for (j, effect) in enumerate(effect_keys)
        output_measures[effect] = get_outputs_auc(patients[i].ode_parameters, ones(length(drug_dose_times)).*dose_matrix[j,i], scaling[i], drug)
    end

    for (j, dose) in enumerate(avg_dose_per_gradation)
        output_measures[string(gradation[j])*"avg_effect"] = get_outputs_auc(patients[i].ode_parameters, ones(length(drug_dose_times)).*dose, avg_scale[i], drug)
    end

    min_dose = dose_matrix[1,i]
    max_dose = dose_matrix[end,i]
    random_doses = min_dose .+ (max_dose - min_dose) .* rand(length(drug_dose_times))
    output_measures["random_effect"] = get_outputs_auc(patients[i].ode_parameters, random_doses, scaling[i], drug)
    output_measures["random"] = get_outputs_auc(patients[i].ode_parameters, rand(length(drug_dose_times)).*single_max, avg_scale[i], drug)
    output_measures["none"] = get_outputs_auc(patients[i].ode_parameters, zeros(length(drug_dose_times)), scaling[i], drug)
    output_measures["min"] = get_outputs_auc(patients[i].ode_parameters, ones(length(drug_dose_times)).*single_min, scaling[i], drug)
    output_measures["half"] = get_outputs_auc(patients[i].ode_parameters, ones(length(drug_dose_times)).*((single_max-single_min)/2), scaling[i], drug)
    output_measures["max"] = get_outputs_auc(patients[i].ode_parameters, ones(length(drug_dose_times)).*single_max, scaling[i], drug)
    # Create a new Patient instance
    # patients[i] = Patient(i, population[:, i], scaling[i], zeros(length(drug_dose_times)), output_measures)
    #recover the optimals
    output_measures["optimal"] = get_outputs_auc(patients[i].ode_parameters, patients[i].output_measures["optimal"].doses, scaling[i], drug)

    # patients[i].output_measures = output_measures
    repaired_patients[i] = Patient2(i, patients[i].ode_parameters, scaling[i], output_measures)
end

#serialize the repaired patients
filename = "results/optim/$(drug)_dms_100patients_2AUCloss_2024-02-23_richoutputsmeasures.jls"

open(filename, "w") do file
    serialize(file, repaired_patients)
end