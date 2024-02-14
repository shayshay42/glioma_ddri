filename = "results/optim/rg_probmask_test_ADAM_finitediff_temp2_lr0.1_2024-02-08_400patients.jls"
patients = deserialize(open(filename, "r"))

#


# Extracting measures and parameters
losses = [patient.output_measures["optimal"].loss for patient in patients] 
ftvs = [patient.output_measures["optimal"].ftv for patient in patients]    
drug_aucs = [patient.output_measures["optimal"].drug_auc for patient in patients]
tumor_aucs = [patient.output_measures["optimal"].tumor_auc for patient in patients]

if drug == "rg"
    cl2_values = [patient.ode_parameters[indexin(["Cl2"], param_order)[1]] for patient in patients]
    vpla_values = [patient.ode_parameters[indexin(["Vpla"], param_order)[1]] for patient in patients]
    criteria = [("Loss", losses), ("FTV", ftvs), ("Drug AUC", drug_aucs), ("Tumor AUC", tumor_aucs), ("Cl2", cl2_values), ("Vpla", vpla_values)]
else
    kel_values = [patient.ode_parameters[indexin(["kel"], param_order)[1]] for patient in patients]
    k23_values = [patient.ode_parameters[indexin(["k23"], param_order)[1]] for patient in patients]
    criteria = [("Loss", losses), ("FTV", ftvs), ("Drug AUC", drug_aucs), ("Tumor AUC", tumor_aucs), ("kel", kel_values), ("k23", k23_values)]
end

opt = reduce(hcat, [patient.output_measures["optimal"].doses for patient in patients])
discrete_palette = reverse(cgrad(:bilbao, 20, categorical=true))
for (criterion_name, criterion_values) in criteria
    # Sorting the patients based on the criterion
    sorted_indices = sortperm(criterion_values)
    sorted_opts = opt[:, sorted_indices]

    # Generate heatmap
    heatmap_title = "Sorted by " * criterion_name
    heatmap_filename_png = "./results/drug_matrix/$(drug)_heatmap_" * lowercase(criterion_name) * "_.png"
    heatmap_filename_svg = "./results/drug_matrix/$(drug)_heatmap_" * lowercase(criterion_name) * "_.svg"
    # cgrad(:matter, 10, categorical = true)palette(:acton10, 5)
    # Creating the heatmap
    heatmap_plot = heatmap(sorted_opts', color=discrete_palette, xlabel="Day", ylabel="Patient", title=heatmap_title, size=(600, 1200))

    # Adding annotations
    num_patients = size(sorted_opts, 2)
    highest_value_patient = num_patients  # Last in sorted array
    lowest_value_patient = 1  # First in sorted array

    annotate!(12, highest_value_patient+7, text("High " * criterion_name, 8, :right, :top))
    annotate!(12, lowest_value_patient-7, text("Low " * criterion_name, 8, :right, :bottom))

    # Save the heatmap
    savefig(heatmap_plot, heatmap_filename_png)
    savefig(heatmap_plot, heatmap_filename_svg)
end