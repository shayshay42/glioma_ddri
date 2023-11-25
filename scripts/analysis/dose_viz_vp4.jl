using DifferentialEquations, LinearAlgebra, ModelingToolkit
using Random, Serialization
using Plots
using StatsPlots
using Colors, ColorSchemes

gr()

drug = "gdc"
#loads the PK/PD model
include("../../models/$(drug)_pkpd2.jl")
#loads t"../..g schedule and amount
include("../../models/$(drug)_dosing2.jl")
#loads t"../.. parameters
include("../../models/$(drug)_params.jl")
#loads t"../..ty functions
include("../../utilities/utils.jl")
# include("../../assets/generate_vp.jl")
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
    minimal_tumor::Float64
    maximal_tumor::Float64
    optimal_doses::Vector{Float64}
    output_measures::Dict{String, ConditionalOutput}
end
filename = "./assets/$(drug)_patients_struct_v5.jls" # Update this to the correct path
patients = deserialize(filename)

# Extracting measures and parameters
losses = [patient.output_measures["Optim"].loss for patient in patients] 
ftvs = [patient.output_measures["Optim"].ftv for patient in patients]    
drug_aucs = [patient.output_measures["Optim"].drug_auc for patient in patients]
tumor_aucs = [patient.output_measures["Optim"].tumor_auc for patient in patients]

if drug == "rg"
    cl2_values = [patient.ode_parameters[indexin(["Cl2"], param_order)[1]] for patient in patients]
    vpla_values = [patient.ode_parameters[indexin(["Vpla"], param_order)[1]] for patient in patients]
    criteria = [("Loss", losses), ("FTV", ftvs), ("Drug AUC", drug_aucs), ("Tumor AUC", tumor_aucs), ("Cl2", cl2_values), ("Vpla", vpla_values)]
else
    kel_values = [patient.ode_parameters[indexin(["kel"], param_order)[1]] for patient in patients]
    k23_values = [patient.ode_parameters[indexin(["k23"], param_order)[1]] for patient in patients]
    criteria = [("Loss", losses), ("FTV", ftvs), ("Drug AUC", drug_aucs), ("Tumor AUC", tumor_aucs), ("kel", kel_values), ("k23", k23_values)]
end


discrete_palette = cgrad(:imola, 10, categorical=true)
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


# # Sorting patients by a specific measure (e.g., losses)
# sorted_indices_by_loss = sortperm(losses)
# sorted_opt_by_loss = opt[:, sorted_indices_by_loss]

# # Sort by FTV
# sorted_indices_by_ftv = sortperm(ftvs)
# sorted_opt_by_ftv = opt[:, sorted_indices_by_ftv]

# # Sort by drug AUC
# sorted_indices_by_drug_auc = sortperm(drug_aucs)
# sorted_opt_by_drug_auc = opt[:, sorted_indices_by_drug_auc]

# # Sort by tumor AUC
# sorted_indices_by_tumor_auc = sortperm(tumor_aucs)
# sorted_opt_by_tumor_auc = opt[:, sorted_indices_by_tumor_auc]

# # Sort by Cl2
# sorted_indices_by_cl2 = sortperm(cl2_values)
# sorted_opt_by_cl2 = opt[:, sorted_indices_by_cl2]

# # Sort by Vpla
# sorted_indices_by_vpla = sortperm(vpla_values)
# sorted_opt_by_vpla = opt[:, sorted_indices_by_vpla]


# # Heatmap sorted by loss
# heatmap_loss = heatmap(sorted_opt_by_loss', aspect_ratio=:auto, color=:viridis, xlabel="Patient", ylabel="Day", title="Dosage Heatmap Sorted by Loss")

# # Heatmap sorted by FTV
# heatmap_ftv = heatmap(sorted_opt_by_ftv', aspect_ratio=:auto, color=:viridis, xlabel="Patient", ylabel="Day", title="Dosage Heatmap Sorted by FTV")

# # Heatmap sorted by drug AUC
# heatmap_drug_auc = heatmap(sorted_opt_by_drug_auc', aspect_ratio=:auto, color=:viridis, xlabel="Patient", ylabel="Day", title="Dosage Heatmap Sorted by Drug AUC")

# # Heatmap sorted by tumor AUC
# heatmap_tumor_auc = heatmap(sorted_opt_by_tumor_auc', aspect_ratio=:auto, color=:viridis, xlabel="Patient", ylabel="Day", title="Dosage Heatmap Sorted by Tumor AUC")

# # Heatmap sorted by Cl2
# heatmap_cl2 = heatmap(sorted_opt_by_cl2', aspect_ratio=:auto, color=:viridis, xlabel="Patient", ylabel="Day", title="Dosage Heatmap Sorted by Cl2")

# # Heatmap sorted by Vpla
# heatmap_vpla = heatmap(sorted_opt_by_vpla', aspect_ratio=:auto, color=:viridis, xlabel="Patient", ylabel="Day", title="Dosage Heatmap Sorted by Vpla")


# # Save the heatmaps as PNG and SVG
# savefig(heatmap_loss, "./heatmap_loss_subset.png")
# savefig(heatmap_loss, "./heatmap_loss_subset.svg")

# savefig(heatmap_ftv, "./heatmap_ftv_subset.png")
# savefig(heatmap_ftv, "./heatmap_ftv_subset.svg")

# savefig(heatmap_drug_auc, "./heatmap_drug_auc_subset.png")
# savefig(heatmap_drug_auc, "./heatmap_drug_auc_subset.svg")

# savefig(heatmap_tumor_auc, "./heatmap_tumor_auc_subset.png")
# savefig(heatmap_tumor_auc, "./heatmap_tumor_auc_subset.svg")

# savefig(heatmap_cl2, "./heatmap_cl2_subset.png")
# savefig(heatmap_cl2, "./heatmap_cl2_subset.svg")

# savefig(heatmap_vpla, "./heatmap_vpla_subset.png")
# savefig(heatmap_vpla, "./heatmap_vpla_subset.svg")

# function binned_colormap(data, colormap_name, num_bins)
#     # Create bins
#     bins = range(minimum(data), stop=maximum(data), length=num_bins+1)
    
#     # Assign each value to a bin
#     binned_data = cut(data, bins, labels=false)

#     # Get the colors from the colormap
#     cmap = get(ColorSchemes, colormap_name)
#     colors = [get_color(cmap, i/num_bins) for i in 1:num_bins]

#     # Assign colors based on bins
#     color_assignment = [colors[bin] for bin in binned_data]

#     return color_assignment
# end

# function generate_and_save_heatmap(sorted_opts, criterion_values, indices, title, filename)
#     # Binning the colors based on the criterion values
#     color_bins = binned_colormap(criterion_values, :tab10, 10)

#     # Creating the heatmap
#     heatmap_plot = heatmap(1:size(sorted_opts')[2], 1:size(sorted_opts')[1], sorted_opts', 
#                            color=color_bins[indices],
#                            aspect_ratio=:auto, xlabel="Patient", ylabel="Day",
#                            title=title)

#     # Annotating for high/low values
#     annotate_high_low!(heatmap_plot, criterion_values, indices)

#     # Saving the heatmap
#     savefig(heatmap_plot, "./$(filename).png")
#     savefig(heatmap_plot, "./$(filename).svg")
# end

# # List of criteria and their respective titles and filenames
# criteria = [("loss", losses, sorted_indices_by_loss, "Dosage Heatmap Sorted by Loss", "heatmap_loss_subset"),
#             ("ftv", ftvs, sorted_indices_by_ftv, "Dosage Heatmap Sorted by FTV", "heatmap_ftv_subset"),
#             ("drug_auc", drug_aucs, sorted_indices_by_drug_auc, "Dosage Heatmap Sorted by Drug AUC", "heatmap_drug_auc_subset"),
#             ("tumor_auc", tumor_aucs, sorted_indices_by_tumor_auc, "Dosage Heatmap Sorted by Tumor AUC", "heatmap_tumor_auc_subset"),
#             ("cl2", cl2_values, sorted_indices_by_cl2, "Dosage Heatmap Sorted by Cl2", "heatmap_cl2_subset"),
#             ("vpla", vpla_values, sorted_indices_by_vpla, "Dosage Heatmap Sorted by Vpla", "heatmap_vpla_subset")]

# # Looping through each criterion to generate and save heatmaps
# for (criterion_name, criterion_values, indices, title, filename) in criteria
#     sorted_opts = opt[:, indices]
#     generate_and_save_heatmap(sorted_opts, criterion_values, indices, title, filename)
# end
