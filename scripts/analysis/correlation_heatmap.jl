using Plots, DifferentialEquations, MAT, Statistics, Serialization, ModelingToolkit
using StatsBase
pyplot()

include("../../utilities/utils.jl")

drug = "gdc"
include("../../models/$(drug)_pkpd2.jl")
include("../../models/$(drug)_dosing2.jl")
include("../../models/$(drug)_params.jl")

# population = deserialize("./assets/$(drug)_vp4.jls")
# min_tumor = deserialize("./assets/$(drug)_min_tumor_max_dose.jls")

# #deserialize optimum .jls from results
# optimum = deserialize("./results/optimum/optima_$(drug)_lfbgs_auc_oct27_dms.jls")
# #mkae it into a Float64 array
# opt = Array{Float64, 2}(optimum[1:end-2,:])
     # for heatmap plotting
gr()            
if drug == "rg"
    drug_params = ["Cl2", "ka2", "Vpla", "Q", "Vtis"]
else
    drug_params = ["ka2", "V2", "kel", "k12", "k21"]
end
TMZ_params = ["Cl1", "k23", "ka1"]

indices = indexin([TMZ_params;drug_params], param_order)
labels = [TMZ_params;drug_params]
# select GR backend for Plots
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
# Load the vector of 1000 Patient structs from a file (assumed to be serialized)
filename = "./assets/$(drug)_patients_struct_v5.jls" # Update this to the correct path
patients = deserialize(filename) # or use load(filename) for JLD2

# Compute correlations
conditions = ["Optim", "Min", "Quarter", "Half", "Max", "Random"]
correlation_matrices = Dict{String, Matrix{Float64}}()

for condition in conditions
    param_values = hcat([p.ode_parameters[indices] for p in patients]...)
    loss_values = [p.output_measures[condition].loss for p in patients]
    ftv_values = [p.output_measures[condition].ftv for p in patients]
    drug_auc_values = [p.output_measures[condition].drug_auc for p in patients]
    tumor_auc_values = [p.output_measures[condition].tumor_auc for p in patients]
    optimal_dose_sums = [sum(p.optimal_doses) for p in patients]

    # Calculate correlations and handle NaNs
    correlation_matrix = hcat(
        cor(param_values', loss_values),
        cor(param_values', ftv_values),
        cor(param_values', drug_auc_values),
        cor(param_values', tumor_auc_values),
        cor(param_values', optimal_dose_sums)
    )
    
    correlation_matrices[condition] = correlation_matrix
end

# Before you start plotting, define the consistent color limits for your heatmaps
# You may need to adjust `min_color_val` and `max_color_val` according to your data range.
min_color_val = -1.0 # Assuming the correlation can be between -1 and 1
max_color_val = 1.0

# Make sure the results directory exists
results_dir = "./results/heatmaps/$(drug)"
isdir(results_dir) || mkdir(results_dir)

outs = ["Loss", "Final Tumor Volume", "Drug Exposure", "Tumor Burden", "Dose Sum"]
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
    save_heatmap(matrix, labels, svg_filename, condition)
end






# pyplot()

# outs = ["Loss", "Final Tumor Volume", "Drug Exposure", "Tumor Burden", "Dose Sum"]
# function create_sexy_heatmap(matrix, labels, title)
#     # Define the size of the plot and the font size
#     default(size=(800, 600), fontfamily="DejaVu Sans", legendfontsize=10, guidefontsize=12, tickfontsize=10, titlefontsize=14)

#     # Create the heatmap with a modern color palette and additional styling
#     heatmap = Plots.heatmap(
#         outs, labels, matrix, 
#         clims=(min_color_val, max_color_val), 
#         color=:viridis, # This is a modern, perceptually-uniform color map
#         aspect_ratio=:equal, # Keeps cells square
#         xticks=(1:length(outs), outs), 
#         yticks=(1:length(labels), labels), 
#         xrotation=45, # Rotate x-axis labels for better legibility
#         title=title, # Add a title to the plot
#         titlefontsize=16, # Adjust title font size
#         xlabel="Outputs",
#         ylabel="Parameters",
#         size=(600, 500)
#     )
    
#     # Adding colorbar with label
#     # colorbar!(title="Correlation")

#     return heatmap
# end


# # Assuming you have already defined your `correlation_matrices` and `labels` as before
# # Create a sexy heatmap for the first condition as an example
# first_condition = first(conditions)
# first_matrix = correlation_matrices[first_condition]
# first_title = "$(drug) - $(first_condition)"


# # Your heatmap code here
# heatmap(["Loss", "Final Tumor Volume", "Drug Exposure", "Tumor Burden", "Dose Sum"],
#         labels, first_matrix, 
#         clims=(min_color_val, max_color_val), 
#         color=:bam, # This is a modern, perceptually-uniform color map
#         aspect_ratio=:equal)

# # If the xrotation keyword is not working, you can try modifying the xticks attribute directly
# xticks = get(current(), :xticks)
# set!(xticks.rotation = 45) # This syntax may vary depending on the package and backend


# heatmap(["Loss", "Final Tumor Volume", "Drug Exposure", "Tumor Burden", "Dose Sum"],
#         labels, first_matrix,
#         clims=(min_color_val, max_color_val),
#         color=:viridis,
#         aspect_ratio=:equal,
#         xticks=(1:length(outs), outs),  # Ensure all x-tick labels are to be shown
#         xrotation=45)  
# # Here is how you might rotate the x-axis labels
# PyPlot.xticks(rotation=45)
# # Create the heatmap
# sexy_heatmap = create_sexy_heatmap(first_matrix, labels, first_title)

# # To show the heatmap in a REPL or Jupyter Notebook
# display(sexy_heatmap)
