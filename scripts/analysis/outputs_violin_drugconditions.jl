using DifferentialEquations, LinearAlgebra, ModelingToolkit
using Random, Serialization
using Plots
using StatsPlots

pyplot()

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
# Load the vector of 1000 Patient structs from a file (assumed to be serialized)
filename = "./assets/$(drug)_patients_struct_v5.jls" # Update this to the correct path
patients = deserialize(filename) # or use load(filename) for JLD2

# Define the conditions
conditions = ["Min", "Quarter", "Half", "Max", "Optim", "Random"]
colors = [:red, :green, :blue, :orange, :purple, :yellow] # Define a color for each condition

# Initialize dictionaries to store the data
loss_data = Dict(cond => [] for cond in conditions)
ftv_data = Dict(cond => [] for cond in conditions)
drug_auc_data = Dict(cond => [] for cond in conditions)
tumor_auc_data = Dict(cond => [] for cond in conditions)

# Extract data from patients
for patient in patients
    for cond in conditions
        push!(loss_data[cond], patient.output_measures[cond].loss)
        push!(ftv_data[cond], patient.output_measures[cond].ftv)
        push!(drug_auc_data[cond], patient.output_measures[cond].drug_auc)
        push!(tumor_auc_data[cond], patient.output_measures[cond].tumor_auc)
    end
end

# Function to plot histograms
function plot_histograms(data_dict, title)
    p = plot(legend=:topright, title=title)
    for (i, cond) in enumerate(conditions)
        histogram!(p, data_dict[cond], bins=30, alpha=0.5, label=cond, color=colors[i])
    end
    return p
end
# Function to create a violin plot for a given metric
# Define a pastel color palette
# Define a pastel color palette
pastel_colors = [:lightpink, :lightgreen, :lightgoldenrodyellow, :lavender, :lightblue]
# Modern color palette and their darker versions
modern_colors = [:dodgerblue, :coral, :mediumseagreen, :orchid, :goldenrod, :steelblue]
darker_colors = [:darkblue, :darkred, :darkgreen, :darkorchid, :darkgoldenrod, :darkslateblue]

# Function to create a combined box and violin plot with adjusted aesthetics
function create_combined_plot(metric_data, title)
    p = plot(title=title, legend=false, grid=false)

    for (i, cond) in enumerate(conditions)
        # Create box plot with reduced width and darker fill
        # Adjust the width parameter here to control the box width
        box_width = 0.05 # Adjust this value as needed
        # boxplot!(p, [cond], metric_data[cond], width=box_width, color=darker_colors[i], linecolor=:black, fillalpha=0.3, outliers_marker=:asterisk, outliers_color=:red, label=false)

        # Overlay with violin plot
        violin!(p, [cond], metric_data[cond], color=modern_colors[i], alpha=0.7, label=false)
    end

    return p
end
# Create combined plots for each metric
# Create combined plots for each metric
p1 = create_combined_plot(loss_data, "Loss Distribution")
p2 = create_combined_plot(ftv_data, "FTV Distribution")
p3 = create_combined_plot(drug_auc_data, "Drug AUC Distribution")
p4 = create_combined_plot(tumor_auc_data, "Tumor AUC Distribution")

# Function to save plot in multiple formats
function save_plot(plot, base_filename)
    savefig(plot, "./results/violin/$(drug)" * base_filename * ".png")
    savefig(plot, "./results/violin/$(drug)" * base_filename * ".svg")
end

# Save the plots
save_plot(p1, "Loss_Distribution")
save_plot(p2, "FTV_Distribution")
save_plot(p3, "Drug_AUC_Distribution")
save_plot(p4, "Tumor_AUC_Distribution")


function create_combined_plot(metric_data, title; log_scale=false)
    if log_scale
        p = plot(title=title, legend=false, grid=false, yscale=:log10)
    else
        p = plot(title=title, legend=false, grid=false)
    end

    for (i, cond) in enumerate(conditions)
        # Overlay with violin plot
        violin!(p, [cond], metric_data[cond], color=modern_colors[i], alpha=0.7, label=false)
    end

    return p
end

p3 = create_combined_plot(drug_auc_data, "Drug AUC Distribution", log_scale=true)
save_plot(p3, "Drug_AUC_Distribution")
