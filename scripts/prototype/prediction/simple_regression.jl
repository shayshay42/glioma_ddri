using DifferentialEquations, LinearAlgebra, ModelingToolkit
using Random, Serialization
using Plots
using StatsPlots
# using Colors, ColorSchemes

gr()

drug = "rg"
include("../../../src/setup.jl")
include("../../../utilities/utils.jl")

#loads the model equations, dosing event functions and the parameter values
include("../../../model/$(drug)_pkpd2.jl")
include("../../../model/$(drug)_dosing2.jl")
include("../../../model/$(drug)_params.jl")
include("../../../scripts/setup/init_integrate.jl")

using DifferentialEquations, LinearAlgebra, ModelingToolkit

using Random, Serialization
using SciMLSensitivity
using Zygote, Optimization, OptimizationOptimisers, OptimizationOptimJL
using Logging

using StaticArrays
using ForwardDiff
using Flux
using BenchmarkTools

using Dates
today_date = Dates.today()
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
filename = "/Users/shayanhajhashemi/Documents/glioma_ddri/assets/gdc_patients_struct_v5.jls"#"../../../../assets/$(drug)_patients_struct_v4.jls" # Update this to the correct path
patients = deserialize(filename)

using DataFrames

# Assuming `patients` is an array of `Patient`
ode_params = [p.ode_parameters for p in patients]
optimal_doses = [p.optimal_doses for p in patients]

# Convert to DataFrame for easier manipulation
df = DataFrame(ode_params=ode_params, optimal_doses=optimal_doses)

using GLM

# Example: Flatten the data or select specific parameters for the regression
# This is a placeholder transformation. Adapt it to your specific needs.
df[:flattened_params] = flatten_each_row(df[:ode_params])
df[:flattened_doses] = flatten_each_row(df[:optimal_doses])

# Linear regression model
model = lm(@formula(flattened_doses ~ flattened_params), df)

using Plots

# Calculate predicted doses from the model
# Note: You will need to adapt this to your specific data structure and model
predicted_doses = predict(model, df)

# Assuming `optimal_doses` is a vector of observed doses
# (You may need to adjust this to match the structure of your `df`)
observed_doses = flatten_each_row(df[:optimal_doses])

# Plot observed vs. predicted doses
scatter(observed_doses, predicted_doses, label="Observed vs. Predicted", 
        xlabel="Observed Doses", ylabel="Predicted Doses", title="Linear Regression Fit")

# Add a line representing perfect predictions for reference
plot!([minimum(observed_doses), maximum(observed_doses)], 
      [minimum(observed_doses), maximum(observed_doses)], 
      label="Perfect Fit", color=:red, linestyle=:dash)

# Display the plot
display(plot)
