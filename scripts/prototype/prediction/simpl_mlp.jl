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

using DataFrames

# Assuming `patients` is an array of `Patient`
ode_params = [p.ode_parameters for p in patients]
optimal_doses = [p.optimal_doses for p in patients]

# Convert to DataFrame for easier manipulation
df = DataFrame(ode_params=ode_params, optimal_doses=optimal_doses)

using Pkg
Pkg.add("Flux")
Pkg.add("MLDataUtils")

using Flux, MLDataUtils

# Flatten the optimal doses if they are not already in a 1D format
# This step might need adjustments based on your actual data structure
optimal_doses_flattened = reduce(vcat, optimal_doses)

# Ensure ODE parameters are in a suitable format for training
# Assuming each set of parameters is already a flat vector
X = hcat(ode_params...)'
Y = hcat(optimal_doses_flattened...)'

# Splitting the data (70% train, 15% validation, 15% test)
(train_indices, val_test_indices) = splitobs(collect(1:size(X, 1)), at=0.7)
(val_indices, test_indices) = splitobs(val_test_indices, at=0.5)

X_train, Y_train = X[train_indices, :], Y[train_indices, :]
X_val, Y_val = X[val_indices, :], Y[val_indices, :]
X_test, Y_test = X[test_indices, :], Y[test_indices, :]

# Neural network model
model = Chain(
    Dense(size(X, 2), 32, relu),  # First layer with ReLU activation
    Dense(32, 32, relu),           # Second layer with ReLU activation
    Dense(32, size(Y, 2))          # Output layer
)

# Loss function (mean squared error)
loss(x, y) = Flux.mse(model(x), y)

# Optimizer (you can adjust learning rate and other parameters as needed)
optimizer = ADAM(0.001)

# Data for training (converting to batches might be necessary for larger datasets)
data = [(X_train, Y_train)]

# Training the model
epochs = 100  # Adjust the number of epochs as necessary
for epoch = 1:epochs
    Flux.train!(loss, params(model), data, optimizer)
    
    # Optionally, evaluate on validation set periodically
    if epoch % 10 == 0
        val_loss = loss(X_val, Y_val)
        println("Epoch: $epoch, Validation Loss: $val_loss")
    end
end

test_loss = loss(X_test, Y_test)
println("Test Loss: $test_loss")

