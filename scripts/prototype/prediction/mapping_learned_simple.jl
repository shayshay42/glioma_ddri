
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




using LinearAlgebra

# Correctly transform the data for linear regression
# Assuming `ode_params` and `optimal_doses` need to be flattened or transformed appropriately

# Create new columns in df for flattened or transformed data
# Ensure these transformations match your data's needs
df[!, :flattened_params] = [reduce(vcat, p) for p in df[!, :ode_params]]  # Adjust as needed
df[!, :flattened_doses] = [reduce(vcat, d) for d in df[!, :optimal_doses]]  # Adjust as needed

# 1000-element Vector{Vector{Float64}}
# Assuming you flatten or process the data into a single vector per patient if needed


flattened_params = df[!, :flattened_params]
flattened_doses = df[!, :flattened_doses]

# Assuming all vectors have the same length, use the first element to determine sizes
num_patients = length(flattened_params)
num_params = length(first(flattened_params))
num_doses = length(first(flattened_doses))

# Initialize matrices
X = Matrix{Float64}(undef, num_patients, num_params)
Y = Matrix{Float64}(undef, num_patients, num_doses)

# Fill in the matrices
for i in 1:num_patients
    X[i, :] = flattened_params[i]
    Y[i, :] = flattened_doses[i]
end

X = hcat(ones(num_patients), X)  # Add a column of ones as the first column

using MultivariateStats

# λ is the regularization parameter; you might need to adjust this based on your data
λ = 0.01

# Add λI to X'X before inversion
β_ridge = (X'X + λ * I) \ X'Y

#-----------------MLP-----------------

using Flux, MLDataUtils

# Flatten the optimal doses if they are not already in a 1D format
# This step might need adjustments based on your actual data structure
# optimal_doses_flattened = reduce(vcat, optimal_doses)

# Ensure ODE parameters are in a suitable format for training
# Assuming each set of parameters is already a flat vector
X = hcat(ode_params...)'
Y = hcat(optimal_doses...)'

# Splitting the data (70% train, 15% validation, 15% test)
(train_indices, val_test_indices) = splitobs(collect(1:size(X, 1)), at=0.7)
(val_indices, test_indices) = splitobs(val_test_indices, at=0.5)

X_train, Y_train = X[train_indices, :], Y[train_indices, :]
X_val, Y_val = X[val_indices, :], Y[val_indices, :]
X_test, Y_test = X[test_indices, :], Y[test_indices, :]

X_train = Float32.(X_train')
X_val = Float32.(X_val')
X_test = Float32.(X_test')
Y_train = Float32.(Y_train')
Y_val = Float32.(Y_val')
Y_test = Float32.(Y_test')

# Neural network model
# Neural network model
model = Chain(
    Dense(size(X, 2), 4, relu),  # First layer with ReLU activation
    Dense(4, size(Y, 2))          # Output layer
)

model1 = Chain(
    Dense(size(X, 2), size(Y, 2))          # Output layer
)

model2 = Chain(
    Dense(size(X, 2), 32, relu),  # First layer with ReLU activation
    Dense(32, size(Y, 2))          # Output layer
)

model3 = Chain(
    Dense(size(X, 2), 32, relu),  # First layer with ReLU activation
    Dense(32, 32, relu),
    Dense(32, size(Y, 2))          # Output layer
)

# Loss function (mean squared error)

# Optimizer (you can adjust learning rate and other parameters as needed)


# Data for training (converting to batches might be necessary for larger datasets)

data = [(X_train, Y_train)]

# Training the model
epochs = 10000 # Adjust the number of epochs as necessary
predictions_per_model = Dict{String, Array{Float32, 2}}()
for m in [model, model1, model2, model3]
    val_losses = []
    train_losses = []
    loss(x, y) = Flux.mse(m(x), y)
    optimizer = ADAM(0.001)
    for epoch = 1:epochs
        Flux.train!(loss, Flux.params(m), data, optimizer)

        # Optionally, evaluate on validation set periodically
        if epoch % 10 == 0
            val_loss = loss(X_val, Y_val)
            push!(val_losses, val_loss)
            println("Epoch: $epoch, Validation Loss: $val_loss")
            train_loss = loss(X_train, Y_train)
            # println("Epoch: $epoch, Training Loss: $train_loss")
            push!(train_losses, train_loss)
        end
    end
    predictions_per_model[string(m)] = m(X_test)
end


# Select a patient for visualization
patient_idx = 1

# Create a plot for each model's predictions vs actual doses
p = plot(title="Actual vs. Predicted Doses for Patient $patient_idx", legend=:outertopright)
for (name, preds) in predictions_per_model
    plot!(p, preds[patient_idx, :], label=name)
end
plot!(p, Y_test[patient_idx, :], label="Actual Doses", lw=3)

# Display the plot
display(p)

# Create a histogram of prediction errors for each model
p_errors = plot(title="Distribution of Prediction Errors", xlabel="Error", ylabel="Frequency", legend=:outertopright)
for (name, preds) in predictions_per_model
    errors = vec(Y_test .- preds)
    histogram!(p_errors, errors, bins=50, alpha=0.5, label=name, normed=true)
end

# Add ridge regression errors to the same plot
errors_ridge = vec(Y_test .- Y_pred_ridge)
histogram!(p_errors, errors_ridge, bins=50, alpha=0.5, label="Ridge Regression", normed=true)

# Display the plot
display(p_errors)






# Assuming you have the trained MLP model and the Ridge Regression coefficients β_ridge
# Ensure X_test includes the intercept term for Ridge Regression
X_test_ridge = vcat(ones(1,size(X_test, 2)), X_test)  # Add intercept term

# Predictions using Ridge Regression
Y_pred_ridge = X_test_ridge' * β_ridge

# Predictions using MLP
# Y_pred_mlp = model(X_test)  # Ensure X_test is correctly preprocessed for the MLP model

Y_pred_mlp = predictions_per_model[string(model)]
Y_pred_mlp1 = predictions_per_model[string(model1)]
Y_pred_mlp2 = predictions_per_model[string(model2)]
Y_pred_mlp3 = predictions_per_model[string(model3)]


# Calculate MSE for both models
mse_ridge = mean((Y_test - Y_pred_ridge').^2)
mse_mlp = mean((Y_test - Y_pred_mlp).^2)
mse_mlp1 = mean((Y_test - Y_pred_mlp1).^2)
mse_mlp2 = mean((Y_test - Y_pred_mlp2).^2)
mse_mlp3 = mean((Y_test - Y_pred_mlp3).^2)

println("MSE Ridge Regression: $mse_ridge")
println("MSE MLP: $mse_mlp")
println("MSE MLP1: $mse_mlp1")
println("MSE MLP2: $mse_mlp2")
println("MSE MLP3: $mse_mlp3")

# create a bar plot of the values
bar(["Ridge Regression", "MLP", "MLP1", "MLP2", "MLP3"], [mse_ridge, mse_mlp, mse_mlp1, mse_mlp2, mse_mlp3], 
    title="Mean Squared Error Comparison", xlabel="Model", ylabel="MSE", legend=false)


    # ... [previous code remains unchanged]


# Initialize a dictionary to store the MSE of each model
mse_values = Dict()
mae_values = Dict()
r2_values = Dict()

# Calculate MSE, MAE, and R² for the Ridge Regression model
mse_ridge = mean((Y_test - Y_pred_ridge').^2)
mae_ridge = mean(abs.(Y_test - Y_pred_ridge'))
r2_ridge = 1 - sum((Y_test - Y_pred_ridge').^2) / sum((Y_test .- mean(Y_test)).^2)
mse_values["Ridge Regression"] = mse_ridge
mae_values["Ridge Regression"] = mae_ridge
r2_values["Ridge Regression"] = r2_ridge

# Print the metrics
println("Ridge Regression - MSE: $mse_ridge, MAE: $mae_ridge, R²: $r2_ridge")

# Calculate MSE, MAE, and R² for each MLP model and store in the dictionary
for (name, preds) in predictions_per_model
    mse = mean((Y_test - preds).^2)
    mae = mean(abs.(Y_test - preds))
    r2 = 1 - sum((Y_test - preds).^2) / sum((Y_test .- mean(Y_test)).^2)
    mse_values[name] = mse
    mae_values[name] = mae
    r2_values[name] = r2
    
    # Print the metrics
    println("$name - MSE: $mse, MAE: $mae, R²: $r2")
end

# Create a bar plot for the MSE of each model
bar_names = collect(keys(mse_values))
bar_values = collect(values(mse_values))
bar(bar_names, bar_values, title="Mean Squared Error Comparison", xlabel="Model", ylabel="MSE", legend=false, rotation=10, size=(1080, 700))

# Create a bar plot for the MAE of each model
bar_names = collect(keys(mae_values))
bar_values = collect(values(mae_values))
bar(bar_names, bar_values, title="Mean Absolute Error Comparison", xlabel="Model", ylabel="MAE", legend=false, rotation=10, size=(1080, 700))

# Create a bar plot for the R² of each model
bar_names = collect(keys(r2_values))
bar_values = collect(values(r2_values))
bar(bar_names, bar_values, title="R² Comparison", xlabel="Model", ylabel="R²", legend=false, rotation=10, size=(1080, 700))





using Flux, Statistics
using Plots
using Distributions

models = [model, model1, model2, model3]
model_names = ["4", "0", "32", "2x32", "Ridge"]

# Initialize arrays to store MSEs for each model
mse_distributions = Dict(name => Float64[] for name in model_names)

# Number of training repetitions
num_repetitions = 10  # Adjust based on your computational budget

# Training loop for each model
for (i, m) in enumerate(models)
    println("Training $(model_names[i])")
    for n in 1:num_repetitions
        # Reset model parameters if needed
        Flux.reset!(m)
        epochs = 5000 # Adjust the number of epochs as necessary
        val_losses = []
        train_losses = []
        loss(x, y) = Flux.mse(m(x), y)
        optimizer = ADAM(0.001)
        for epoch = 1:epochs
            Flux.train!(loss, Flux.params(m), data, optimizer)
    
            # Optionally, evaluate on validation set periodically
            if epoch % 10 == 0
                val_loss = loss(X_val, Y_val)
                push!(val_losses, val_loss)
                println("Epoch: $epoch, Validation Loss: $val_loss")
                train_loss = loss(X_train, Y_train)
                # println("Epoch: $epoch, Training Loss: $train_loss")
                push!(train_losses, train_loss)
            end
        end
        predictions_per_model[string(m)] = m(X_test)
        y_pred = m(X_test)
        mse = mean((Y_test - y_pred).^2)

        # Store MSE
        push!(mse_distributions[model_names[i]], mse)
    end
end

# Statistical analysis
p_values = Dict()
for i in 1:length(models)-1
    for j in i+1:length(models)
        p_value = ttest_ind(mse_distributions[model_names[i]], mse_distributions[model_names[j]]).pvalue
        p_values["$(model_names[i]) vs $(model_names[j])"] = p_value
    end
end

# Visualization
bar_plot_data = [[(mean(mse), std(mse)) for mse in values(mse_distributions)];(mse_ridge,0.0)]
bar_heights = first.(bar_plot_data)
bar_errors = last.(bar_plot_data)

bar(model_names, bar_heights, yerr=bar_errors, legend=false)
for (i, height) in enumerate(bar_heights)
    text(model_names[i], height, string(round(height, digits=2)), valign=:bottom)
end

# Optionally, annotate plot with p-values
# ...

display(plot)




using Plots
using Statistics
using HypothesisTests

# Example data (use your actual mse_distributions and p_values)
model_names = ["Model", "Model1", "Model2", "Model3"]
mse_averages = [mean(mse) for mse in values(mse_distributions)]
mse_stds = [std(mse) for mse in values(mse_distributions)]

# Plotting the bar graph with error bars
p = bar(model_names, mse_averages, yerr=mse_stds, legend=false, title="MSE Comparison",
        ylabel="Mean Squared Error", xlabel="Model", color=:lightblue)

# Annotating p-values
for ((model1, model2), p_value) in p_values
    # Calculate positions
    idx1 = findfirst(isequal(model1), model_names)
    idx2 = findfirst(isequal(model2), model_names)
    y_position = max(mse_averages[idx1], mse_averages[idx2]) + max(mse_stds...) * 0.1  # Adjust as necessary

    # Drawing a line for visual linking
    plot!([idx1, idx2], [y_position, y_position], line=(:black, :dash, 1), legend=false)

    # Adding p-value text
    annotate!((mean([idx1, idx2]), y_position, text("p=$(round(p_value, digits=4))", 8, :center, :bottom)))
end

display(p)


