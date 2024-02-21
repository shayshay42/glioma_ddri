using ForwardDiff
using Distributions
using LinearAlgebra

# Set of possible values and parameters
set = [0:10...]
num = 2
logits = rand(length(set), num)
lr = 0.1
max_iter = 1000
temperature = 1
lower_bound = set[1]
upper_bound = set[end]

# Gradient descent loop
for i in 1:max_iter
    # Gumbel-Softmax sampling
    # gumbel_noise = -log.(-log.(rand(length(logits))))
    # sample_pmf = exp.((logits .+ gumbel_noise) ./ temperature) # Applying softmax
    # sample_pmf /= sum(sample_pmf)
    
    # # Compute weighted sum as the sample x
    # x = sum(sample_pmf .* set)
    
    # # Compute loss
    # loss = (x - 5)^2

    # Compute gradient
    grad = ForwardDiff.gradient(λ -> begin
                                        gumbel_noise = -log.(-log.(rand(length(set), num)))
                                        sample_pmf = exp.((λ + gumbel_noise) ./ temperature)
                                        # sample_pmf = exp.((λ) ./ temperature)
                                        sample_pmf ./= sum(sample_pmf, dims=1)
                                        x = sample_pmf' * set
                                        sum(abs2, x' - [4 5])
                                     end, logits)

    # Update logits
    logits -= lr * grad
    # Apply constraints
    for x in logits
        if x < lower_bound
            x = lower_bound
        elseif x > upper_bound
            x = upper_bound
        end
    end
end

# Convert final logits to PMF
final_pmf = exp.(logits ./ temperature)
final_pmf ./= sum(final_pmf, dims=1)

# Select the value with the highest probability
selected_values = repeat(set,1,num)[argmax(final_pmf, dims=1)]

println("Selected value: ", selected_values)
println("Probability distribution: ", final_pmf)
println("Loss at selected value: ", (selected_values - [4 5])^2)


mvpmf = Array([[0.1,0.2,0.7], [0.3,0.5,0.2]])
set = [1,2,3]
mvset = fill(set, 2)
dot_product = dot(mvpmf, mvset)

# Create a 2x3 matrix where the rows sum to one
matrix1 = [0.0 1.0 0.0;
           0.0 0.0 1.0]

# Create a matrix with duplicate [1, 2, 3]
matrix2 = [1 2 3]

# Perform the dot product
# Transpose matrix2 to match dimensions for dot product
result = matrix1' * matrix2

#sum them 
sum(result, dims=1)