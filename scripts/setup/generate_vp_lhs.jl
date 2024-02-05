using Random
using Distributions  # Required for distributions like LogNormal
using Statistics

"""
    generate_virtual_population_LHS(TMZ_params, drug_params, param_values, param_order, num_patients, seed)

# Description
Generate a population matrix with pharmacokinetic parameters for a specified number of patients.
This function will use the provided seed for random number generation to ensure reproducibility.
This uses the Latin Hypercube Sampling (LHS) method to generate the virtual population.

# Arguments
- `TMZ_params`: Dictionary of TMZ parameters containing name as keys and a tuple of mu, sigma, multiplier.
- `drug_params`: Dictionary of other drug parameters.
- `param_values`: Base values for each parameter.
- `param_order`: Order of parameters.
- `num_patients`: Number of patients.
- `seed`: Seed for random number generation.

# Returns
A matrix representing the generated population parameters.

"""
function generate_virtual_population_LHS(TMZ_params, drug_params, param_values, param_order, num_patients, seed)
    fractions=[eps(); [0:1/num_patients:1...][2:end-1]; 1-eps()]
    # Set the random seed for reproducibility
    Random.seed!(seed)

    # Initialize the population matrix
    population = repeat(param_values, 1, num_patients)
    
    # Extract parameter names from the dictionaries
    TMZ_param_names = collect(keys(TMZ_params))
    drug_param_names = collect(keys(drug_params))

    # Find the indices of parameters in the order specified by param_order
    TMZ_indices = indexin(TMZ_param_names, param_order)
    drug_indices = indexin(drug_param_names, param_order)

    # Populate the TMZ parameters in the population matrix
    for (i, param) in enumerate(TMZ_param_names)
        # Generate samples for each parameter (SCALED TO MATCH UNITS)
        inv_cdfs = quantile.(LogNormal(log(TMZ_params[param].mu), TMZ_params[param].sigma), fractions)
        samples = shuffle([rand(Uniform(inv_cdfs[i], inv_cdfs[i+1])) for i in 1:N].* TMZ_params[param].mult)
        # Assign samples to the appropriate row in the population matrix
        population[TMZ_indices[i], :] = samples
    end

    # Populate the other drug parameters in the population matrix
    for (i, param) in enumerate(drug_param_names)
        # Special case for parameter "Q" using LogitNormal distribution
        if param == "Q"
            inv_cdfs = quantile.(LogitNormal(logit(drug_params[param].mu), drug_params[param].sigma), fractions)
            samples = shuffle([rand(Uniform(inv_cdfs[i], inv_cdfs[i+1])) for i in 1:N].* drug_params[param].mult)
        else
            # Generate samples for each parameter using LogNormal distribution
            inv_cdfs = quantile.(LogNormal(log(drug_params[param].mu), drug_params[param].sigma), fractions)
            samples = shuffle([rand(Uniform(inv_cdfs[i], inv_cdfs[i+1])) for i in 1:N].* drug_params[param].mult)
        end
        # Assign samples to the appropriate row in the population matrix
        population[drug_indices[i], :] = samples
    end

    # Return the generated population matrix
    return population
end

# function simple_latin_hypercube_sampling(μ, σ, N)
#     # Generate N samples from the Normal(μ,σ) distribution by getting the inverse CDF at evenly spaced quantiles
#     points = [eps(); [0:1/N:1...][2:end-1]; 1-eps()]
#     inv_cdfs = quantile(Normal(μ, σ), points)
    
#     # Sample in the paired ranges from the inverse CDFs uniformly
#     samples = [rand(Uniform(inv_cdfs[i], inv_cdfs[i+1])) for i in 1:N]
    
#     return shuffle(samples)
# end

# # Use the function
# samples = simple_latin_hypercube_sampling(0, 1,1000)
# scatter(zeros(length(samples)), samples)
