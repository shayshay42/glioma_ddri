using Random
using Distributions  # Required for distributions like LogNormal

"""
    generate_virtual_population(TMZ_params, drug_params, param_values, param_order, num_patients, seed)

# Description
Generate a population matrix with pharmacokinetic parameters for a specified number of patients.
This function will use the provided seed for random number generation to ensure reproducibility.

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
function generate_virtual_population(TMZ_params, drug_params, param_values, param_order, num_patients, seed)
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
        # Generate samples for each parameter
        sample = rand(LogNormal(log(TMZ_params[param].mu), TMZ_params[param].sigma), num_patients) .* TMZ_params[param].mult
        # Assign samples to the appropriate row in the population matrix
        population[TMZ_indices[i], :] = sample
    end

    # Populate the other drug parameters in the population matrix
    for (i, param) in enumerate(drug_param_names)
        # Special case for parameter "Q" using LogitNormal distribution
        if param == "Q"
            sample = rand(LogitNormal(logit(drug_params[param].mu), drug_params[param].sigma), num_patients) .* drug_params[param].mult
        else
            # Generate samples for each parameter using LogNormal distribution
            sample = rand(LogNormal(log(drug_params[param].mu), drug_params[param].sigma), num_patients) .* drug_params[param].mult
        end
        # Assign samples to the appropriate row in the population matrix
        population[drug_indices[i], :] = sample
    end

    # Return the generated population matrix
    return population
end
