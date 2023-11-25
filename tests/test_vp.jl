using Pkg
Pkg.activate(".")
Pkg.instantiate()
include("../utilities/utils.jl")

using Random 
using Distributions
using HypothesisTests
using Statistics


include("../scripts/setup/generate_vp.jl")
include("../assets/pk_params.jl")


function test_virtual_population_generation(drug = "rg")
    # Parameters for the test
    num_patients = 1000  # Large enough sample size for statistical tests
    seed = 123

    include("../model/$(drug)_params.jl")

    drug_params = eval(Symbol(uppercase(drug) * "_params"))

    # Combining all parameter dictionaries
    all_params = merge(TMZ_params, drug_params)

    # Generate population
    population = generate_virtual_population(TMZ_params, drug_params, ode_params, param_order, num_patients, seed)

    # Check if the dimensions are correct
    if size(population, 1) != length(param_order) || size(population, 2) != num_patients
        return "Dimension mismatch: Expected dimensions (length(param_order), num_patients), got $(size(population))"
    end

    feedback = []

    # Check each parameter distribution
    for (i, param) in enumerate(param_order)
        if param in keys(all_params)
            mu, sigma, mult = all_params[param]
            expected_dist = param == "Q" ? LogitNormal(logit(mu), sigma) : LogNormal(log(mu), sigma)
            # ks_test = HypothesisTests.ExactOneSampleKSTest(population[i, :]./mult, expected_dist)
            # p_value = pvalue(ks_test)
            ad_test = HypothesisTests.OneSampleADTest(population[i, :]./mult, expected_dist)
            p_value = pvalue(ad_test)
            if p_value < 0.05
                push!(feedback, "Distribution mismatch for parameter '$param' (p-value: $p_value)")
            else
                push!(feedback, "Parameter '$param' matches expected distribution (p-value: $p_value)")
            end
        end
    end

    # If all tests pass
    if isempty(feedback)
        return ""
    else
        return join(feedback, "\n")
    end
end

println(test_virtual_population_generation("rg"))