using Pkg
Pkg.activate(".")
Pkg.instantiate()
include("../utilities/utils.jl")

using Random 
using Distributions
using HypothesisTests
using Statistics

include("../scripts/setup/generate_vp.jl")
include("../scripts/setup/precompute_scale.jl")
include("../assets/pk_params.jl")


function test_scale_precompute(drug)
    # Parameters for the test
    num_patients = 1000  # Large enough sample size for statistical tests
    seed = 123

    include("../model/$(drug)_pkpd2.jl")
    include("../model/$(drug)_dosing2.jl")
    include("../model/$(drug)_params.jl")
    include("../scripts/setup/init_integrate.jl")

    drug_params = eval(Symbol(uppercase(drug) * "_params"))

    # Combining all parameter dictionaries
    all_params = merge(TMZ_params, drug_params)

    # Generate population
    population = generate_virtual_population(TMZ_params, drug_params, ode_params, param_order, num_patients, seed)

    max_dose_min_tumor, min_dose_max_tumor = compute_loss_scaling(population, doses, fill(min_drug_dosage/length(doses), length(doses)))

    nb_patients = size(population, 2)
    # Validation checks (adjust as necessary)
    if length(max_dose_min_tumor) != nb_patients || length(min_dose_max_tumor) != nb_patients
        error("Output array sizes are incorrect.")
    end

    if any(isnan.(max_dose_min_tumor)) || any(isnan.(min_dose_max_tumor))
        error("NaN values found in output.")
    end

    # Validation check: max tumor size should be larger than min tumor size
    if any(max_dose_min_tumor .> min_dose_max_tumor)
        error("Validation failed: Max tumor size is not consistently larger than min tumor size.")
    end

    "Test passed: Scaling precompute function output is valid and max tumor size is larger than min tumor size"
end

# Call the test function with a specific drug
test_result = test_scale_precompute("rg")  # Replace with actual drug name
println(test_result)