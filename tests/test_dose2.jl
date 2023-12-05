using Pkg
Pkg.activate(".")
# Pkg.instantiate()
include("../utilities/utils.jl")

#holds the NLMEM parameters for the PK model
include("../assets/pk_params.jl")

#inlcude the functions that are used in setting up the struct
include("../scripts/setup/generate_vp.jl")
include("../scripts/setup/compute_dose_2.jl")

drug = "rg"

#loads the model equations, dosing event functions and the parameter values
include("../model/$(drug)_pkpd2.jl")
include("../model/$(drug)_params.jl")

function test_dose_matrix(num_patients, seed, drug)

    drug_params = eval(Symbol(uppercase(drug) * "_params"))
    # Generate population
    population = generate_virtual_population(TMZ_params, drug_params, ode_params, param_order, num_patients, seed)
    @info "Computing dose matrix..."
    ic50_fractional_doses = compute_dose2(population, drug_params)
    return ic50_fractional_doses
end

doses = test_dose_matrix(1000, 123, drug)

using Plots
 
#make a heat map of the doses
heatmap(doses')