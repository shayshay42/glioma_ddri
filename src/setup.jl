using Pkg
Pkg.activate(".")
Pkg.instantiate()

#add the utils function mainly logit and erelu used by functions in this file
include("../utilities/utils.jl")
#holds the NLMEM parameters for the PK model
include("../assets/pk_params.jl")

#specify the struct that holds the patient data
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

#inlcude the functions that are used in setting up the struct
include("../scripts/setup/generate_vp.jl")
include("../scripts/setup/precompute_scale.jl")

function generate_patients_struct(num_patients, seed, drug)
    patients = Vector{Patient}(undef, num_patients)

    #loads the model equations, dosing event functions and the parameter values
    include("../model/$(drug)_pkpd2.jl")
    include("../model/$(drug)_dosing2.jl")
    include("../model/$(drug)_params.jl")
    include("../scripts/setup/init_integrate.jl")

    drug_params = eval(Symbol(uppercase(drug) * "_params"))

    # Generate population
    population = generate_virtual_population(TMZ_params, drug_params, ode_params, param_order, num_patients, seed)

    max_dose_min_tumor, min_dose_max_tumor = compute_loss_scaling(population, doses, fill(min_drug_dosage/length(doses), length(doses)))


NEED TO ADD A FUNCTION TO COMPUTE THE TRAJECTORY AND OUTPUTS FOR ALL DOSING CONDITIONS THAT ARENT THE OPTIMAL ONES (THIS IS WHERE IC50 FRACTIONS COME IN SO MOVE TO PROTOTYPE)



    # put ODE parameters in the struct
    for i in 1:num_patients
        patients[i].idx = i
        patients[i].ode_parameters = population[:, i]
        patients[i].minimal_tumor = max_dose_min_tumor[i]
        patients[i].maximal_tumor = min_dose_max_tumor[i]
    end

    
end