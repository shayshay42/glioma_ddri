using Pkg
Pkg.activate(".")
Pkg.instantiate()
include("../utilities/utils.jl")

include("../scripts/setup/generate_vp.jl")
include("../assets/pk_params.jl")

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

function generate_patients_struct(num_patients, seed, drug)
    patients = Vector{Patient}(undef, num_patients)

    include("../model/$(drug)_params.jl")

    drug_params = eval(Symbol(uppercase(drug) * "_params"))

    # Generate population
    population = generate_virtual_population(TMZ_params, drug_params, ode_params, param_order, num_patients, seed)

    # put ODE parameters in the struct
    for i in 1:num_patients
        patients[i].idx = i
        patients[i].ode_parameters = population[:, i]
    end

    
end