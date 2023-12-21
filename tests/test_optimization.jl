include("../src/setup.jl")

using DifferentialEquations, LinearAlgebra, ModelingToolkit

using Random, Serialization
using SciMLSensitivity
using Zygote, Optimization, OptimizationOptimisers, OptimizationOptimJL
using Logging

using StaticArrays
using ForwardDiff
using Flux
using BenchmarkTools

using Dates
today_date = Dates.today()

struct optimum_result
    dose_vector::Vector{Float64}
    iter_num::Int64
    loss_trajectory::Vector{Float64}
end

function compute_optimal_doses(patients, drug)
    drug_dose_times = eval(Symbol(drug * "_dosetimes"))
    #order patients vector by the .idx property of the elements
    sort!(patients, by=x->x.idx)
    @info "Starting optimization for $(uppercase(drug))-model ..."
    optima = Vector{optimum_result}(undef, length(patients))
    patients_optim = []

    adtype = Optimization.AutoForwardDiff()
    opt = Optim.LBFGS()
    t_iter=100

    # Initialize an Atomic counter
    completed_patients = Threads.Atomic{Int64}(0)
    Threads.@threads for i in 1:num_patients
        println("Processing patient $i, Completed patients: $(completed_patients[])")
        #retrive element in patients vector with .idx property equal to i
        patient = patients[findfirst(x->x.idx==i, patients)]
        ode_p = patient.ode_parameters
        scaler = patient.scaling
        init_doses = patient.output_measures["0.99999effect"].doses
        # ode_params = [population[:,i];scaling[i]]
    
        # gradient descent optimization of the dose amounts
            # Redefine the optimization function for the current patient
        optf = Optimization.OptimizationFunction((x, _) -> loss(x, ode_p, scaler), adtype)
        optprob = Optimization.OptimizationProblem(optf, init_doses, lb=ones(length(drug_dose_times)).*min_tested, ub=ones(length(drug_dose_times)).*max_tested)
        callback, losses, iter_ref = create_callback(t_iter,verbose=false, animate=false, progress_bar=false, saving=false, early_stopping=true)
        res = Optimization.solve(optprob, opt, maxiters=t_iter, callback=callback)
        optima[i] = optimum_result(res.u, iter_ref[], losses)
        patient.output_measures["optimal"] = get_outputs(ode_p, res.u, scaler, drug)
        
        # Increment and print the number of patients that have been completed
        Threads.atomic_add!(completed_patients, 1)
    
        # Save the optima array every 2 patients
        if (completed_patients[] % 50 == 0 && completed_patients[] > 0)
            temp_filename = "../results/optima_$(drug)_$(completed_patients[])_$(today_date).jls"
            open(temp_filename, "w") do file
                serialize(file, optima)
            end
            println("Saved optima after patient $(completed_patients[])")
            @info "Saved optima after patient $(completed_patients[])"
        end
        push!(patients_optim, patient) 
        return optima, patients_optim
    end
end

num_patients = 20
seed = 123
drug = "rg"

#loads the model equations, dosing event functions and the parameter values
include("../model/$(drug)_pkpd2.jl")
include("../model/$(drug)_dosing2.jl")
include("../model/$(drug)_params.jl")
include("../scripts/setup/init_integrate.jl")

patients = generate_patients_struct(num_patients, seed, drug)

optima, patients = compute_optimal_doses(patients, drug)
