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
    ode_p = 
    @info "Starting optimization for $(uppercase(drug))-model ..."
    optima = Vector{optimum_result}(undef, length(patients))

    adtype = Optimization.AutoForwardDiff()
    optf = Optimization.OptimizationFunction((x,p,s) -> loss(x, p, s), adtype)
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
        optprob = Optimization.OptimizationProblem(optf, init_doses, ode_p, scaler, lb=ones(length(drug_dose_times)).*min_tested, ub=ones(length(drug_dose_times)).*max_tested)
        callback, losses, iter_ref = create_callback(t_iter,verbose=false, animate=false, progress_bar=false, saving=false, early_stopping=true)
        res = Optimization.solve(optprob, opt, maxiters=t_iter, callback=callback)
        optima[i] = optimum_result(res.u, iter_ref[], losses)
        
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
    end
end