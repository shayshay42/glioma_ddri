using Pkg
prefix = pwd()#"/scratch/c/craigm/popsicle/gbm/DosageOptimization.jl"
Pkg.activate(prefix)
Pkg.resolve()
Pkg.instantiate()

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

drug = "rg"
filename = prefix*"/results/optima_$(drug)_lfbgs_auc_maxinit_$(today_date)_dms.jls"
#loads the PK/PD model
include(prefix*"/models/$(drug)_pkpd2.jl")
#loads the dosing schedule and amount
include(prefix*"/models/$(drug)_dosing2.jl")
#loads the PK/PD parameters
include(prefix*"/models/$(drug)_params.jl")
#loads the utility functions
include(prefix*"/utilities/utils.jl")
# include("../../assets/generate_vp.jl")
population = deserialize(prefix*"/assets/$(drug)_vp4.jls")
tumor_min_scaling = deserialize(prefix*"/assets/$(drug)_max_dose_min_tumor_v4.jls")
tumor_max_scaling = deserialize(prefix*"/assets/$(drug)_min_dose_max_tumor_v4.jls")
#make a voctor of tuples of the form (x,y)
scaling = [collect(s) for s in zip(tumor_min_scaling, tumor_max_scaling)]
# drug_min_scaling = deserialize(prefix*"/assets/$(drug)_min_dose_auc.jls")
# drug_max_scaling = deserialize(prefix*"/assets/$(drug)_max_dose_auc.jls")
drug_min_scaling = 20.0*avg_human_surface_area*length(rg_dosetimes)#min_drug_dosage
drug_max_scaling = 1800.0*avg_human_surface_area*length(rg_dosetimes)#sum(doses)

num_patients = size(population, 2)

@info "Loaded Packages and Snippets."

u0 = zeros(Float64, 9)
u0[1] = 17.7
tspan = (0.0, end_time + 7.0) .* hours
dosetimes = eval(Symbol("$(drug)_dosetimes"))
# time_dose_map = Dict(zip(dosetimes, Int64.(1:length(dosetimes))))
ode_params = [population[:,1]; scaling[1]]#; time_dose_map]

# p = zeros(Float64, length(ode_params) + length(doses))
# p[1:length(ode_params)] .= ode_params
# p[length(ode_params)+1:end] .= doses
p = [ode_params; doses]

prob = ODEProblem(pk_pd!, u0, tspan, p)
#sol = solve(prob, saveat=0.1)
#plotter(sol)

sys = modelingtoolkitize(prob)
simp = structural_simplify(sys)
prob_jac = ODEProblem(simp, u0, tspan, [population[:,1];scaling[1];doses], jac=true)
#solj = solve(prob_jac, saveat=0.1, callback=hit)
#plotter(solj)

function loss(θ, ode_params)
    int_sol = solve(prob_jac, nothing, p=[ode_params;θ], callback=hit, sensealg=nothing)
    cell = (int_sol[end,end]-ode_params[end-1])/(ode_params[end]-ode_params[end-1])
    drug = (sum(abs,θ)-drug_min_scaling)/(drug_max_scaling-drug_min_scaling)
    return  cell + drug
end

#TESTING THE SCALING OF THE LOSS FUNCTION

# doses = fill(1800.0*avg_human_surface_area, length(rg_dosetimes))
# # doses = fill(20.0*avg_human_surface_area, length(rg_dosetimes))
# int_sol = solve(prob_jac, nothing, p=[ode_params;doses], callback=hit, sensealg=nothing)
# int_sol[end,end]
# ode_params[end-1]
# ode_params[end]
# loss(doses, ode_params)
#rng = Random.default_rng()
#Random.seed!(rng, 1)
#dose_init = rand(length(doses)).*dose_amount

# scaler = Array{Float64, 2}(undef,  1000, 1)
# Threads.@threads for i in 1:1000
#     println("\rPatient: $i")
#     ode_params = [population[:,i]; scaling[i]]
#     int_sol = solve(prob_jac, nothing, p=[ode_params;doses], callback=hit, sensealg=nothing)
#     scaler[i] = int_sol[end,end]
# end
# serialize(prefix*"/assets/maxrg&TMZ_ftv_scaler.jls", scaler)
algo=nothing
sense=nothing
#θ = dose_init

loss(doses, ode_params)
#ForwardDiff.gradient(doses -> loss(doses, ode_params), doses)

#adtype = Optimization.AutoForwardDiff()
#optf = Optimization.OptimizationFunction((x, p) -> loss(x, p), adtype)
#optprob = Optimization.OptimizationProblem(optf, doses, ode_params, lb=zeros(length(doses)), ub=doses)
#t_iter=100
#callback, losses, iter, anim = create_callback(t_iter,verbose=false, animate=true, progress_bar=true, saving=false, early_stopping=true)
#opt = Optim.LBFGS()
#res = Optimization.solve(optprob, opt, maxiters=t_iter, callback=callback)
#gif(anim, "results/$(drug)_auc_boxconstraint_test_lfbgs.gif", fps=10)
#sole = solve(prob_jac, nothing, p=[ode_params;res.u], callback=hit, sensealg=nothing)
#plotter(sole)

#optimization
@info "Starting optimization"

#optima = zeros(length(doses)+1, num_patients)
optima = Array{Any,2}(undef, length(doses)+2, num_patients)

adtype = Optimization.AutoForwardDiff()
optf = Optimization.OptimizationFunction((x,p) -> loss(x, p), adtype)
opt = Optim.LBFGS()
t_iter=100

# Initialize an Atomic counter
completed_patients = Threads.Atomic{Int64}(0)

# Generate a random permutation of indices from 1 to 1000 and select the first 500
random_indices = randperm(length(1:num_patients))[1:500]
# Now loop over these randomly selected indices
# Threads.@threads for i in random_indices
Threads.@threads for i in 1:num_patients
    println("Processing patient $i, Completed patients: $(completed_patients[])")
    ode_params = [population[:,i];scaling[i]]

    # gradient descent optimization of the dose amounts
    optprob = Optimization.OptimizationProblem(optf, doses, ode_params, lb=zeros(length(doses)), ub=doses)
    callback, losses, iter_ref = create_callback(t_iter,verbose=false, animate=false, progress_bar=false, saving=false, early_stopping=true)
    res = Optimization.solve(optprob, opt, maxiters=t_iter, callback=callback)
    optima[:, i] = [res.u;iter_ref[];[losses]]
    
    # Increment and print the number of patients that have been completed
    Threads.atomic_add!(completed_patients, 1)

    # Save the optima array every 2 patients
    if (completed_patients[] % 50 == 0 && completed_patients[] > 0)
        temp_filename = prefix*"/results/optima_$(drug)_$(completed_patients[])_$(today_date).jls"
        open(temp_filename, "w") do file
            serialize(file, optima)
        end
        println("Saved optima after patient $(completed_patients[])")
        @info "Saved optima after patient $(completed_patients[])"
    end
end

open(filename, "w") do file
    serialize(file, optima)
end;
