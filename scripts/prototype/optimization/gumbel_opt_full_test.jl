include("../../../src/setup.jl")
#loads the model equations, dosing event functions and the parameter values
drug = "rg"

include("../../../model/$(drug)_pkpd2.jl")
include("../../../model/$(drug)_dosing2.jl")
include("../../../model/$(drug)_params.jl")
include("../../../scripts/setup/init_integrate.jl")
include("../../../utilities/utils.jl")


using DifferentialEquations, LinearAlgebra, ModelingToolkit

using Random, Serialization
using SciMLSensitivity
using Zygote, Optimization, OptimizationOptimisers, OptimizationOptimJL
using Logging

using StaticArrays
using ForwardDiff
using FiniteDiff
using Flux
using BenchmarkTools

using Dates
today_date = Dates.today()

struct optimum_result
    dose_vector::Vector{Float64}
    iter_num::Int64
    loss_trajectory::Vector{Float64}
end

const mult_num = 20
const dose_mults = [0:mult_num...]
const unit_dose = single_max/mult_num
const possible_doses = [0:unit_dose:single_max...]
const temp = 2
#should schedule the temperature to decrease over time (phasing)

# an issue arises in setting the pill constant dose to vary with the maximum amount of the drug given based on the highest effect
# to keep the loss function form going outside of 0-2 range we have to scale the loss not by effect but by plausible min and max doses
# otherwise can ask a different question: we know some patients are more sensitive since their max effect will be reached at a lower dose
# But now we want to know what the optimal way to dose them over time is. No need to use multiples of the unit dose that end up above the maximal effect dose
# could keep the distribution 0 after a certain multiple closest to the max effect dose
# or in reality thte optimla dose may be allowed to go above the maximal effect dose? Since maximal effect dose is computed as the dose when given 18 times reaching the maximal effect
# but one a single day the dose should be allowed to reach higher than that.
# solution/decision: away from varying constant dose based on the patient. But new problem is ensuring that the loss is still defined in this scenario

function loss(λ, ode_params, scaling)
    # const_dose = scaling[3]/mult_num
    drug_min_scaling = scaling[4]*num_dose_times
    drug_max_scaling = scaling[3]*num_dose_times

    # gumbel_noise = -log.(-log.(rand(length(set), num)))
    gumbel_noise = -log.(-log.(rand(length(possible_doses), num_dose_times)))
    sample_pmf = exp.((λ + gumbel_noise) ./ temp)
    sample_pmf ./= sum(sample_pmf, dims=1)
    # θ = ((sample_pmf' * dose_mults) .* const_dose)'
    θ = (sample_pmf' * possible_doses)'

    sol_temp = solve(prob, Rodas4P2(), p=[ode_params..., θ...], callback=hit, sensealg=nothing)
    cell = (sol_temp[end,end] - scaling[1])/(scaling[2]-scaling[1])
    drug = (sum(abs,θ)-drug_min_scaling)/(drug_max_scaling-drug_min_scaling)
    loss = cell + drug
    return loss
end

function prob_to_dose(prob)
    final_pmf = exp.(prob ./ temp)
    final_pmf ./= sum(final_pmf, dims=1)
    doses = (repeat(possible_doses,1,num_dose_times)[argmax(final_pmf, dims=1)])
    return doses
end

num_patients = 400
seed = 123

patients = generate_patients_struct(num_patients, seed, drug)

drug_dose_times = eval(Symbol(drug * "_dosetimes"))
#order patients vector by the .idx property of the elements
sort!(patients, by=x->x.idx)
@info "Starting optimization for $(uppercase(drug))-model ..."
optima = Vector{optimum_result}(undef, length(patients))
patients_optim = []

t_iter=2

# Initialize an Atomic counter
completed_patients = Threads.Atomic{Int64}(0)
Random.seed!(seed)

lr=0.1

function adam!(logits, grad, m,v, i, lr=0.1, beta1=0.9, beta2=0.999, epsilon=1e-8)
    m = beta1 * m + (1 - beta1) * grad
    v = beta2 * v + (1 - beta2) * grad.^2

    m_hat = m / (1 - beta1^i)
    v_hat = v / (1 - beta2^i)

    logits -= lr * m_hat ./ (sqrt.(v_hat) .+ epsilon)
end


Threads.@threads for i in 1:3#length(patients)
    println("Processing patient $i, Completed patients: $(completed_patients[])")
    #retrive element in patients vector with .idx property equal to i
    patient = patients[findfirst(x->x.idx==i, patients)]
    ode_p = patient.ode_parameters
    scaler = patient.scaling

    init_logits = ones(mult_num+1, num_dose_times)
    logits = deepcopy(init_logits)
    
    # loss_values = losses = []

    m = zeros(size(logits))
    v = zeros(size(logits))

    for j in 1:t_iter
        @info "Iter: $j"
        grad = ForwardDiff.gradient(λ -> loss(λ, ode_p, scaler), logits)
        adam!(logits, grad, m, v, j)
    end
    converted_doses = prob_to_dose(logits)
    optima[i] = optimum_result(converted_doses, t_iter, [])
    patient.output_measures["optimal"] = get_outputs(ode_p, converted_doses, scaler, drug)
    
    # Increment and print the number of patients that have been completed
    Threads.atomic_add!(completed_patients, 1)

    # Save the optima array every 2 patients
    if (completed_patients[] % 50 == 0 && completed_patients[] > 0)
        temp_filename = "results/optim/$(drug)_probmask_test_ADAM_finitediff_temp$(temp)_lr$(lr)_$(today_date)_completion__$(completed_patients[]).jls"
        #  "../results/optima_$(drug)_$(completed_patients[])_$(today_date).jls"
        open(temp_filename, "w") do file
            serialize(file, optima)
        end
        println("Saved optima after patient $(completed_patients[])")
        @info "Saved optima after patient $(completed_patients[])"
    end
    push!(patients_optim, patient) 
end
    
#save the patients vector that contains the optimized patient structs
open("results/optim/$(drug)_probmask_test_ADAM_finitediff_temp$(temp)_lr$(lr)_$(today_date)_$(num_patients)patients.jls", "w") do file
    serialize(file, patients_optim)
end