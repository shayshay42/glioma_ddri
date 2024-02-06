#NOTES:
# no need for lower and upper bounds since the logits can take any real value
# the maximum dose per patient is constant now and no longer effect based
# can make the optimization for each patient a fucntion of their maximal effect dose
# ADAM saved the days
# optimum is max_tested dose hypothesis: since the scaling is for a very high effect


include("../../../src/setup.jl")
include("../../../utilities/utils.jl")

drug = "rg"
#loads the model equations, dosing event functions and the parameter values
include("../../../model/$(drug)_pkpd2.jl")
include("../../../model/$(drug)_dosing2.jl")
include("../../../model/$(drug)_params.jl")
include("../../../scripts/setup/init_integrate.jl")


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

const mult_num = 100
const dose_mults = [0:5:single_max...]
const temp = 1

function loss(λ, ode_params, scaling)
    # dose_unit = scaling[3]/mult_num
    dose_unit = (single_max)/mult_num
    drug_min_scaling = scaling[4]*num_dose_times
    drug_max_scaling = scaling[3]*num_dose_times

    gumbel_noise = -log.(-log.(rand(length(dose_mults), num_dose_times)))
    sample_pmf = exp.((λ + gumbel_noise) ./ temp)
    sample_pmf ./= sum(sample_pmf, dims=1)
    θ = ((sample_pmf' * dose_mults) .* dose_unit)'
    # println(θ)

    sol_temp = solve(prob, Rodas4P2(), p=[ode_params..., θ...], callback=hit, sensealg=nothing)
    cell = (sol_temp[end,end] - scaling[1])/(scaling[2]-scaling[1])
    drug = (sum(abs,θ)-drug_min_scaling)/(drug_max_scaling-drug_min_scaling)
    loss = cell + drug
    return loss
end

num_patients = 10
seed = 123

patients = generate_patients_struct(num_patients, seed, drug)

# adtype = Optimization.AutoForwardDiff()
# opt = Optim.LBFGS()
t_iter=100
i=7
patient = patients[findfirst(x->x.idx==i, patients)]
ode_p = patient.ode_parameters
scaler = patient.scaling
# init_doses = ones(length(drug_dose_times)).*minimum([patient.output_measures["0.99999effect"].doses[1], max_tested])
# init_doses = ones(length(drug_dose_times)).*maximum([patient.output_measures["0.99999effect"].doses[1], max_tested])
init_logits = ones(length(dose_mults), num_dose_times)
# ode_params = [population[:,i];scaling[i]]

#some preliminary checks
loss(init_logits, ode_p, scaler)
plotter_gif(init_logits, 0.0)
Random.seed!(123)
plotter_gif(rand(length(dose_mults), num_dose_times), 0.0)

# gradient descent optimization of the dose amounts
#     # Redefine the optimization function for the current patient
# optf = Optimization.OptimizationFunction((x, _) -> loss(x, ode_p, scaler), adtype)
# # lower_bound = ones(length(drug_dose_times)).*0.0#min_tested
# # upper_bound = ones(length(drug_dose_times)).*patient.output_measures["0.99999effect"].doses[1]#(max_tested*1.15)#
# optprob = Optimization.OptimizationProblem(optf, init_logits)#, lb=lower_bound, ub=upper_bound)
# callback, losses, iter_ref, anim = create_callback(t_iter,verbose=true, animate=true, progress_bar=false, saving=false, early_stopping=true)
# res = Optimization.solve(optprob, opt, maxiters=t_iter, callback=callback)
# gif(anim, "results/$(drug)_probmask_test_lbfgs.gif", fps=10)

function prob_to_dose(prob)
    final_pmf = exp.(prob ./ temp)
    final_pmf ./= sum(final_pmf, dims=1)
    # println(argmax(final_pmf, dims=1))
    doses = (repeat(dose_mults,1,num_dose_times)[argmax(final_pmf, dims=1)]).*(max_tested/mult_num)
    # dose = ((prob' * dose_mults) .* dose_unit)'
end


lr = 0.1
logits = deepcopy(init_logits)

anim = Animation()
loss_values = losses = []

early_stopping=true
window = 5

lr_schedule = false
decay_rate = 0.1

adam = true
beta1 = 0.9
beta2 = 0.999
epsilon = 1e-8
m = zeros(size(logits))
v = zeros(size(logits))

for i in 1:t_iter
    Random.seed!(123)
    current_loss = loss(logits, ode_p, scaler)
    push!(losses, current_loss)
    @info "Epoch: $i, Loss: $current_loss, LR: $lr"
    # println(argmax(logits, dims=1))
    println(prob_to_dose(logits))
    frame(anim, plotter_gif(logits, current_loss))

    if early_stopping
        if length(loss_values) > window+1
            avg_change = mean([abs(loss_values[end-i] - loss_values[end-(i+1)]) for i in 0:window])
            if avg_change < 1e-5
                println("Early stopping triggered!")
                break
            end
        end
    end

    Random.seed!(123)
    grad = ForwardDiff.gradient(λ -> loss(λ, ode_p, scaler), logits)
    if adam
        m = beta1 * m + (1 - beta1) * grad
        v = beta2 * v + (1 - beta2) * grad.^2
    
        m_hat = m / (1 - beta1^i)
        v_hat = v / (1 - beta2^i)
    
        logits -= lr * m_hat ./ (sqrt.(v_hat) .+ epsilon)
    else
        logits -= lr * grad
    end

    if lr_schedule
        lr /= (1 + decay_rate * i)
    end
end
gif(anim, "results/$(drug)_probmask_test_ADAM_temp$(temp)_lr$(lr).gif", fps=10)
plot(losses, label="Loss", xlabel="Epoch", ylabel="Loss", title="Loss vs Epoch")
savefig("results/$(drug)_probmask_losses_test_ADAM_temp$(temp)_lr$(lr).png")