include("../../../src/setup.jl")
#loads the model equations, dosing event functions and the parameter values
drug = "rg"

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

const mult_num = 20
const dose_mults = [0:mult_num...]
const unit_dose = single_max/mult_num
const possible_doses = [0:unit_dose:single_max...]
const temp = 1

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
end

num_patients = 10
seed = 123

patients = generate_patients_struct(num_patients, seed, drug)

drug_dose_times = eval(Symbol(drug * "_dosetimes"))
#order patients vector by the .idx property of the elements
sort!(patients, by=x->x.idx)
@info "Starting optimization for $(uppercase(drug))-model ..."
optima = Vector{optimum_result}(undef, length(patients))
patients_optim = []

t_iter=100

# Initialize an Atomic counter
# completed_patients = Threads.Atomic{Int64}(0)
Random.seed!(seed)


i=1
println("Processing patient $i, Completed patients: $(completed_patients[])")
#retrive element in patients vector with .idx property equal to i
patient = patients[findfirst(x->x.idx==i, patients)]
ode_p = patient.ode_parameters
scaler = patient.scaling
# init_doses = ones(length(drug_dose_times)).*minimum([patient.output_measures["0.99999effect"].doses[1], max_tested])
# init_doses = ones(length(drug_dose_times)).*maximum([patient.output_measures["0.99999effect"].doses[1], max_tested])
init_logits = ones(mult_num+1, num_dose_times)


anim=Animation()
loss_values = losses = []


callback, losses, iter_ref = create_callback(t_iter,verbose=false, animate=true, progress_bar=false, saving=false, early_stopping=true)
logits = deepcopy(init_logits)

lr=0.1


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


gif(anim, "results/$(drug)_probmask_test_ADAM_temp$(temp)_lr$(lr)_$(today_date).gif", fps=10)
plot(losses, label="Loss", xlabel="Epoch", ylabel="Loss", title="Loss vs Epoch")
savefig("results/$(drug)_probmask_losses_test_ADAM_temp$(temp)_lr$(lr)_$(today_date).png")