function pk_pd!(du,u,p,t)#(du::Vector{Float64}, u::Vector{Float64}, p::Vector{Float64}, t::Float64)
    # gamma_1, psi, C0, D0, r, K, BW, IC50_1, Imax_1, IC50_2, gamma_2, Imax_2, xi, VD1, Cl1, k23, ka1, k32, Cl2, ka2, Vpla, Q, Vtis = p[1:length(ode_params)]
    gamma_1,psi,C0,D0,r,K,BW,IC50_1,Imax_1,IC50_2,gamma_2,Imax_2,xi,VD1,Cl1,k23,ka1,k32,Cl2,ka2,Vpla,Q,Vtis = p[1:p_num]#length(ode_params)]
    # for name in keys(p)
    #     eval(:($name = nt.$name))
    # end
    xi = IC50_1/IC50_2
    C, D, AbsTMZ, PlaTMZ, CSFTMZ, AbsRG, PlaRG, TisRG, pRGauc, cAUC = u

    dAbsTMZ = -ka1 * AbsTMZ
    dPlaTMZ = ka1 * AbsTMZ - (Cl1 / VD1) * PlaTMZ - k23 * PlaTMZ + k32 * CSFTMZ
    dCSFTMZ = k23 * PlaTMZ - k32 * CSFTMZ

    dAbsRG = -ka2 * AbsRG
    dPlaRG = ka2 * AbsRG - (Cl2 / Vpla) * PlaRG + (Q / Vtis) * TisRG - (Q / Vpla) * PlaRG
    dTisRG = -(Q / Vtis) * TisRG + (Q / Vpla) * PlaRG
   
    cCSFTMZ = CSFTMZ / 140
    cPlaRG = PlaRG / (1000 * Vpla)

    #domain error
    cCSFTMZ = erelu(cCSFTMZ)
    cPlaRG = erelu(cPlaRG)
    C = erelu(C)

    # pi1 = psi * IC50_1
    exp1 = (cCSFTMZ /(psi*IC50_1))^gamma_1
    exp2 = (cPlaRG  /(psi*IC50_2))^gamma_2
    exp12 = exp1 * exp2

    Imax_sum = Imax_1 + Imax_2
    E_num = Imax_1 * exp1 + Imax_2 * exp2 + Imax_sum * exp12 - Imax_1 * Imax_2 * exp12
    E_den = exp1 + exp2 + exp12 + 1
    E = E_num / E_den

    t = (-log(log(C/K) / log(C0/K)) / r) + 72
    fun = K * (C0/K) ^ exp(-r * t)
    # delta = (E * fun) / (72 * C)

    dC = C * r * log(K / C) - (E * fun / 72)#delta * C
    dD = (E * fun / 72)#delta * C

    dpRGauc = PlaRG
    dcAUC = C
    
    du .= [dC, dD, dAbsTMZ, dPlaTMZ, dCSFTMZ, dAbsRG, dPlaRG, dTisRG, dpRGauc, dcAUC]
end
# create a dictionary of the state and their index
states = OrderedDict(zip(["C", "D", "AbsTMZ", "PlaTMZ", "CSFTMZ", "AbsRG", "PlaRG", "TisRG", "PlaRGAUC", "cAUC"], 1:10))

function generate_perturbation_amounts(state::Int, state_dim::Int, amount::Vector{Float64})
    perturbation_amounts = zeros(state_dim, length(amount))
    perturbation_amounts[state, :] = amount
    return perturbation_amounts
end

function generate_perturbation_data(states::Int, perturb_times::Vector{Vector{Float64}}, perturbation_amounts::Vector{Matrix{Float64}})
    all_time_points = vcat(perturb_times...)
    unique_time_points = sort(unique(all_time_points))
    perturbation_matrix = zeros(states, length(unique_time_points))
    for (i, time_points) in enumerate(perturb_times)
        for (j, time) in enumerate(time_points)
            index = findfirst(isequal(time), unique_time_points)
            perturbation_matrix[:, index] += perturbation_amounts[i][:, j]
        end
    end
    return unique_time_points, perturbation_matrix
end

include("../../../utilities/utils.jl")
include("../../../model/rg_params.jl")

#time frames
hours = 24.0
end_time = 28.0*5.0
end_treat = 42.0

#dose amounts
avg_human_surface_area = 1.7 #m^2
tmz_treat_dose = 75.0*avg_human_surface_area
tmz_adjuv_dose = 150.0*avg_human_surface_area
dose_amount = 1920.0*avg_human_surface_area #1800
const max_tested = dose_amount
const single_max = max_tested*1.15

tmz_treat_dosetimes = spaced_list(end_treat,1.0,0.0,0.0).*hours
tmz_adjuv_dosetimes = spaced_list(end_time,5.0,23.0,end_treat+28.0).*hours
rg_dosetimes = spaced_list(end_time-1.0,18.0,10.0,0.0).*hours #its the opposite18 days off and 10days on
const num_dose_times = length(rg_dosetimes)

doses = ones(length(rg_dosetimes)).*dose_amount

inject_times = sort(unique([rg_dosetimes;tmz_treat_dosetimes;tmz_adjuv_dosetimes]));

function adjust_dose(x)
    return min(dose_amount,relu(x))
end

# Convert arrays to sets for O(1) lookup
tmz_treat_dosetimes_set = Set(tmz_treat_dosetimes)
tmz_adjuv_dosetimes_set = Set(tmz_adjuv_dosetimes)
rg_dosetimes_set = Set(rg_dosetimes)

perturbation_amounts = [generate_perturbation_amounts(states["AbsTMZ"], length(states), ones(length(tmz_treat_dosetimes)).*tmz_treat_dose), generate_perturbation_amounts(states["AbsTMZ"], length(states), ones(length(tmz_adjuv_dosetimes)).*tmz_adjuv_dose), generate_perturbation_amounts(states["AbsRG"], length(states), doses)]
perturb_times = [tmz_treat_dosetimes, tmz_adjuv_dosetimes, rg_dosetimes]
perturb_times, perturb_matrix = generate_perturbation_data(length(states), perturb_times, perturbation_amounts)

# Create a ComponentArray for parameters
params = ComponentArray(θ=ode_params, perturb_times=perturb_times, perturb_matrix=perturb_matrix)

using DifferentialEquations, Plots, LinearAlgebra, Zygote, SciMLSensitivity, ForwardDiff, FiniteDiff, Flux

function affect_dose!(integrator)
    index = findfirst(isequal(integrator.t), integrator.p.perturb_times)
    integrator.u += integrator.p.perturb_matrix[:, index]
end

hit = PresetTimeCallback(inject_times, affect_dose!)

using ModelingToolkit

u0 = zeros(10)
u0[1] = C0
tspan = (0,end_time+7).*hours
# p = [ode_params..., default_scaling..., doses...]
# p = [ode_params..., doses...]
# p = [ode_params; default_scaling; doses]
prob = ODEProblem(pk_pd!,u0,tspan,params)
sys = modelingtoolkitize(prob)
sys = structural_simplify(sys)
func = ODEFunction(sys, jac=true)
prob = ODEProblem(func, u0, tspan, params)

solution = solve(prob, Rodas5P(), callback=hit)

plot(solution)

function loss(ρ, params)
    params.perturb_matrix = ρ
    prob_temp = remake(prob, p=params)
    solve(prob_temp, Rodas5P(), callback=hit)[states["cAUC"], end]
end

function loss(ρ)
    prob_temp = remake(prob, p=ρ)
    solve(prob_temp, Rodas5P(), callback=hit)[end, end]
end

loss(params)

# Example of gradient computation with respect to μ
zygote_gradient_ρ = Zygote.gradient(x -> loss(x, params), params.perturb_matrix)

forwarddiff_gradient_ρ = ForwardDiff.gradient(loss, params)

finitediff_gradient_ρ = FiniteDiff.finite_difference_gradient(x -> loss(x, params), params.perturb_matrix)
