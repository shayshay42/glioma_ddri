using Optim, Optimization, ForwardDiff, DifferentialEquations, StaticArrays, OptimizationOptimisers, OptimizationOptimJL
using ComponentArrays

function compute_dose_optim(population, gradation)

    function dose_affect!(integrator)
        # SciMLBase.set_proposed_dt!(integrator, 0.01)
        integrator.u[1] += integrator.p.dose
    end
    # cb = PeriodicCallback(dose_affect!, 24.0, initial_affect=true)
    event_times = collect(1:18).*24.0
    cb = PresetTimeCallback(event_times, dose_affect!)

    tspan = (0.0, event_times[end]+(10.0*24.0))

    p = ComponentArray(θ=ode_params, dose=1000.0)
    u0 = [0.0, 0.0, 0.0]
    # tspan = (0.0, 18.0).*24.0
    function drug_pk!(du, u, p, t)
        gamma_1,psi,C0,D0,r,K,BW,IC50_1,Imax_1,IC50_2,gamma_2,Imax_2,xi,VD1,Cl1,k23,ka1,k32,Cl2,ka2,Vpla,Q,Vtis = p.θ
        AbsRG, PlaRG, TisRG = u
        dAbsRG = -ka2 * AbsRG
        dPlaRG = ka2 * AbsRG - (Cl2 / Vpla) * PlaRG + (Q / Vtis) * TisRG - (Q / Vpla) * PlaRG
        dTisRG = -(Q / Vtis) * TisRG + (Q / Vpla) * PlaRG
        du .= [dAbsRG, dPlaRG, dTisRG]
    end
    prob = ODEProblem(drug_pk!, u0, tspan, p)
    sys = modelingtoolkitize(prob)
    sys = structural_simplify(sys)
    func = ODEFunction(sys, jac=true)
    prob = ODEProblem(func, u0, tspan, p)

    function objective(dose, pk_values, ic502, volume, gamma2, imax2, psi, mult)
        p=ComponentArray(θ=pk_values, dose=dose[1])
        tmp_prob = remake(prob, p=p)
        tmp_sol = solve(tmp_prob, callback=cb)# abstol=1e-10, reltol=1e-10,dtmax=1)#, alg_hints=[:stiff])
        # peak_PlaRG = maximum(tmp_sol[2, :])  # Assuming PlaRG is the second variable
        cPlaRG_sol = tmp_sol[2,:] ./ (1000 * volume)
        cPlaRG_sol = erelu.(cPlaRG_sol)
        exp2_sol = (cPlaRG_sol ./(psi*ic502)).^gamma2
        E_sol = (exp2_sol.*imax2) ./ (exp2_sol .+ 1)
        return abs2(maximum(E_sol) - mult) #trajectories don't cross so we can be sure that the amount is minimal for reaching this peak
    end

    mults = gradation
    avg_human_surface_area=1.7
    # dosage guess for each multiplier of IC50
    initial_guess = gradation.* (10000*avg_human_surface_area)
    # initialize array to hold patient doses for each multiplier
    patient_doses = zeros(length(mults), size(population, 2))
    minimal_doses = zeros(length(mults), size(population, 2))
    Threads.@threads for i in 1:size(population, 2)
        println("\rPatient: $i")
        # pk_param_names = collect(keys(pk_params))
        # pk_indices = indexin(pk_param_names, param_order)
        # pk_param_values= population[pk_indices, :]
        ic502 = population[indexin(["IC50_2"], param_order), i][1]
        gamma2 = population[indexin(["gamma_2"], param_order), i][1]
        imax2 = population[indexin(["Imax_2"], param_order), i][1]
        psi=1
        if "Vpla" in param_order
            volume = population[indexin(["Vpla"], param_order), i][1]
        elseif "V2" in param_order
            volume = population[indexin(["V2"], param_order), i][1]
        end
        # Optimization for each multiplier
        opts = Optim.Options(g_tol = 1e-6, f_tol = 1e-6, iterations = 1000)
        for (j,mult) in enumerate(mults)
            initial_value = [initial_guess[j]]
            # result = optimize(x -> objective(x, population[:, i], ic502, volume, gamma2, imax2, psi, mult), [initial_guess[j]], LBFGS(), opts)
            result = optimize(x -> objective(x, population[:, i], ic502, volume, gamma2, imax2, psi, mult), initial_value, opts)
            
            patient_doses[j, i] = result.minimizer[1]
            minimal_doses[j, i] = result.minimum
        end
    end
    result = (doses_matrix=patient_doses, mins=minimal_doses )
    return result
end


#TESTING


# function dose_affect!(integrator)
#     # SciMLBase.set_proposed_dt!(integrator, 0.01)
#     integrator.u[1] += integrator.p.dose
# end
# # cb = PeriodicCallback(dose_affect!, 24.0, initial_affect=true)
# event_times = collect(1:18).*24.0
# cb = PresetTimeCallback(event_times, dose_affect!)

# tspan = (0.0, event_times[end]+(10.0*24.0))

# p = ComponentArray(θ=ode_params, dose=1000.0)
# u0 = [0.0, 0.0, 0.0]
# # tspan = (0.0, 18.0).*24.0
# function drug_pk!(du, u, p, t)
#     gamma_1,psi,C0,D0,r,K,BW,IC50_1,Imax_1,IC50_2,gamma_2,Imax_2,xi,VD1,Cl1,k23,ka1,k32,Cl2,ka2,Vpla,Q,Vtis = p.θ
#     AbsRG, PlaRG, TisRG = u
#     dAbsRG = -ka2 * AbsRG
#     dPlaRG = ka2 * AbsRG - (Cl2 / Vpla) * PlaRG + (Q / Vtis) * TisRG - (Q / Vpla) * PlaRG
#     dTisRG = -(Q / Vtis) * TisRG + (Q / Vpla) * PlaRG
#     du .= [dAbsRG, dPlaRG, dTisRG]
# end
# prob = ODEProblem(drug_pk!, u0, tspan, p)
# sys = modelingtoolkitize(prob)
# sys = structural_simplify(sys)
# func = ODEFunction(sys, jac=true)
# prob = ODEProblem(func, u0, tspan, p)

# function objective(dose, pk_values, ic502, volume, gamma2, imax2, psi, mult)
#     p=ComponentArray(θ=pk_values, dose=dose[1])
#     tmp_prob = remake(prob, p=p)
#     tmp_sol = solve(tmp_prob, callback=cb)# abstol=1e-10, reltol=1e-10,dtmax=1)#, alg_hints=[:stiff])
#     # peak_PlaRG = maximum(tmp_sol[2, :])  # Assuming PlaRG is the second variable
#     cPlaRG_sol = tmp_sol[2,:] ./ (1000 * volume)
#     # cPlaRG_sol = erelu.(cPlaRG_sol)
#     exp2_sol = (cPlaRG_sol ./(psi*ic502)).^gamma2
#     E_sol = (exp2_sol.*imax2) ./ (exp2_sol .+ 1)
#     return abs2(maximum(E_sol) - mult) #trajectories don't cross so we can be sure that the amount is minimal for reaching this peak
# end

# mults = [1e-5, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99]
# avg_human_surface_area=1.7
# # dosage guess for each multiplier of IC50
# initial_guess = [20.0, 100.0, 500.0, 1000.0, 2000.0, 3000.0].*avg_human_surface_area # Adjust this as needed
# # initialize array to hold patient doses for each multiplier

# population  = vp_param_matrix
# patient_doses = zeros(length(mults), size(population, 2))
# minimal_doses = zeros(length(mults), size(population, 2))
# i = 1
# println("\rPatient: $i")
# # pk_param_names = collect(keys(pk_params))
# # pk_indices = indexin(pk_param_names, param_order)
# # pk_param_values= population[pk_indices, :]
# ic502 = population[indexin(["IC50_2"], param_order), i][1]
# gamma2 = population[indexin(["gamma_2"], param_order), i][1]
# imax2 = population[indexin(["Imax_2"], param_order), i][1]
# psi=1
# if "Vpla" in param_order
#     volume = population[indexin(["Vpla"], param_order), i][1]
# elseif "V2" in param_order
#     volume = population[indexin(["V2"], param_order), i][1]
# end
# # Optimization for each multiplier
# j=1
# mult = 0.9

# objective(initial_guess[j], population[:, i], ic502, volume, gamma2, imax2, psi, mult)
# using ForwardDiff
# forwarddiff_gradient_ρ = ForwardDiff.gradient(x -> objective(x, population[:, i], ic502, volume, gamma2, imax2, psi, mult), [initial_guess[j]])
# # do lbfgs descent
# # Initial guess for the optimization

# # Setting up and running the optimization
# opts = Optim.Options(g_tol = 1e-6, f_tol = 1e-6, iterations = 1000)
# function optim_objective(x, grad)
#     if length(grad) > 0
#         g = ForwardDiff.gradient(x -> objective(x, population[:, i], ic502, volume, gamma2, imax2, psi, mult), x)
#         grad .= g
#     end
#     return objective(x, population[:, i], ic502, volume, gamma2, imax2, psi, mult)
# end
# # result = optimize(x -> objective(x, population[:, i], ic502, volume, gamma2, imax2, psi, mult), [initial_guess[j]], LBFGS(), opts)
# result = optimize(x -> objective(x, population[:, i], ic502, volume, gamma2, imax2, psi, mult), [1000.0], opts)
# result.res

# # Extracting the optimized dose
# optimal_dose = Optim.minimizer(result)
# println("Optimal dose: ", optimal_dose)

# # You may also want to check the optimization status and value
# println("Optimization status: ", Optim.converged(result))
# println("Minimum found: ", Optim.minimum(result))


# initial_value = [initial_guess[j]]
# lower_bound = [0.0]
# upper_bound = [Inf]

# optf = Optimization.OptimizationFunction((dose, p) -> objective(dose, p, ic502, volume, gamma2, imax2, psi, mult), Optimization.AutoForwardDiff())
# optprob = OptimizationProblem(optf, initial_value, (population[:, i], ic502[i], volume[i], gamma2[i], imax2[i], psi, mult), lb=lower_bound, ub=upper_bound)
# optimizer = Optim.LBFGS()
# result = solve(optprob, optimizer, maxiters=10)

# patient_doses[j, i] = result.minimizer[1]
# minimal_doses[j, i] = result.minimum
