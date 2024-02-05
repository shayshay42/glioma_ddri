using Pkg
Pkg.activate(".")
Pkg.instantiate()

#holds the NLMEM parameters for the PK model
include("../assets/pk_params.jl")
include("../utilities/utils.jl")

#inlcude the functions that are used in setting up the struct
include("../scripts/setup/generate_vp.jl")
include("../scripts/setup/compute_dose.jl")
include("../scripts/setup/compute_dose_2.jl")

function simulate_to_check(amount_per_mult_per_patient, population, drug, pk_params, mults=[0.01, 0.1, 0.5, 1.0, 1.5, 2.0])
    function dose_affect!(integrator)
        # SciMLBase.set_proposed_dt!(integrator, 0.01)
        integrator.u[1] += integrator.p[end]
    end

    event_times = collect(1:18).*24.0
    cb = PresetTimeCallback(event_times, dose_affect!)
    
    tspan = (0.0, event_times[end]+(10.0*24.0))

    p = rand(6)
    u0 = [0.0, 0.0, 0.0]
    # tspan = (0.0, 18.0).*24.0
    prob = ODEProblem(drug_pk!, u0, tspan, p)
    function output(dose, pk_values, ic50, volume, mult)
        tmp_prob = remake(prob, p=[pk_values...,dose])
        tmp_sol = solve(tmp_prob, callback=cb)# abstol=1e-10, reltol=1e-10,dtmax=1)#, alg_hints=[:stiff])
        peak_PlaRG = maximum(tmp_sol[2, :])  # Assuming PlaRG is the second variable
        return abs(peak_PlaRG - mult*(ic50*(volume*1000))) #trajectories don't cross so we can be sure that the amount is minimal for reaching this peak
    end
    residuals = zeros(length(mults), size(amount_per_mult_per_patient, 2))
    for (j,mult) in enumerate(mults)
        for (i,dose) in enumerate(amount_per_mult_per_patient[j,:])
            # println("\rPatient: $i")
            pk_param_names = collect(keys(pk_params))
            pk_indices = indexin(pk_param_names, param_order)
            pk_param_values= population[pk_indices, :]
            ic50 = population[indexin(["IC50_2"], param_order), :]
            if "Vpla" in param_order
                volume = population[indexin(["Vpla"], param_order), :]
            elseif "V2" in param_order
                volume = population[indexin(["V2"], param_order), :]
            end
            # Optimization for each multiplier
            # for (j,mult) in enumerate(mults)
            # result = optimize(x -> objective(x, pk_param_values[:,i],ic50[i],volume[i],mult), [initial_guess[j]], LBFGS())
            # input_dose = Optim.minimizer(result)
            # patient_doses[j, i] = input_dose[1]
            residuals[j,i] = output(dose, pk_param_values[:,i],ic50[i],volume[i],mult)
            # end
        end
    end
    return residuals
end

"""
`test_dose_matrix(num_patients, seed, drug, method)`

Compute the dose matrix for a given number of patients using a specified method.

# Arguments
- `num_patients`: Number of patients to simulate.
- `seed`: Random seed for population generation.
- `drug`: The drug for which the dose matrix is being computed.
- `method`: The method used for computation, e.g., 'Optimization' or 'BVP'.

# Returns
- Dose matrix and, if applicable, minimal objectives.
"""
function test_dose_matrix(num_patients, seed, drug, method)
    # Constants for methods
    OPTIMIZATION = "optimization"
    BVP = "bvp"

    # Load the model and parameters
    include("../model/$(drug)_pkpd2.jl")
    include("../model/$(drug)_params.jl")

    drug_params = eval(Symbol(uppercase(drug) * "_params"))

    # Generate population
    population = generate_virtual_population(TMZ_params, drug_params, ode_params, param_order, num_patients, seed)
    @info "Computing dose matrix via $method..."

    # Compute dose matrix
    ic50_fractional_doses = nothing
    minimal_objectives = nothing
    if method == OPTIMIZATION
        ic50_fractional_doses, minimal_objectives = compute_dose(population, drug_params)
    elseif method == BVP
        ic50_fractional_doses = compute_dose2(population, drug_params)
    else
        error("Unknown method: $method")
    end

    @info "Done computing dose matrix. Performing a quick test..."
    residuals = simulate_to_check(ic50_fractional_doses, population, drug, drug_params)

    # Check residuals
    if !all(residuals .< 1e-8)
        @info "Haven't reached the peak with $method"
    end
    
    # Package results in a NamedTuple for consistency
    results = (doses = ic50_fractional_doses, residuals = residuals)
    if method == OPTIMIZATION
        results = (doses = ic50_fractional_doses, residuals = residuals, objectives = minimal_objectives)
    end
    return results
    # return method == OPTIMIZATION ? (ic50_fractional_doses, residuals, minimal_objectives) : ic50_fractional_doses, residuals
end

optim_results_bfgs = test_dose_matrix(5, 123, "rg", "optimization")
optim_results.doses
sqrt.(optim_results.objectives)
optim_results.residuals

bvp_results = test_dose_matrix(5, 123, "rg", "bvp")
bvp_results.doses
bvp_results.residuals

#check the difference between the optim and the bvp
# doses_diff = doses_optim .- doses_bvp
diff = optim_results.doses .- bvp_results.doses
#evaluate if bellow a 1e-8
abs.(diff) .< 1e-8

# using Plots

# #make a heat map of the doses
# heatmap(doses')


using Printf

function compare_matrices(optim_results, bvp_results, threshold=1e-8, precision=4)
    # Extract doses and residuals from both methods
    doses_optim = optim_results.doses
    residuals_optim = optim_results.residuals
    doses_bvp = bvp_results.doses
    residuals_bvp = bvp_results.residuals

    # Compute differences and check against the threshold
    diff = abs.(doses_optim - doses_bvp)
    below_threshold = diff .< threshold

    function print_matrix(mat, title, precision=4)
        println(title)
        for i in 1:size(mat, 1)
            for j in 1:size(mat, 2)
                # Directly format each element
                formatted_value = round(mat[i, j], digits=precision)
                @printf("%8.4f ", formatted_value)
            end
            println()
        end
    end
    
    # Print matrices
    println("Comparing matrices from 'optimization' and 'bvp' methods:")
    print_matrix(doses_optim, "Doses - Optimization Method:")
    print_matrix(doses_bvp, "Doses - BVP Method:")
    print_matrix(residuals_optim, "Residuals - Optimization Method:")
    print_matrix(residuals_bvp, "Residuals - BVP Method:")
    print_matrix(diff, "Absolute Differences in Doses:")
    print_matrix(below_threshold, "Differences Below Threshold ($threshold):")
end
# Example usage
optim_results = test_dose_matrix(5, 123, "rg", "optimization")
bvp_results = test_dose_matrix(5, 123, "rg", "bvp")

compare_matrices(optim_results, bvp_results)

compare_matrices(optim_results_bfgs, bvp_results)



function simulate_to_plot(amount_per_mult_per_patient, population, drug, pk_params, mults=[0.01, 0.1, 0.5, 1.0, 1.5, 2.0])
    function dose_affect!(integrator)
        # SciMLBase.set_proposed_dt!(integrator, 0.01)
        integrator.u[1] += integrator.p[end]
    end

    event_times = collect(1:18).*24.0
    cb = PresetTimeCallback(event_times, dose_affect!)
    
    tspan = (0.0, event_times[end]+(10.0*24.0))
    timepoints = collect(0:0.1:tspan[2])

    p = rand(6)
    u0 = [0.0, 0.0, 0.0]
    # tspan = (0.0, 18.0).*24.0
    prob = ODEProblem(drug_pk!, u0, tspan, p)
    function output(dose, pk_values, ic50, volume, mult)
        tmp_prob = remake(prob, p=[pk_values...,dose])
        tmp_sol = solve(tmp_prob, callback=cb, saveat=timepoints)# abstol=1e-10, reltol=1e-10,dtmax=1)#, alg_hints=[:stiff])
        # peak_PlaRG = maximum(tmp_sol[2, :])  # Assuming PlaRG is the second variable
        return tmp_sol #trajectories don't cross so we can be sure that the amount is minimal for reaching this peak
    end
    sols = zeros(length(mults), size(amount_per_mult_per_patient, 2), length(timepoints))
    for (j,mult) in enumerate(mults)
        for (i,dose) in enumerate(amount_per_mult_per_patient[j,:])
            # println("\rPatient: $i")
            pk_param_names = collect(keys(pk_params))
            pk_indices = indexin(pk_param_names, param_order)
            pk_param_values= population[pk_indices, :]
            ic50 = population[indexin(["IC50_2"], param_order), :]
            if "Vpla" in param_order
                volume = population[indexin(["Vpla"], param_order), :]
            elseif "V2" in param_order
                volume = population[indexin(["V2"], param_order), :]
            end
            # Optimization for each multiplier
            # for (j,mult) in enumerate(mults)
            # result = optimize(x -> objective(x, pk_param_values[:,i],ic50[i],volume[i],mult), [initial_guess[j]], LBFGS())
            # input_dose = Optim.minimizer(result)
            # patient_doses[j, i] = input_dose[1]
            sol = output(dose, pk_param_values[:,i],ic50[i],volume[i],mult)
            sols[j,i,:] = sol[2,:]
            # end
        end
    end
    return timepoints, sols
end

times, sols = simulate_to_plot(optim_results.doses, population, "rg", TMZ_params)

using Plots
data_matrix = sols[1,1,:]
timepoints = sols ./24

title = "Tumor Trajectories of Virtual Population \n Treated with TMZ and GDC max"
means = vec(mean(data_matrix, dims=1))

# Plot the trajectories with a label for the legend
plot(timepoints, data_matrix'[:,1]     , alpha=0.3, linecolor=:lightgray, grid=false, label="Tumor Trajectories", legend=:topleft)
plot!(timepoints, data_matrix'[:,2:end], alpha=0.3, linecolor=:lightgray, grid=false, label=false)
# Plot the mean and standard deviation with labels for the legend
plot!(timepoints, means, 
    ribbon=(vec(std(data_matrix, dims=1)), vec(std(data_matrix, dims=1))), 
    fillalpha=0.5, linewidth=3, linecolor=:pink, fillcolor=:lightblue, 
    label="Mean ± SD")

# Set the title and axis labels
title!(title)
xlabel!("Time (days)")
ylabel!("Tumor Volume (mL)")

savefig("./results/$(drug)_model_TMZonly.png")
savefig("./results/$(drug)_model_TMZonly.svg")
