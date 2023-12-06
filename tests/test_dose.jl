using Pkg
Pkg.activate(".")
Pkg.instantiate()

#holds the NLMEM parameters for the PK model
include("../assets/pk_params.jl")

#inlcude the functions that are used in setting up the struct
include("../scripts/setup/generate_vp.jl")
include("../scripts/setup/compute_dose.jl")

function simulate_to_check(amount_per_mult_per_patient, drug ,mults=[0.01, 0.1, 0.5, 1.0, 1.5, 2.0])
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
        return abs2(peak_PlaRG - mult*(ic50*(volume*1000))) #trajectories don't cross so we can be sure that the amount is minimal for reaching this peak
    end
    residuals = zeros(length(mults), size(amount_per_mult_per_patient, 2))
    for (j,mult) in enumerate(mults)
        for (i,dose) in enumerate(amount_per_mult_per_patient[j,:])
            println("\rPatient: $i")
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
            residuals[j,i] = output(input_dose[1], pk_param_values[:,i],ic50[i],volume[i],mult)
            # end
        end
    end
    return residuals
end

function test_dose_matrix(num_patients, seed, drug)

    #loads the model equations, dosing event functions and the parameter values
    include("../model/$(drug)_pkpd2.jl")
    include("../model/$(drug)_params.jl")

    drug_params = eval(Symbol(uppercase(drug) * "_params"))
    # Generate population
    population = generate_virtual_population(TMZ_params, drug_params, ode_params, param_order, num_patients, seed)
    @info "Computing dose matrix..."
    ic50_fractional_doses = compute_dose(population, drug_params)
    residuals = simulate_to_check(ic50_fractional_doses, drug)
    # check that all residuals are sufficiently low (i.e. the dose is sufficient to reach the peak) otherwise print Boundary

    if !all(residuals .< 1e-10)
        @info "haven't reached the peak"
    end
    
    return ic50_fractional_doses
end

doses = test_dose_matrix(1000, 123, "rg")

using Plots

#make a heat map of the doses
heatmap(doses')