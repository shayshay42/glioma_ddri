using DifferentialEquations, Optim, OptimizationOptimisers

function compute_dose(population, pk_params)

    function dose_affect!(integrator)
        # SciMLBase.set_proposed_dt!(integrator, 0.01)
        integrator.u[1] += integrator.p[end]
    end
    # cb = PeriodicCallback(dose_affect!, 24.0, initial_affect=true)
    event_times = collect(1:18).*24.0
    cb = PresetTimeCallback(event_times, dose_affect!)
    
    tspan = (0.0, event_times[end]+(10.0*24.0))

    p = rand(6)
    u0 = [0.0, 0.0, 0.0]
    # tspan = (0.0, 18.0).*24.0
    prob = ODEProblem(drug_pk!, u0, tspan, p)

    function objective(dose, pk_values, ic50, volume, mult)
        tmp_prob = remake(prob, p=[pk_values...,dose[1]])
        tmp_sol = solve(tmp_prob, callback=cb)# abstol=1e-10, reltol=1e-10,dtmax=1)#, alg_hints=[:stiff])
        peak_PlaRG = maximum(tmp_sol[2, :])  # Assuming PlaRG is the second variable
        return abs2(peak_PlaRG - mult*(ic50*(volume*1000))) #trajectories don't cross so we can be sure that the amount is minimal for reaching this peak
    end

    mults = [0.01, 0.1, 0.5, 1.0, 1.5, 2.0]
    avg_human_surface_area=1.7
    # dosage guess for each multiplier of IC50
    initial_guess = [20.0, 100.0, 500.0, 1000.0, 2000.0, 3000.0].*avg_human_surface_area # Adjust this as needed
    # initialize array to hold patient doses for each multiplier
    patient_doses = zeros(length(mults), size(population, 2))
    minimal_doses = zeros(length(mults), size(population, 2))
    Threads.@threads for i in 1:size(population, 2)
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
        for (j,mult) in enumerate(mults)
            # Assuming initial_guess[j] is a finite value within the bounds
            initial_value = [initial_guess[j]]  # Replace this with your actual initial guess
            # Lower and upper bounds
            lower_bound = [0.0]  # Lower bound
            upper_bound = [Inf]  # Upper bound
            # Optimizer
            optimizer = Fminbox(BFGS())
            # Optimization call
            result = optimize(x -> objective(x, pk_param_values[:,i], ic50[i], volume[i], mult), 
                            lower_bound, upper_bound, initial_value, optimizer)#, Optim.Options(iterations=10))

            input_dose = Optim.minimizer(result)
            patient_doses[j, i] = input_dose[1]
            minimal_doses[j, i] = Optim.minimum(result)
            # if j == 6
            #     println("Patient $i: $result")
            # end
        end
    end
    return patient_doses, minimal_doses
end