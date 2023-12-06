#want to instead of the whole set just get 1 dose to reflec this
using DifferentialEquations, Opti

function compute_dose(population, pk_params)
    
    p = rand(6)
    u0 = [0.0, 0.0, 0.0]
    tspan = (0.0, 3.0).*24.0
    prob = ODEProblem(drug_pk!, u0, tspan, p)

    function objective(dose, pk_values, ic50, volume, mult)
        tmp_prob = remake(prob, u0=[dose, 0.0, 0.0], p=pk_values)
        tmp_sol = solve(tmp_prob)
        peak_PlaRG = maximum(tmp_sol[2, :])  # Assuming PlaRG is the second variable
        return abs2(peak_PlaRG - mult*(ic50*(volume*1000))) #trajectories don't cross so we can be sure that the amount is minimal for reaching this peak
    end

    mults = [0.01, 0.1, 0.5, 1.0, 1.5, 2.0]
    avg_human_surface_area=1.7
    # dosage guess for each multiplier of IC50
    initial_guess = [20.0, 100.0, 500.0, 1000.0, 2000.0, 3000.0].*avg_human_surface_area # Adjust this as needed
    # initialize array to hold patient doses for each multiplier
    patient_doses = zeros(length(mults), size(population, 2))
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
            result = optimize(x -> objective(x, pk_param_values[:,i],ic50[i],volume[i],mult), [initial_guess[j]], LBFGS())
            input_dose = Optim.minimizer(result)
            patient_doses[j, i] = input_dose[1]
        end
    end
    return patient_doses
end