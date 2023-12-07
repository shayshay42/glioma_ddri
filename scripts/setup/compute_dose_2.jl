using DifferentialEquations, Optim
using BoundaryValueDiffEq

function compute_dose2(population, pk_params)

    function dose_affect!(integrator)
        integrator.u[1] += integrator.u[4]
    end
    # cb = PeriodicCallback(dose_affect!, 24.0, initial_affect=true)

    event_times = collect(1:18).*24.0
    cb = PresetTimeCallback(event_times, dose_affect!)
    tspan = (0.0, event_times[end]+(10.0*24.0))
    
    p = rand(6)
    u0 = [0.0, 0.0, 0.0, 5000.0]
    prob = ODEProblem(drug_pk_dose!, u0, tspan, p)

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
        pk_param_values= population[pk_indices, i]
        ic50 = population[indexin(["IC50_2"], param_order), i][1]
        if "Vpla" in param_order
            volume = population[indexin(["Vpla"], param_order), i][1]
        elseif "V2" in param_order
            volume = population[indexin(["V2"], param_order), i][1]
        end
        # Boundary condition function

        for (j,mult) in enumerate(mults)
            function bc_drug_pk!(residual, sol, p, t)
                residual[1] = maximum(sol[2,:]) - ((ic50 * (1000 * volume))*mult) # Ensure the second state matches ic50 amount conversion at end time
                residual[2] = sol[1][2] # Ensure the seocnd state is 0 at the beginning
                residual[3] = sol[1][3] # Ensure the third state is 0 at the beginning
                residual[4] = sol[1][1]
            end
            u0_temp = [0.0, 0.0, 0.0, initial_guess[j]]
            p_temp = pk_param_values
            prob_temp = BVProblem(drug_pk_dose!, bc_drug_pk!, u0_temp, tspan, p_temp)
            sol_temp = solve(prob_temp, Shooting(KenCarp4()), callback=cb)
            input_dose = sol_temp[4, 1]
            patient_doses[j, i] = input_dose
        end
    end
    return patient_doses
end