using DifferentialEquations, Optim, ModelingToolkit, BoundaryValueDiffEq


"""
    compute_dose_effect(population, pk_params)
    works on a smaller pkpd that only contains the drug it self and the cells not hte TMZ

    the output is a named tuple with the following fields:
    grad: the gradation of the effect
    doses_matrix: the doses for each patient for each gradation
"""
function compute_dose_effect(population, gradation)

    function dose_affect!(integrator)
        SciMLBase.set_proposed_dt!(integrator, 0.1)
        integrator.u[2] += integrator.u[end] #indices folow the order of the state
    end
    # cb = PeriodicCallback(dose_affect!, 24.0, initial_affect=true)

    event_times = collect(1:18).*24.0
    cb = PresetTimeCallback(event_times, dose_affect!)
    tspan = (0.0, event_times[end]+(10.0*24.0))
 
    p = ode_params
    u0 = [17.7, 0.0, 0.0, 0.0, 5000.0]
    prob = ODEProblem(simple_pkpd!, u0, tspan, p)
    sys = modelingtoolkitize(prob)
    sys = structural_simplify(sys)
    func = ODEFunction(sys, jac=true)
    prob = ODEProblem(func, u0, tspan, p)
    # sol = solve(prob_jac, Rodas4P2(), callback=cb)

    # gradation = [0.00001, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99]
    # gradation = [0.1, 0.25, 0.5, 0.75, 0.9]
    # avg_human_surface_area=1.7
    # dosage guess for each multiplier of IC50
    # initial_guess = [20.0, 100.0, 500.0, 1000.0, 2000.0, 3000.0].*avg_human_surface_area # Adjust this as needed
    # initialize array to hold patient doses for each multiplier
    patient_doses = zeros(length(gradation), size(population, 2))
    unstable_patients = Dict()
    Threads.@threads for i in 1:size(population, 2)
        println("\rPatient: $i")
        ode_p = population[:, i]
        # pk_param_names = collect(keys(pk_params))
        # pk_indices = indexin(pk_param_names, param_order)
        # pk_param_values= population[pk_indices, i]
        ic502 = population[indexin(["IC50_2"], param_order), i][1]
        gamma2 = population[indexin(["gamma_2"], param_order), i][1]
        imax2 = population[indexin(["Imax_2"], param_order), i][1]
        psi=1
        if "Vpla" in param_order
            volume = population[indexin(["Vpla"], param_order), i][1]
        elseif "V2" in param_order
            volume = population[indexin(["V2"], param_order), i][1]
        end
        # Boundary condition function
        for (j,mult) in enumerate(gradation)
            function bc_drug_pk!(residual, sol, p, t)
                cPlaRG_sol = sol[3,:] ./ (1000 * volume)
                cPlaRG_sol = erelu.(cPlaRG_sol)
                exp2_sol = (cPlaRG_sol ./(psi*ic502)).^gamma2
                E_sol = (exp2_sol.*imax2) ./ (exp2_sol .+ 1)
                residual[1] = C0 - sol[1][1] # Ensure that the initial condition is fixed
                residual[2] = sol[1][2] # Ensure that the initial condition is fixed to 0
                residual[3] = sol[1][3]
                residual[4] = sol[1][4]
                residual[5] = maximum(E_sol) - mult # Ensure the plasma compartment matches ic50 amount conversion at end time
            end
            u0_temp = [17.7, 0.0, 0.0, 0.0, 1000.0] #1000.0 as an initial guess for the dose
            prob_temp = BVProblem(simple_pkpd!, bc_drug_pk!, u0_temp, tspan, ode_p)
            try
                sol_temp = solve(prob_temp, Shooting(Rodas4P2()), callback=cb)
                input_dose = sol_temp[1][5]
                patient_doses[j, i] = input_dose
                if sol_temp.retcode != 1 || sol_temp.retcode != :Success
                    # Store the patient index and return code
                    unstable_patients[i] = sol_temp.retcode
                end
            catch e
                println("Patient $i unstable at $mult")
                unstable_patients[i] = e
            end
        end
    end
    result = (grad=gradation, doses_matrix=patient_doses, retcodes=unstable_patients)
    return result
end