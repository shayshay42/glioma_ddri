# filename = "results/optim/rg_probmask_test_ADAM_finitediff_temp2_lr0.1_2024-02-08_400patients.jls"
filename = "results/dms_optim/rg_probmask_2AUCloss_ADAM_finitediff_temp2_lr0.1_2024-02-21_100patients.jls"

using Dates, StatsPlots, Plots, Serialization, OrderedCollections
today = Dates.today()

drug = "rg"
#
include("../../model/$(drug)_pkpd2.jl")
include("../../model/$(drug)_dosing2.jl")
include("../../model/$(drug)_params.jl")
include("../../scripts/setup/init_integrate.jl")
include("../../assets/pk_params.jl")
include("../../scripts/setup/generate_vp_lhs.jl")
include("../../scripts/setup/compute_dose_bvp.jl")
include("../../scripts/setup/precompute_scale.jl")
include("../../scripts/setup/compute_outputs.jl")
include("../../src/setup.jl")
include("../../utilities/utils.jl")

# filename = "results/optim/rg_probmask_2AUCloss_test_ADAM_finitediff_temp2_lr0.1_2024-02-17_400patients.jls"
patients = deserialize(open(filename, "r"))

function compute_dose_effect(patients, gradation = [1e-5, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99999])

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
    patient_doses = zeros(length(gradation), size(patients, 2))
    unstable_patients = Dict()
    Threads.@threads for i in 1:length(patients)
        println("\rPatient: $i")
        ode_p = patients[i].ode_parameters
        # pk_param_names = collect(keys(pk_params))
        # pk_indices = indexin(pk_param_names, param_order)
        # pk_param_values= population[pk_indices, i]
        ic502 = patients[i].ode_parameters[indexin(["IC50_2"], param_order)][1]
        gamma2 = patients[i].ode_parameters[indexin(["gamma_2"], param_order)][1]
        imax2 = patients[i].ode_parameters[indexin(["Imax_2"], param_order)][1]
        psi=1
        if "Vpla" in param_order
            volume = patients[i].ode_parameters[indexin(["Vpla"], param_order)][1]
        elseif "V2" in param_order
            volume = patients[i].ode_parameters[indexin(["V2"], param_order)][1]
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
                # sol_temp = solve(prob_temp, Shooting(Rodas4P2()), callback=cb)
                sol_temp = solve(prob_temp, callback=cb)
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

patient_doses = compute_dose_effect(patients)
#show the unstable patients
println(patient_doses.retcodes)





dose_matrix = patient_doses.doses_matrix
scaling = 
#compute the metrics from the doses computed
Threads.@threads for i in 1:length(patients)
    @info "Completing Patient $i"
    output_measures = OrderedDict()

    for (j, effect) in enumerate(effect_keys)
        output_measures[effect] = get_outputs(population[:,i], ones(length(drug_dose_times)).*doses_matrix[j,i], scaling[i], drug)
    end

    for (j, dose) in enumerate(avg_dose_per_gradation)
        output_measures[string(gradation[j])*"avg_effect"] = get_outputs(population[:,i], ones(length(drug_dose_times)).*dose, avg_scale[i], drug)
    end

    min_dose = doses_matrix[1,i]
    max_dose = doses_matrix[end,i]
    random_doses = min_dose .+ (max_dose - min_dose) .* rand(length(drug_dose_times))
    output_measures["random"] = get_outputs(population[:,i], random_doses, scaling[i], drug)
    output_measures["none"] = get_outputs(population[:,i], zeros(length(drug_dose_times)), scaling[i], drug)
    output_measures["max"] = get_outputs(population[:,i], ones(length(drug_dose_times)).*single_max, scaling[i], drug)
    # Create a new Patient instance
    # patients[i] = Patient(i, population[:, i], scaling[i], zeros(length(drug_dose_times)), output_measures)
    patients[i] = Patient(i, population[:, i], scaling[i], output_measures)
end