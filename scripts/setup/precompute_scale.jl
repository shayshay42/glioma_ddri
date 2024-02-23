function compute_loss_scaling(population)
    nb_patients = size(population, 2)
    max_dose_min_tumor = Vector{Float64}(undef,  nb_patients)
    min_dose_max_tumor = Vector{Float64}(undef,  nb_patients)

    for (j, dosage) in enumerate([ones(num_dose_times).*(single_max), zeros(num_dose_times)])
        Threads.@threads for i in 1:nb_patients
            # println("\rPatient: $i")
            p = [population[:, i]..., dosage...]
            p_prob = remake(prob, p=p)
            p_sol = solve(p_prob, Rodas4P2(), callback=hit)#, abstol=1e-10, reltol=1e-10,dtmax=1)#, alg_hints=[:stiff])

            sols = Array(p_sol)

            cAUC = sols[end,end]

            if j == 1
                max_dose_min_tumor[i] = cAUC
            else
                min_dose_max_tumor[i] = cAUC
            end
        end
    end
    return max_dose_min_tumor, min_dose_max_tumor
end


# same as above except different dose for each patient (effect based)

function compute_loss_scaling_effect_based(population, max_effect_dose_per_patient, min_effect_dose_per_patient, drug_dosetimes)
    nb_patients = size(population, 2)
    max_dose_min_tumor = Vector{Float64}(undef,  nb_patients)
    min_dose_max_tumor = Vector{Float64}(undef,  nb_patients)
    for (j, condition) in enumerate(["min dose", "max dose"])
        Threads.@threads for i in 1:nb_patients
            # println("\rPatient: $i")
            if condition == "min dose"
                dosage = min_effect_dose_per_patient[i]
            else
                dosage = max_effect_dose_per_patient[i]
            end
            p = [population[:, i]..., ones(length(drug_dosetimes)).*dosage...]
            p_prob = remake(prob, p=p)
            p_sol = solve(p_prob, Rodas4P2(), callback=hit)#, abstol=1e-10, reltol=1e-10,dtmax=1)#, alg_hints=[:stiff])

            sols = Array(p_sol)

            cAUC = sols[end,end]

            if condition == "min dose"
                min_dose_max_tumor[i] = cAUC
            else
                max_dose_min_tumor[i] = cAUC
            end
        end
    end
    return max_dose_min_tumor, min_dose_max_tumor
end


function compute_loss_scaling_AUC(vp_param_matrix, maxi_doses, mini_doses, drug)
    nb_patients = size(vp_param_matrix, 2)
    max_dose_min_tumor = Vector{Float64}(undef,  nb_patients)
    min_dose_max_tumor = Vector{Float64}(undef,  nb_patients)
    max_dose_max_drug = Vector{Float64}(undef,  nb_patients)
    min_dose_min_drug = Vector{Float64}(undef,  nb_patients)

    for (j, dosage) in enumerate([maxi_doses, mini_doses])
        Threads.@threads for i in 1:nb_patients
            # println("\rPatient: $i")
            p = [vp_param_matrix[:, i]..., dosage...]
            p_prob = remake(prob, p=p)
            p_sol = solve(p_prob, Rodas4P2(), callback=hit)#, abstol=1e-10, reltol=1e-10,dtmax=1)#, alg_hints=[:stiff])

            sols = Array(p_sol)

            cAUC = sols[states["cAUC"],end]
            plasmadrugAUC = sols[states["Pla"*uppercase(drug)*"AUC"],end]

            if j == 1
                max_dose_min_tumor[i] = cAUC
                max_dose_max_drug[i] = plasmadrugAUC
            else
                min_dose_max_tumor[i] = cAUC
                min_dose_min_drug[i] = plasmadrugAUC
            end
        end
    end
    #change to component array
    return max_dose_min_tumor, min_dose_max_tumor, max_dose_max_drug, min_dose_min_drug
end