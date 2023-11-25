function compute_loss_scaling(population, max_doses, min_doses)
    nb_patients = size(population, 2)
    max_dose_min_tumor = Vector{Float64}(undef,  nb_patients)
    min_dose_max_tumor = Vector{Float64}(undef,  nb_patients)

    for (j, dosage) in enumerate([max_doses, min_doses])
        Threads.@threads for i in 1:nb_patients
            # println("\rPatient: $i")
            p = [population[:, i];1.0;1.0;dosage]
            p_prob = remake(prob_jac, p=p)
            p_sol = solve(p_prob, callback=hit)#, abstol=1e-10, reltol=1e-10,dtmax=1)#, alg_hints=[:stiff])

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