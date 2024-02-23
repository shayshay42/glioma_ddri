function get_outputs(parameters, doses, scaling, drug_name)
    drug_min_scaling = scaling[4]*num_dose_times
    drug_max_scaling = scaling[3]*num_dose_times
    prob_temp = remake(prob, p=[parameters..., doses...])
    sol_temp = solve(prob_temp, Rodas4P2(), callback=hit, saveat=1)
    cell = (sol_temp[end,end] - scaling[1])/(scaling[2]-scaling[1])
    drug = (sum(abs,doses)-drug_min_scaling)/(drug_max_scaling-drug_min_scaling)
    loss = cell + drug

    ftv = sol_temp[states["C"],end]
    cell_auc = sol_temp[states["cAUC"],end]
    drug_compartment = "Pla"*uppercase(drug_name)
    drug_auc = trapezoidal_rule(sol_temp.t, sol_temp[states[drug_compartment],:])

    trajectory = hcat(sol_temp[states["C"],:], sol_temp[states[drug_compartment],:])
    
    return ConditionalOutput(doses, loss, ftv, drug_auc, cell_auc, trajectory)
end

function get_outputs_auc(parameters, doses, scaling, drug_name)
    drug_min_scaling = scaling[4]*num_dose_times
    drug_max_scaling = scaling[3]*num_dose_times
    prob_temp = remake(prob, p=[parameters..., doses...])
    sol_temp = solve(prob_temp, Rodas4P2(), callback=hit, saveat=1)

    ftv = sol_temp[states["C"],end]
    cell_auc = sol_temp[states["cAUC"],end]
    drug_compartment = "Pla"*uppercase(drug_name)
    drug_auc = sol_temp[states[drug_compartment*"AUC"],end]

    cell = (cell_auc - scaling[1])/(scaling[2]-scaling[1])
    drug = (drug_auc - scaling[4])/(scaling[3]-scaling[4])
    loss = cell + drug
    
    return ConditionalOutput(doses, loss, ftv, drug_auc, cell_auc)
end