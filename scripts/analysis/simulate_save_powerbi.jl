# Ensure that two arguments have been provided
if length(ARGS) != 2
    println("Usage: julia args_example.jl <drug_name> <directory_path>")
    exit(1)
end
directory = ARGS[2]
using Pkg
prefix = directory#pwd()#"/scratch/c/craigm/popsicle/gbm/DosageOptimization.jl"
Pkg.activate(prefix)
Pkg.instantiate()



drug = "gdc" #ARGS[1]

using Plots, DifferentialEquations, MAT, Statistics, Serialization, ModelingToolkit
using HDF5
pyplot()

include("../../utilities/utils.jl")

include("../../models/$(drug)_pkpd2.jl")
include("../../models/$(drug)_dosing2.jl")
include("../../models/$(drug)_params.jl")

# population = deserialize("./assets/$(drug)_vp3.jls")
# population = deserialize("./assets/$(drug)_vp4_subset.jls")


population = deserialize("./assets/$(drug)_vp5.jls")
optimum = deserialize("./results/optimum/optima_$(drug)_lfbgs_auc_maxinit_2023-11-17_dms.jls")
opt = Array{Float64, 2}(optimum[1:end-2,:])
loss_bit = optimum[end,:]

# opt = deserialize("./assets/$(drug)_opt_subset.jls")
# min_tumor = deserialize("./assets/$(drug)_min_tumor_max_dose.jls")

tumor_min_scaling = deserialize("./assets/$(drug)_max_dose_min_tumor_v5.jls")
tumor_max_scaling = deserialize("./assets/$(drug)_min_dose_max_tumor_v5.jls")


# population = deserialize("./assets/$(drug)_vp4.jls")
# optimum = deserialize("./results/optimum/optima_$(drug)_lfbgs_auc_maxinit_2023-11-10_dms.jls")
# function find_populated_columns(matrix)
#     populated_indices = []
#     for col in 1:size(matrix, 2)
#         column_populated = false
#         for row in 1:size(matrix, 1)
#             if isassigned(matrix, row, col) && matrix[row, col] != nothing # Update the condition as needed
#                 column_populated = true
#                 break
#             end
#         end
#         if column_populated
#             push!(populated_indices, col)
#         end
#     end
#     return populated_indices
# end
# populated_columns = find_populated_columns(optimum)
# selected_patients = population[:,populated_columns]
# population = selected_patients

# opt = Array{Float64, 2}(optimum[1:end-2,populated_columns])
# loss_bit = optimum[end,populated_columns]

# tumor_min_scaling = deserialize("./assets/$(drug)_max_dose_min_tumor_v4.jls")[populated_columns]
# tumor_max_scaling = deserialize("./assets/$(drug)_min_dose_max_tumor_v4.jls")[populated_columns]

#make a vector of tuples of the form (x,y)
scaling = [collect(s) for s in zip(tumor_min_scaling, tumor_max_scaling)]

drug_min_scaling = min_drug_dosage
drug_max_scaling = sum(doses)


u0 = [17.7,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
tspan = (0,end_time+7).*hours
ode_params = [ode_params; 1.0; 0.1]
p = [ode_params;doses]
prob = ODEProblem(pk_pd!,u0,tspan,p)
sys = modelingtoolkitize(prob)
sys = structural_simplify(sys)
prob_jac = ODEProblem(sys, u0, tspan, p, jac=true)
sol = solve(prob_jac, callback=hit, saveat=0.5)

nb_patients = size(population, 2)
rng = Random.default_rng()
Random.seed!(rng, 1)

struct ConditionalOutput
    loss::Float64
    ftv::Float64
    drug_auc::Float64
    tumor_auc::Float64
    trajectory::Array{Float64,2}
end

struct Patient
    idx::Int
    ode_parameters::Vector{Float64}
    minimal_tumor::Float64
    maximal_tumor::Float64
    optimal_doses::Vector{Float64}
    output_measures::Dict{String, ConditionalOutput}
end

# Function to get doses based on specification
function get_dose(specification, patient_idx)
    if specification == "Min"
        return fill(min_drug_dosage/90, 90)  # Assuming min_dose is defined
    elseif specification == "Quarter"
        return doses ./ 4  # Assuming max_dose is defined
    elseif specification == "Half"
        return doses ./ 2
    elseif specification == "Max"
        return doses
    elseif specification == "Random"
        return rand(min_drug_dosage/90:doses[1], 90)
    elseif specification == "Optim"
        return opt[:,patient_idx]
    end
end

function loss(θ, ode_params)
    int_sol = solve(prob_jac, nothing, p=[ode_params;θ], callback=hit, sensealg=nothing)
    cell = (int_sol[end,end]-ode_params[end-1])/(ode_params[end]-ode_params[end-1])
    drug = (sum(abs,θ)-drug_min_scaling)/(drug_max_scaling-drug_min_scaling)
    return  cell + drug
end

patients = Vector{Patient}(undef, nb_patients)
Threads.@threads for i in 1:nb_patients
    println("\rPatient: $i")
    out = Dict{String, ConditionalOutput}()
    for dose_spec in ["Optim", "Min", "Quarter", "Half", "Max", "Random"]
        println("\rDose Spec: $dose_spec")
        dosing = get_dose(dose_spec, i)
        min_tumor = scaling[i][1]
        max_tumor = scaling[i][2]
        ode_p = [population[:, i];min_tumor;max_tumor]
        p = [ode_p;dosing]
        p_prob = remake(prob_jac, p=p)
        p_sol = solve(p_prob, callback=hit, saveat=1)#, abstol=1e-10, reltol=1e-10,dtmax=1)#, alg_hints=[:stiff])
        sols = Array(p_sol)
        c_solutions = hcat(sols[1,:],p_sol.t)
        drug_auc = trapezoidal_rule(p_sol.t, p_sol[7,:])
        ftv = sols[1,end]
        tumor_auc = sols[end,end]
        loss = ((tumor_auc-min_tumor)/(max_tumor-min_tumor) 
                + (sum(abs,dosing)-drug_min_scaling)/(drug_max_scaling-drug_min_scaling))
        #check that the computed loss is a Float64 and between 0.0 and 2.0
        if !((typeof(loss) <: Float64) && (0.0 <= loss <= 2.0))
            println("Loss is not a Float64 or is not between 0.0 and 2.0")
            println("Loss: $loss")
            println("Patient: $i")
            println("Dose Spec: $dose_spec")
            println("Dose: $dosing")
            println("Tumor AUC: $tumor_auc")
            println("Drug AUC: $drug_auc")
            println("FTV: $ftv")
            println("Min Tumor: $min_tumor")
            println("Max Tumor: $max_tumor")
            println("Optimum: $opt[:,i]")
            println("Optimum Loss: $(loss_bit[i])")
            println("Optimum Tumor AUC: $(sol[end,end])")
        end
        out[dose_spec] = ConditionalOutput(loss, ftv, drug_auc, tumor_auc, c_solutions)
    end
    patients[i] = Patient(i, population[:, i], scaling[i][1], scaling[i][2], opt[:,i], out)
end
#save the patient vector
serialize("./assets/$(drug)_patients_struct_subset.jls", patients)
using CSV, DataFrames

# 1. Patients Table
function generate_patients_csv(patients, param_order, filename="$(drug)_patients.csv")
    # Initialize an empty DataFrame with desired columns
    cols = [:patient_idx, :minimal_tumor, :maximal_tumor]
    append!(cols, Symbol.(param_order))
    for i in 1:90
        push!(cols, Symbol("optimal_dose_", i))
    end
    df_patients = DataFrame(;[col => Float64[] for col in cols]...)
    
    for patient in patients
        row = [patient.idx, patient.minimal_tumor, patient.maximal_tumor]
        # Add ODE parameters
        append!(row, patient.ode_parameters)
        # Add optimal doses-
        append!(row, patient.optimal_doses)
        push!(df_patients, row)
    end
    
    CSV.write(filename, df_patients)
end

# 2. ConditionalOutput Table
function generate_conditional_output_csv(patients, filename="$(drug)_outputs.csv")
    rows = []
    for patient in patients
        for (measure, output) in patient.output_measures
            push!(rows,(patient_idx = patient.idx, 
                        dose_label = measure,
                        loss = output.loss,
                        ftv = output.ftv,
                        drug_auc = output.drug_auc,
                        tumor_auc = output.tumor_auc))
        end
    end
    df_conditional_outputs = DataFrame(rows)
    CSV.write(filename, df_conditional_outputs)
end

# 3. Trajectories Table
function generate_trajectory_csv(patients, filename="$(drug)_trajectories.csv")
    rows = []
    for patient in patients
        for (measure, output) in patient.output_measures
            for i in 1:size(output.trajectory,1)
                push!(rows, (patient_idx=patient.idx,
                             dose_label=measure,
                             time_point=output.trajectory[i, 2],  # Assuming second dimension is time
                             value=output.trajectory[i, 1]))
            end
        end
    end
    df_trajectories = DataFrame(rows)
    CSV.write(filename, df_trajectories)
end

# 4. Loss Table
function generate_loss_csv(patients, filename="$(drug)_losses.csv")
    rows = []
    for i in 1:length(patients)
        init_loss = loss_bit[i][1]
        final_loss = loss_bit[i][end]
        many = length(loss_bit[i])
        loss_auc = sum(loss_bit[i])
        stopped_early = many < 100
        push!(rows,(patient_idx=i,
                    init_loss=init_loss,
                    final_loss=final_loss,
                    count=many,
                    loss_auc = loss_auc,
                    early_stop = stopped_early))
    end
    df_loss = DataFrame(rows)
    CSV.write(filename, df_loss)
end

# Call the functions to generate CSVs
generate_patients_csv(patients, param_order, "$(drug)_patients_subset.csv")
generate_conditional_output_csv(patients,   "$(drug)_outputs_subset.csv")
generate_trajectory_csv(patients  ,   "$(drug)_trajectories_subset.csv")
generate_loss_csv(patients      ,   "$(drug)_losses_subset.csv")