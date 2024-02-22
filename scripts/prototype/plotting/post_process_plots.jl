# filename = "results/optim/rg_probmask_test_ADAM_finitediff_temp2_lr0.1_2024-02-08_400patients.jls"
filename = "results/optim/rg_probmask_2AUCloss_2_test_ADAM_finitediff_temp2_lr0.1_2024-02-19_400patients.jls"

using Dates, StatsPlots, Plots, Serialization, OrderedCollections
today = Dates.today()

drug = "rg"
#
include("./model/$(drug)_pkpd2.jl")
include("./model/$(drug)_dosing2.jl")
include("./model/$(drug)_params.jl")
include("./scripts/setup/init_integrate.jl")
include("./assets/pk_params.jl")
include("./scripts/setup/generate_vp_lhs.jl")
include("./scripts/setup/compute_dose_bvp.jl")
include("./scripts/setup/precompute_scale.jl")
include("./scripts/setup/compute_outputs.jl")
include("./src/setup.jl")
include("./utilities/utils.jl")

# filename = "results/optim/rg_probmask_2AUCloss_test_ADAM_finitediff_temp2_lr0.1_2024-02-17_400patients.jls"
patients = deserialize(open(filename, "r"))



# Extracting measures and parameters
losses = [patient.output_measures["optimal"].loss for patient in patients] 
ftvs = [patient.output_measures["optimal"].ftv for patient in patients]    
drug_aucs = [patient.output_measures["optimal"].drug_auc for patient in patients]
tumor_aucs = [patient.output_measures["optimal"].tumor_auc for patient in patients]

if drug == "rg"
    pars = keys(RG_params)
    criteria = [("Loss", losses), ("FTV", ftvs), ("Drug AUC", drug_aucs), ("Tumor AUC", tumor_aucs)]
    for par in pars
        par_values = [patient.ode_parameters[indexin([par], param_order)[1]] for patient in patients]
        push!(criteria, (par, par_values))
    end
    # cl2_values = [patient.ode_parameters[indexin(["Cl2"], param_order)[1]] for patient in patients]
    # vpla_values = [patient.ode_parameters[indexin(["Vpla"], param_order)[1]] for patient in patients]
    # criteria = [("Loss", losses), ("FTV", ftvs), ("Drug AUC", drug_aucs), ("Tumor AUC", tumor_aucs), ("Cl2", cl2_values), ("Vpla", vpla_values)]
else
    kel_values = [patient.ode_parameters[indexin(["kel"], param_order)[1]] for patient in patients]
    k23_values = [patient.ode_parameters[indexin(["k23"], param_order)[1]] for patient in patients]
    criteria = [("Loss", losses), ("FTV", ftvs), ("Drug AUC", drug_aucs), ("Tumor AUC", tumor_aucs), ("kel", kel_values), ("k23", k23_values)]
end

opt = reduce(hcat, [patient.output_measures["optimal"].doses for patient in patients])
discrete_palette = reverse(cgrad(:bilbao, 20, categorical=true))
slate = zeros(Int(rg_dosetimes[end]/24)+1, length(patients))
for (criterion_name, criterion_values) in criteria
    # Sorting the patients based on the criterion
    sorted_indices = sortperm(criterion_values)
    sorted_opts = opt[:, sorted_indices]
    slate[Int.(rg_dosetimes./24).+1,:] = sorted_opts

    # Generate heatmap
    heatmap_title = "Sorted by " * criterion_name
    heatmap_filename_png = "./results/drug_matrix/$(drug)_heatmap_2AUC_2__$(today)_" * lowercase(criterion_name) * "_.png"
    heatmap_filename_svg = "./results/drug_matrix/$(drug)_heatmap_2AUC_2__$(today)_" * lowercase(criterion_name) * "_.svg"
    # cgrad(:matter, 10, categorical = true)palette(:acton10, 5)
    # Creating the heatmap
    heatmap_plot = heatmap(slate', color=discrete_palette, xlabel="Day", ylabel="Patient", title=heatmap_title, size=(600, 1200))

    # Adding annotations
    num_patients = size(sorted_opts, 2)
    highest_value_patient = num_patients  # Last in sorted array
    lowest_value_patient = 1  # First in sorted array

    annotate!(12, highest_value_patient+7, text("High " * criterion_name, 8, :right, :top))
    annotate!(12, lowest_value_patient-7, text("Low " * criterion_name, 8, :right, :bottom))

    # Save the heatmap
    savefig(heatmap_plot, heatmap_filename_png)
    savefig(heatmap_plot, heatmap_filename_svg)
end


# gradation=[1e-5, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0-1e-5]
# patient_doses = compute_dose_effect(population, gradation)












function print_checks(idx)
    println("DRUG AUC FROM, $(patients[idx].output_measures["none"].drug_auc) TO  $(patients[idx].output_measures["max"].drug_auc) COMPARE TO at opt  $(patients[idx].output_measures["optimal"].drug_auc)")
    println("TUMOR AUC FROM,  $(patients[idx].output_measures["none"].tumor_auc) TO  $(patients[idx].output_measures["max"].tumor_auc) COMPARE TO at opt  $(patients[idx].output_measures["optimal"].tumor_auc)")
    println("LOSS FROM,  $(patients[idx].output_measures["none"].loss) TO  $(patients[idx].output_measures["max"].loss) COMPARE TO at opt  $(patients[idx].output_measures["optimal"].loss)")
end


#find patient with negative loss
for i in 1:length(patients)
    if patients[i].output_measures["optimal"].loss < 0
        println("Patient $i has negative loss")
        print_checks(i)
    end
end

relu(x::Real) = max(zero(x), x)
function loss(θ, ode_params, scaling)
    # const_dose = scaling[3]/mult_num
    drug_min_scaling = scaling[4]#*num_dose_times
    drug_max_scaling = scaling[3]#*num_dose_times

    # # gumbel_noise = -log.(-log.(rand(length(set), num)))
    # gumbel_noise = -log.(-log.(rand(length(possible_doses), num_dose_times)))
    # sample_pmf = exp.((λ + gumbel_noise) ./ temp)
    # sample_pmf ./= sum(sample_pmf, dims=1)
    # # θ = ((sample_pmf' * dose_mults) .* const_dose)'
    # θ = (sample_pmf' * possible_doses)'

    sol_temp = solve(prob, Rodas4P2(), p=[ode_params..., θ...], callback=hit, sensealg=nothing)

    println("[end-1, end]", sol_temp[end-1,end] )

    cell = (sol_temp[end,end] - scaling[1])/(scaling[2]-scaling[1])
    drug = (sol_temp[end-1,end]-drug_min_scaling)/(drug_max_scaling-drug_min_scaling)
    loss = cell + drug
    return loss
end

idx=1
scaler = [patients[idx].output_measures["max"].tumor_auc, patients[idx].output_measures["none"].tumor_auc, patients[idx].output_measures["max"].drug_auc, patients[idx].output_measures["none"].drug_auc]
loss(patients[idx].output_measures["optimal"].doses, patients[idx].ode_parameters, patients[idx].scaling)
loss(patients[idx].output_measures["optimal"].doses, patients[idx].ode_parameters, scaler)
print_checks(idx)


sol_temp = solve(prob, Rodas4P2(), p=[patients[idx].ode_parameters..., patients[idx].output_measures["optimal"].doses...], callback=hit, sensealg=nothing)
sol_temp[states["cAUC"],end]
sol_temp[end,end]
sol_temp[states["PlaRGAUC"],end]

drug_auc = trapezoidal_rule(sol_temp.t, sol_temp[states["PlaRG"],:])
tumor_auc = trapezoidal_rule(sol_temp.t, sol_temp[states["C"],:])

cell = (sol_temp[end,end] - scaler[1]) / (scaler[2]-scaler[1])
(sol_temp[end,end] - scaler[1])
(scaler[2]-scaler[1])

drug = (sol_temp[end-1,end] - scaler[4]) / (scaler[3]-scaler[4])

scaler[3] - scaler[4]
(sol_temp[end-1,end] - scaler[4])

println(patients[idx].output_measures["optimal"].doses)




# Assuming `patients` is an array of `Patient` structs as defined
# and `conditions` is not predefined, we need to extract all unique conditions first

# Extract all unique conditions
unique_conditions = Set{String}()
for patient in patients
    for cond in keys(patient.output_measures)
        if  !occursin("avg", cond)
            push!(unique_conditions, cond)
        end
    end
end
conditions = collect(unique_conditions)

conditions = [ "1.0e-5effect"
               ,"1.0e-5avg_effect"
               ,"0.1effect"
               ,"0.1avg_effect"
               ,"0.25effect"
               ,"0.25avg_effect"
               ,"0.5effect"
               ,"0.5avg_effect"
               ,"0.75effect"
               ,"0.75avg_effect"
               ,"0.9effect"
               ,"0.9avg_effect"
               ,"0.99999effect"
               ,"0.99999avg_effect"
               ,"random"
               ,"optimal"]


conditions = [ "1.0e-5avg_effect"
              ,"0.1avg_effect"
              ,"0.25avg_effect"
              ,"0.5avg_effect"
              ,"0.75avg_effect"
              ,"0.9avg_effect"
              ,"0.99999avg_effect"
              ,"random"
              ,"optimal"]



# Now place "random" and "optimal" at the desired positions
# Assuming "random" should be at the very end and "optimal" just before "random"
# Note: This step might be adjusted based on your exact needs for placing "random" and "optimal"
final_sorted_conditions = filter(x -> x != "random" && x != "optimal", sorted_conditions)
push!(final_sorted_conditions, "optimal")
push!(final_sorted_conditions, "random")

final_sorted_conditions


# Initialize dictionaries to store the data for all conditions
loss_data = OrderedDict(cond => Float64[] for cond in conditions)
ftv_data = OrderedDict(cond => Float64[] for cond in conditions)
drug_auc_data = OrderedDict(cond => Float64[] for cond in conditions)
tumor_auc_data = OrderedDict(cond => Float64[] for cond in conditions)

# Extract data from patients for all conditions
for patient in patients
    for cond in conditions
        if haskey(patient.output_measures, cond)
            push!(loss_data[cond], patient.output_measures[cond].loss)
            push!(ftv_data[cond], patient.output_measures[cond].ftv)
            push!(drug_auc_data[cond], patient.output_measures[cond].drug_auc)
            push!(tumor_auc_data[cond], patient.output_measures[cond].tumor_auc)
        end
    end
end

# Function to plot histograms
function plot_histograms(data_dict, title)
    p = plot(legend=:topright, title=title)
    for (i, cond) in enumerate(conditions)
        histogram!(p, data_dict[cond], bins=30, alpha=0.5, label=cond, color=colors[i])
    end
    return p
end
# Function to create a violin plot for a given metric
# Define a pastel color palette
# Define a pastel color palette
pastel_colors = [:lightpink, :lightgreen, :lightgoldenrodyellow, :lavender, :lightblue]
# Modern color palette and their darker versions
modern_colors = [:dodgerblue, :coral, :mediumseagreen, :orchid, :goldenrod, :steelblue]
darker_colors = [:darkblue, :darkred, :darkgreen, :darkorchid, :darkgoldenrod, :darkslateblue]

gr()

# Function to create a combined box and violin plot with adjusted aesthetics
function create_combined_plot(metric_data, title)
    p = plot(title=title, legend=false, grid=false)

    for (i, cond) in enumerate(conditions)
        # Create box plot with reduced width and darker fill
        # Adjust the width parameter here to control the box width
        box_width = 0.05 # Adjust this value as needed
        # boxplot!(p, [cond], metric_data[cond], width=box_width, color=darker_colors[i], linecolor=:black, fillalpha=0.3, outliers_marker=:asterisk, outliers_color=:red, label=false)

        # Overlay with violin plot
        # violin!(p, [cond], metric_data[cond], color=modern_colors[i], alpha=0.7, label=false)
        violin!(p, [cond], metric_data[cond], alpha=0.7, label=false, xrotation=45)
    end

    return p
end
# Create combined plots for each metric
# Create combined plots for each metric
p1 = create_combined_plot(loss_data, "Loss Distribution")
p2 = create_combined_plot(ftv_data, "FTV Distribution")
p3 = create_combined_plot(drug_auc_data, "Drug AUC Distribution")
p4 = create_combined_plot(tumor_auc_data, "Tumor AUC Distribution")

# Function to save plot in multiple formats
function save_plot(plot, base_filename)
    savefig(plot, "./results/violin/$(drug)_2AUC_avg_$(today)_" * base_filename * ".png")
    savefig(plot, "./results/violin/$(drug)_2AUC_avg_$(today)_" * base_filename * ".svg")
end

# Save the plots
save_plot(p1, "Loss_Distribution")
save_plot(p2, "FTV_Distribution")
save_plot(p3, "Drug_AUC_Distribution")
save_plot(p4, "Tumor_AUC_Distribution")


function create_combined_plot(metric_data, title; log_scale=false)
    if log_scale
        p = plot(title=title, legend=false, grid=false, yscale=:log10)
    else
        p = plot(title=title, legend=false, grid=false)
    end

    for (i, cond) in enumerate(conditions)
        # Overlay with violin plot
        violin!(p, [cond], metric_data[cond], alpha=0.7, label=false, xrotation=45)
    end

    return p
end

p5 = create_combined_plot(drug_auc_data, "Drug AUC Distribution", log_scale=true)
p6 = create_combined_plot(tumor_auc_data, "Tumor AUC Distribution", log_scale=true)
p7 = create_combined_plot(ftv_data, "FTV Distribution", log_scale=true)
p8 = create_combined_plot(loss_data, "Loss Distribution", log_scale=true)
save_plot(p5, "log_Drug_AUC_Distribution")
save_plot(p6, "log_Tumor_AUC_Distribution")
save_plot(p7, "log_FTV_Distribution")
save_plot(p8, "log_Loss_Distribution")







#make an interactive histogram of 0.5avgeffect drug_auc

p = plot(legend=:topright, title="Drug AUC Distribution")
histogram!(p, log.(drug_auc_data["0.5avg_effect"]), bins=30, alpha=0.5, label="0.5avg_effect", color=:lightblue)


using Plots

# Function to create a combined violin and scatter plot
function create_violin_scatter_plot(metric_data, title; point_size=1, jitter_width=0.9)
    p = plot(title=title, legend=false, grid=false, yscale=:log10)

    conditions = collect(keys(metric_data)) # Ensure conditions are in the desired order
    n_conditions = length(conditions)
    
    # Prepare colors for better visualization
    colors = palette(:tab10, n_conditions)
    
    for (i, cond) in enumerate(conditions)
        # Calculate positions for scatter points with jitter
        xs = fill(i, length(metric_data[cond])) .+ (rand(length(metric_data[cond])) .- 0.5) .* jitter_width

        # Violin plot for the distribution
        violin!(p, [i], metric_data[cond], label=false, color=colors[i], alpha=0.4)

        # Scatter plot for individual points
        scatter!(p, xs, metric_data[cond], label=false, color=colors[i], markersize=point_size, markerstrokealpha=0, alpha=0.7)
    end

    # Customize axes
    xticks!(p, 1:n_conditions, conditions)
    xaxis!(p, rotation=45, xlabel="Condition")
    yaxis!(p, ylabel="Metric Value")
    
    return p
end

# Example usage with one of your metrics
# Assuming `metric_data` is one of your data dictionaries like `loss_data`, `ftv_data`, etc.
p = create_violin_scatter_plot(drug_auc_data, "drug auc disctributions")
display(p)




using Plots
plotlyjs() # Switch to the PlotlyJS backend for interactivity

# Define a structure to hold metric values and patient indices
struct MetricPoint
    value::Float64
    patient_idx::Int
    jittered_x::Float64
end

# Function to create an interactive combined violin and scatter plot with hover metadata
function create_interactive_violin_scatter_plot(conditions, patients, property, scale; point_size=1.5, jitter_width=0.9)
    if scale == :log10
        p = plot(title="$(property)", legend=false, grid=false, yscale=:log10)
    else
        p = plot(title="$(property)", legend=false, grid=false)
    end

    n_conditions = length(conditions)
    colors = palette(:tab10, n_conditions)
    
    for (i, cond) in enumerate(conditions)
        metric_points = MetricPoint[]
        
        for patient in patients
            if haskey(patient.output_measures, cond)
                value = getproperty(patient.output_measures[cond], property)
                idx = patient.idx
                jittered_x = i + (rand() - 0.5) * jitter_width
                push!(metric_points, MetricPoint(value, idx, jittered_x))
            end
        end
        
        # Sort metric points by jittered x to keep the association correct
        sort!(metric_points, by = m -> m.jittered_x)

        # Prepare hover text
        hover_texts = ["Patient: $(m.patient_idx)" for m in metric_points]
        xs = [m.jittered_x for m in metric_points]
        ys = [m.value for m in metric_points]

        # Check associations
        # for i in 1:length(xs)
        #     println("X: ", xs[i], ", Y: ", ys[i], ", Hover: ", hover_texts[i])
        # end

        # Violin plot for the distribution
        violin!(p, [i], ys, label=false, color=colors[i], alpha=0.5, hovertext="")

        # Interactive scatter plot for individual points with hover information
        scatter!(p, xs, ys, label=false, color=colors[i], markersize=point_size, markerstrokealpha=0, alpha=0.7, 
                 hovertext=hover_texts)
    end

    xticks!(p, 1:n_conditions, conditions)
    xaxis!(p, rotation=45, xlabel="Condition")
    yaxis!(p, ylabel=property)
    
    return p
end

# Example usage
p = create_interactive_violin_scatter_plot(conditions, patients, :loss, nothing)
display(p)

#find index in patients of a patient with idx property = 378
idx = findfirst(x->x.idx==378, patients)

#make all the plots and save them
for output in [:loss, :ftv, :drug_auc, :tumor_auc]
    p = create_interactive_violin_scatter_plot(conditions, patients, output, :log10)
    savefig(p, "./results/violin/$(drug)_$(output)_violin_log_scatter_plot.html")
end







































function death_time(v_init)
    v_init = max(v_init, 0.0)
    log((log(v_init/K)/log(113/K)))/r
 end
 
 # find the volume at a time where drug doages are zero in the blood
 # use a callback to stop the simulation when the amount of u4 adn 7 are both below 1e-9
 function find_volume_at_zero_dose(u, t, integrator)
    tmz_plasma = integrator.u[4]
    drug_plasma = integrator.u[7]
    (tmz_plasma < 1e-8 && drug_plasma < 1e-8) && integrator.t > (end_time+7.0)*hours#drug_dosetimes[end]#*hours
        # terminate!(integrator)
    # end
 end
 function affect!(integrator)
    integrator.terminate!
 end
 
 zero_dose_ = ContinuousCallback(find_volume_at_zero_dose, affect!)
 cbset = CallbackSet(hit, zero_dose_)
 
 u0 = zeros(Float64, 9)
 u0[1] = 17.7
 tspan = (0.0, end_time + 7.0) .* hours
#  ode_params = [ode_params; [0.1,1.0]]
 p = [ode_params; doses]
 
 prob = ODEProblem(pk_pd!, u0, tspan, p)
 sys = modelingtoolkitize(prob)
 sys = structural_simplify(sys)
 prob_jac = ODEProblem(sys, u0, tspan, p, jac=true)
 sol = solve(prob_jac, callback=cbset, saveat=0.5)
 
 
 using Survival
 # Initialize an empty plot
 plt = plot(xlabel="Time (Days)", ylabel="Survival Probability", title="Kaplan-Meier Curve for $(uppercase(drug)) Dosing", grid=false)
 nb_patients = size(patients, 1)
 for dose_spec in conditions
    println("\rDose Spec: $dose_spec")
    death_times = zeros(nb_patients)
    Threads.@threads for i in 1:nb_patients
        println("\rPatient: $i")
        dosing = patients[i].output_measures[dose_spec].doses
        ode_p = patients[i].ode_parameters
        p = [ode_p;dosing]
        p_prob = remake(prob_jac, p=p)
        p_sol = solve(p_prob, callback=cbset, saveat=1)#, abstol=1e-10, reltol=1e-10,dtmax=1)#, alg_hints=[:stiff])
        sols = Array(p_sol)
        # println("Patient $i terminated at $(sols[1,end])")
        # println("Patient $i stopped at time $(p_sol.t[end]) vs $((end_time+7)*24) (hours)")
        #assuming it terminated when the drug was zero_dose_
        death_times[i] = death_time(sols[1,end]) + (end_time+7)*24
    end
 
    # Create an event indicator vector, all ones (1) since all are events
 # Create an event indicator vector (all ones if all are death events)
 #     event_indicators = ones(Int, length(death_times))
 
 #     # Create a DataFrame with the required structure
 #     data = DataFrame(time = death_times, event = event_indicators)
 
 #     # Fit the Kaplan-Meier model
 #     km = fit(KaplanMeier, data.time, data.event)
 
 # # Extract survival probabilities and time points
 #     survival_probabilities = km.survival
 #     # Get unique event times from the original death_times data
 #     unique_times = unique(sort(death_times))
 #     # Ensure that the lengths of unique_times and survival_probabilities are the same
 #     # This is necessary because the Kaplan-Meier method may have fewer points than unique_times
 #     if length(unique_times) != length(survival_probabilities)
 #         unique_times = unique_times[1:length(survival_probabilities)]
 #     end
 
 
 #  # Plot the Kaplan-Meier curve
 #     plot!(plt,unique_times, survival_probabilities, label="$dose_spec")
    sorted_death_times = sort(death_times)
 
    # Initialize variables
    n = length(sorted_death_times)
    survival_probabilities = []
    times = []
 
    # Calculate survival probabilities
    for (i, time) in enumerate(sorted_death_times)
        if i == 1 || time != sorted_death_times[i-1]
            survival_probability = (n - i + 1) / n
            push!(survival_probabilities, survival_probability)
            push!(times, time)
        end
    end
 
    # Plot the Kaplan-Meier-type curve
    plot!(plt, times./24, survival_probabilities, label="$(dose_spec) Doses")
 
 end
 
 display(plt)
 
 savefig(plt, "./results/$(drug)_kaplan_meier_2AUC.svg")
 savefig(plt, "./results/$(drug)_kaplan_meier_2AUC.png")
 
 
