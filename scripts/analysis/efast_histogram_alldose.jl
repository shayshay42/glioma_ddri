using GlobalSensitivity
using Serialization
using ModelingToolkit
using HDF5
using StatsPlots
using Dates

today = Dates.today()

# Define the model
drug = "gdc"
include("../../models/$(drug)_pkpd2.jl")
include("../../models/$(drug)_dosing2.jl")
include("../../models/$(drug)_params.jl")
include("../../utilities/utils.jl")

# population = deserialize("./assets/$(drug)_vp5.jls")
optimum = deserialize("./results/optimum/optima_$(drug)_lfbgs_auc_maxinit_2023-11-17_dms.jls")
opt = Array{Float64, 2}(optimum[1:end-2,:])
loss_bit = optimum[end,:]

#deserialize optimum .jls from results
# optimum = 
# #mkae it into a Float64 array
# opt = Array{Float64, 2}(optimum[1:end-2,:])
#save it in the same directory as a HDF5
# h5write("./results/optimum/optima_$(drug)_lfbgs_auc_oct24_dms.h5", "optimum", opt)

#average across columns to make a vector of average doses per day
avg_doses = mean(opt, dims=2)

u0 = zeros(Float64, 9)
u0[1] = 17.7
tspan = (0.0, end_time + 7.0) .* hours
ode_params = [ode_params; [0.1,1.0]]
p = [ode_params; doses]

prob = ODEProblem(pk_pd!, u0, tspan, p)
sys = modelingtoolkitize(prob)
sys = structural_simplify(sys)
prob_jac = ODEProblem(sys, u0, tspan, p, jac=true)


include("../../assets/bounds4.jl")
println("Loaded the packages and variables...")

dosetimes = eval(Symbol("$(drug)_dosetimes"))

if drug == "rg"
    # Replace bounds for the RG and TMZ parameters
    drug_params = RG_params
    drug_lower = RG_lower
    drug_upper = RG_upper
else
    # Replace bounds for the RG and TMZ parameters
    drug_params = GDC_params
    drug_lower = GDC_lower
    drug_upper = GDC_upper
end
# List of parameter names in the order they are expected
# param_order = ["gamma_1", "psi", "C0", "D0", "r", "K", "BW", "IC50_1", "Imax_1", "IC50_2", "gamma_2", "Imax_2", "xi", "VD1", "Cl1", "k23", "ka1", "k32", "Cl2", "ka2", "Vpla", "Q", "Vtis"]

# TMZ_params = ["VD1", "Cl1", "k23"]
# Parameter values
param_values = ode_params

# Lower and upper bounds
# lb = copy(param_values)
# ub = copy(param_values)


# for (i, param) in enumerate(param_order)
#     if param in drug_params
#         param_index = findfirst(x -> x ==(param), drug_params)
#         lb[i] = drug_lower[param_index]
#         ub[i] = drug_upper[param_index]
#     elseif param in TMZ_params
#         param_index = findfirst(x -> x ==(param), TMZ_params)
#         lb[i] = TMZ_lower[param_index]
#         ub[i] = TMZ_upper[param_index]
#     end
# end

# quarter_dose = 0.25.*doses
# # Add dose to lower and upper bounds
# lb = vcat(lb, quarter_dose)
# ub = vcat(ub, quarter_dose)

# Find index of "BW" in param_order
bw_index = findfirst(x -> x == "BW", param_order)

# # Update lower and upper bounds for "BW"
# lb[bw_index] = 40
# ub[bw_index] = 140

new_p = copy(p)
function outit(p, dose, out)
    # Update the parameters
    indices = indexin([TMZ_params;drug_params;["BW"]], param_order)
    new_p[indices] = p
    new_p[length(ode_params)+1:end] = dose
    tmp_prob = remake(prob_jac, p=new_p)
    sol = solve(tmp_prob, callback=hit, dense=false)#, trajectories=size(p,2))

    if out == "Tumor Burden"
        return sol[end,end]
    elseif out== "Toxicity"
        return trapezoidal_rule(sol.t,sol[7,:])
    end
end

lb = [TMZ_lower;drug_lower;40]
ub = [TMZ_upper;drug_upper;140]

println("Starting eFAST!")
n_samples = 200
n_harmonics = 4
bounds = [[lb[i],ub[i]] for i in 1:length(lb)]


# Function to get doses based on specification
function get_dose(specification)
    if specification == "Min"
        return fill(min_drug_dosage/length(dosetimes), 90)  # Assuming min_dose is defined
    elseif specification == "Quarter"
        return doses ./ 4  # Assuming max_dose is defined
    elseif specification == "Half"
        return doses ./ 2
    elseif specification == "Max"
        return doses
    elseif specification == "Random"
        return rand(min_drug_dosage/length(dosetimes):doses[1], 90)
    elseif specification == "Optimal(mean)"
        return avg_doses
    end
end


function plot_histogram(res_tumor, res_drug, n_harmonics, n_samples, dose_spec, filename, drug)
    # Determine the number based on the drug type
    number = drug == "rg" ? "7112" : "0068"

    bw_index = indexin(["BW"], param_order)
    tmz_pk_params = indexin(TMZ_params, param_order)
    drug_pk_params = indexin(drug_params, param_order)
    params_indices = [tmz_pk_params; drug_pk_params; bw_index]
    params_names = [param_order[i] for i in params_indices]

    # Assuming params_names is the list of all parameter names on the x-axis
    # and it is already in the correct order as the columns of the heatmap
    tmz_params_count = 3  # The first three parameters are TMZ
    drug_params_count = length(params_names) - tmz_params_count - 1  # All except the first three and the last are Drug parameters
    dummy_params_count = 1  # The last parameter is Dummy

    # Data for heatmap
    heatmap_data = [getfield(res_tumor, Symbol("S1"));
                    getfield(res_tumor, Symbol("ST"));
                    getfield(res_drug, Symbol("S1"));
                    getfield(res_drug, Symbol("ST"))]

    p = heatmap(
        heatmap_data,
        xticks = (1:length(params_names), params_names),
        yticks = (1:4, ["Tumor Burden First Order", "Tumor Burden Total Order", "Drug Exposure First Order", "Drug Exposure Total Order"]),
        clims=(0, 1), 
        color=cgrad(:acton, rev=true), 
        colorbar_title="eFAST Index",
        title="Sensitivity to $(uppercase(drug))-$number PK Params, $dose_spec Dose, H=$n_harmonics, N=$n_samples",
        size = (800, 600)
    )


    # # Remove the legend
    # legend = false


    # # Y-axis brackets to indicate "Tumor" and "Drug"
    # y_bracket_top_position = 3.5  # Position for the top y-axis bracket
    # y_bracket_bottom_position = 1.5  # Position for the bottom y-axis bracket

    # # Draw y-axis brackets
    # plot!(p, [0, 5], [y_bracket_top_position, y_bracket_top_position], line = (:black, 2))
    # plot!(p, [0, 5], [y_bracket_bottom_position, y_bracket_bottom_position], line = (:black, 2))

    # # Add caps to the y-axis brackets
    # plot!(p, [0, 0], [y_bracket_top_position - 0.1, y_bracket_top_position + 0.1], line = (:black, 2))
    # plot!(p, [0, 0], [y_bracket_bottom_position - 0.1, y_bracket_bottom_position + 0.1], line = (:black, 2))

    # # X-axis brackets to the top
    # x_bracket_position = 4.5  # Position for the x-axis brackets above the heatmap

    # # Draw x-axis brackets for "TMZ" and "Drug" parameters
    # plot!(p, [0.5, tmz_params_count + 0.5], [x_bracket_position, x_bracket_position], line = (:black, 2))
    # plot!(p, [tmz_params_count + 0.5, length(params_names) - 0.5], [x_bracket_position, x_bracket_position], line = (:black, 2))

    # # Add caps to the x-axis brackets
    # plot!(p, [0.5, 0.5], [x_bracket_position - 0.1, x_bracket_position + 0.1], line = (:black, 2))
    # plot!(p, [tmz_params_count + 0.5, tmz_params_count + 0.5], [x_bracket_position - 0.1, x_bracket_position + 0.1], line = (:black, 2))
    # plot!(p, [length(params_names) - 0.5, length(params_names) - 0.5], [x_bracket_position - 0.1, x_bracket_position + 0.1], line = (:black, 2))

    # # Annotations for "Tumor" and "Drug" on the y-axis
    # annotate!(p, [(5.5, y_bracket_top_position, Plots.text("Tumor", :right, 10)),
    #               (5.5, y_bracket_bottom_position, Plots.text("Drug", :right, 10))])

    # Save the plot
    savefig(p, filename * ".png")
    savefig(p, filename * ".svg")
end

# Loop through dose specifications and outcome interests
for dose_spec in ["Min", "Quarter", "Half", "Max", "Random", "Optimal(mean)"]
    dose = get_dose(dose_spec)
    out_both=[]
    for outcome in ["Tumor Burden", "Toxicity"]
        # Update filename based on dose_spec and outcome_interest
        filename = "./results/efast/$drug/$(today)_$(drug)_$(dose_spec)dose_efast_$(n_samples)samples_$(n_harmonics)harmonics_$(outcome).jls"
        
        # Run eFAST analysis
        efast = gsa((p) -> outit(p, dose, outcome), eFAST(num_harmonics=n_harmonics), bounds, samples=n_samples)
        push!(out_both, efast)
        open(filename, "w") do file
            serialize(file, efast)
        end;
        println("Finished eFAST for $dose_spec dose and $outcome")
    end
    # Plot eFAST results
    filename = "./results/efast/$drug/$(today)_$(drug)_$(dose_spec)dose_efast_$(n_samples)samples_$(n_harmonics)harmonics"
    plot_histogram(out_both[1], out_both[2], n_harmonics, n_samples, dose_spec, filename, drug)
end

# #test by just loading the jls file
# for dose_spec in ["Min"]#, "Quarter", "Half", "Max", "Random", "Optimal(mean)"]
#     dose = get_dose(dose_spec)
#     out_both=[]
#     for outcome in ["Tumor Burden", "Toxicity"]
#         # Update filename based on dose_spec and outcome_interest
#         filename = "./results/efast/$drug/$(today)_$(drug)_$(dose_spec)dose_efast_$(n_samples)samples_$(n_harmonics)harmonics_$(outcome).jls"
        
#         # Run eFAST analysis
#         efast = deserialize(filename)
#         push!(out_both, efast)
#         println("Finished eFAST for $dose_spec dose and $outcome")
#     end
#     # Plot eFAST results
#     filename = "./results/efast/$drug/$(today)_$(drug)_$(dose_spec)dose_efast_$(n_samples)samples_$(n_harmonics)harmonics"
#     plot_histogram(out_both[1], out_both[2], n_harmonics, n_samples, dose_spec, filename, drug)
# end