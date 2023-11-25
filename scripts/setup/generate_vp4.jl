using Pkg
Pkg.activate(".")
Pkg.instantiate()

using Distributions, Random, Serialization, Statistics, ArgParse
# pyplot()

function parse_my_args()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--drug"
            help = "which drug"
            default = "rg"
            arg_type = String
        "--nbr"
            help = "number of patients"
            default = 1000
            arg_type = Int
        "--filename"
            help = "output filename"
            default = "default.jls"
            arg_type = String
        "--seed"
            help = "random seed"
            default = 1234
            arg_type = Int
    end

    return ArgParse.parse_args(s)
end

args = parse_my_args()
drug = args["drug"]
num_patients = args["nbr"]
filename = args["filename"]
seed = args["seed"]

include("../../utilities/utils.jl")

# # Script Arguments:
# # ARG 1: Drug Type ('rg' or 'gdc')
# # ARG 2: Number of Patients (integer)

# # Ensure two arguments are provided
# if length(ARGS) < 2
#     println("Usage: julia script_name.jl [drug_type] [num_patients]")
#     println("drug_type: 'rg' or 'gdc'")
#     println("num_patients: integer value")
#     exit(1)
# end

# # Validate and assign drug type
# drug = ARGS[1]
# if drug ∉ ["rg", "gdc"]
#     println("Invalid drug type. Please specify 'rg' or 'gdc'.")
#     exit(1)
# end

# # Validate and assign number of patients
# try
#     num_patients = parse(Int64, ARGS[2])
# catch
#     println("Invalid number of patients. Please specify an integer.")
#     exit(1)
# end

# drug = "rg"
include("../../model/$(drug)_pkpd2.jl")
include("../../model/$(drug)_dosing2.jl")
include("../../model/$(drug)_params.jl")

Random.seed!(seed)
# include("../assets/bounds4.jl")

#ode parameters and drug must be defined in the script that calls this file in include

param_values = ode_params

# Define parameters for RG
#           1000*60   60    1000 1000*60  1000
# RG_params = ["Cl2", "ka2", "Vpla", "Q", "Vtis"]
# RG_mu = convert(Array{Float64}, [0.00028, 0.00017, 0.0015, 0.00000054, 0.0000044])
# RG_omega = convert(Array{Float64}, [0.25, 0.73, 1.76, 2.25, 5.08])
# RG_mult = [60, 1000*60, 1000, 1000*60, 1000]

# # Define parameters for GDC
# #              60  1000    60    60    60
# GDC_params = ["ka2","V2" ,"kel","k12","k21"]
# GDC_mu = convert(Array{Float64}, [0.0039, 0.077*1e3, 0.029, 1.12, 2.86])
# GDC_omega = convert(Array{Float64}, [0.017, 0.036, 0.3, 0.064, 0.073])
# GDC_mult = [60, 1000, 60, 60, 60]

# # Define parameters for TMZ
# #                    
# TMZ_params = ["Cl1", "k23", "ka1"]
# TMZ_mu = convert(Array{Float64}, [10.0, 7.2 * 10^-4, 5.8])
# TMZ_sigma = convert(Array{Float64}, [0.047, 0.1649, 0.8961])
include("../../assets/pk_params.jl")

# Initialize population
population = repeat(param_values, 1, num_patients)

# Replace values for the TMZ parameters
TMZ_indices = indexin(TMZ_params, param_order)
#findall(x -> x in TMZ_params, param_order)
# TMZ_samples = hcat([rand(Truncated(Normal(TMZ_mu[i], TMZ_sigma[i]), 0.0, Inf), num_patients) for i in 1:length(TMZ_mu)]...)
# population[TMZ_indices, :] = TMZ_samples'

for (i, param) in enumerate(TMZ_params)
    sample = rand(LogNormal(log(TMZ_mu[i]), TMZ_sigma[i]), num_patients) #can be negative?
    # sample = TMZ_mu[i] .* exp.(X)
    population[TMZ_indices[i], :] = sample
end

# function filter_patient(ps)
#     u0 = [17.7,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
#     tspan = (0,end_time+7).*hours
#     p = [ps;doses]
#     prob = ODEProblem(pk_pd!,u0,tspan,p)
#     sol = solve(prob, callback=hit, saveat=0.1)
#     auc_drug = trapezoidal_rule(sol.t, sol[7,:])
#     if auc_drug > 1.0e5
#         return true
#     else
#         return false
#     end
# end

# Switch between drug
if drug == "rg"
    println("made for RG!")
    # Get indices for the RG parameters
    RG_indices = indexin(RG_params, param_order)
    #findall(x -> x in RG_params, param_order)

    # Replace values for the RG parameters
    # RG_samples = hcat([rand(Truncated(LogNormal(RG_mu[i], RG_sigma[i]), 0.0, Inf), num_patients) for i in 1:length(RG_mu)]...)
    

    for (i, param) in enumerate(RG_params)
        if param == "Q"
            sample = rand(LogitNormal(logit(RG_mu[i]), RG_sigma[i]), num_patients) .* RG_mult[i]
        else
            sample =   rand(LogNormal(  log(RG_mu[i]), RG_sigma[i]), num_patients) .* RG_mult[i]
        end
        population[RG_indices[i], :] = sample
    end

    # population[RG_indices, :] = RG_samples'

elseif drug == "gdc"
    println("made for GDC!")
    # Get indices for the GDC parameters
    GDC_indices = indexin(GDC_params, param_order)
    #findall(x -> x in GDC_params, param_order)

    # Replace values for the GDC parameters
    # GDC_samples = hcat([rand(Truncated(LogNormal(GDC_mu[i], GDC_sigma[i]), 0.0, Inf), num_patients) for i in 1:length(GDC_mu)]...)
    for (i, param) in enumerate(GDC_params)
        population[GDC_indices[i], :] = rand(LogNormal(  log(GDC_mu[i]), GDC_sigma[i]), num_patients) .* GDC_mult[i]
    end
    # population[GDC_indices, :] = GDC_samples'
end;

# Assign a default filename or use the third argument if provided
# filename = length(ARGS) >= 3 ? ARGS[3] : "$(drug)_vp.jls"

# Serialize the population
serialize("./assets/$filename", population)

# filter patients
# the following loop replaces parameters for patients that do not meet the criteria
# for example, if the AUC of TMZ is greater than 500000, then the patient is replaced
# with a new patient
# iteration_counter = zeros(num_patients)  # Initialize counter for each patient

# for i in 1:num_patients
#     while filter_patient(population[:, i])
#         iteration_counter[i] += 1  # Update counter for the patient
#         print("\rIteration number for patient $i: $(iteration_counter[i])")  # print iteration number in real-time on the same line
#         if drug == "rg"
#             for (j, param) in enumerate(RG_params)
#                 if param == "Q"
#                     X = rand(Normal(0, RG_omega[j]))
#                     m = logit(RG_mu[j]) + X
#                     sample = inv_logit(m)
#                 else
#                     X = rand(Normal(0, RG_omega[j])) #can be negative?
#                     sample = RG_mu[j] * exp(X)
#                 end
#                 population[RG_indices[j], i] = sample
#             end
#         elseif drug == "gdc"
#             for (j, param) in enumerate(GDC_params)
#                 X = rand(Normal(0, GDC_omega[j])) #can be negative?
#                 sample = GDC_mu[j] * exp(X)
#                 population[GDC_indices[j], i] = sample
#             end
#         end
#     end
#     #show new area under curve
#     u0 = [17.7,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
#     tspan = (0,end_time+7).*hours
#     p = [population[:, i];doses]
#     prob = ODEProblem(pk_pd!,u0,tspan,p)
#     sol = solve(prob, callback=hit, saveat=0.1)
#     auc_drug = trapezoidal_rule(sol.t, sol[7,:])
#     println("\rReplaced Successfully! New AUC: $auc_drug")
#     println("New Parameters: $(population[:, i])")
#     println("Total number of iterations for patient $i: $(iteration_counter[i])")
# end


# for i in 1:num_patients
#     while filter_patient(population[:, i])
#         if drug == "gdc"
#             for (j, param) in enumerate(GDC_params)
#                 X = rand(Normal(0, GDC_omega[j])) #can be negative?
#                 sample = GDC_mu[j] * exp(X)
#                 population[GDC_indices[j], i] = sample
#             end
#         end
#     end
#     println("Replaced Succesfully!")
# end

# serialize the population
# open("./assets/$(drug)_vp.jls", "w") do file
#     serialize(file, population)
# end

# ##################

#       "Cl1"                 "k23"              "ka1"                "ka2"                "V2"               "kel"
#      10.0126            0.00073245           8.68857              0.234253             77.3874             1.76836

# [9.998305012809254, 0.1295329656486665, 2.6845585910553553, 0.003898131004770741, 77.12735197205599, 0.03010744187902203]


# # Define parameters for GDC
# GDC_params = ["ka2","V2" ,"kel","k12","k21"]
# GDC_mu = convert(Array{Float64}, [0.0039, 0.077 * 10^3, 0.029, 1.12, 2.86])
# GDC_omega = convert(Array{Float64}, [0.017, 0.036, 0.3, 0.064, 0.073])

# # Define parameters for TMZ
# TMZ_params = ["VD1", "Cl1", "k23"]
# TMZ_mu = convert(Array{Float64}, [5.8, 10.0, 7.2 * 10^-4])
# TMZ_sigma = convert(Array{Float64}, [0.8961, 0.047, 0.1649])