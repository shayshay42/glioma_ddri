using Pkg
Pkg.activate(".")
Pkg.instantiate()

using DifferentialEquations, Serialization, ModelingToolkit
# pyplot()

function parse_my_args()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--drug"
            help = "which drug"
            default = "rg"
            arg_type = String
        "--filename"
            help = "output filename"
            default = "rg_1000_vp.jls"
            arg_type = String
    end

    return ArgParse.parse_args(s)
end

args = parse_my_args()
drug = args["drug"]
filename = args["filename"]

include("../../utilities/utils.jl")
include("init_integrate.jl")

# drug = "gdc"
# include("../../models/$(drug)_pkpd2.jl")
# include("../../models/$(drug)_dosing2.jl")
# include("../../models/$(drug)_params.jl")

population = deserialize("./assets/$filename")

# u0 = [17.7,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
# tspan = (0,end_time+7).*hours
# ode_params = [ode_params; 1.0; 1.0]
# p = [ode_params;doses]
# prob = ODEProblem(pk_pd!,u0,tspan,p)
# sys = modelingtoolkitize(prob)
# sys = structural_simplify(sys)
# prob_jac = ODEProblem(sys, u0, tspan, p)

nb_patients = size(population, 2)

max_dose_min_tumor = Vector{Float64}(undef,  nb_patients)
min_dose_max_tumor = Vector{Float64}(undef,  nb_patients)

for (j, dosage) in enumerate([doses, fill(min_drug_dosage/length(doses), length(doses))])
    Threads.@threads for i in 1:nb_patients
        println("\rPatient: $i")
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

#save resutls in assets
serialize("./assets/$(drug)_max_dose_min_tumor.jls", max_dose_min_tumor)
serialize("./assets/$(drug)_min_dose_max_tumor.jls", min_dose_max_tumor)
