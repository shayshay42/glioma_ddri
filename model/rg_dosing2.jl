#standard of care (SOC):
#20 - 1800mg /m2
#10 days on 28 day cycle

using Random, DifferentialEquations

include("../utilities/utils.jl")
include("rg_params.jl")

#time frames
hours = 24.0
end_time = 28.0*5.0
end_treat = 42.0

#dose amounts
avg_human_surface_area = 1.7 #m^2
tmz_treat_dose = 75.0*avg_human_surface_area
tmz_adjuv_dose = 150.0*avg_human_surface_area
dose_amount = 1920.0*avg_human_surface_area #1800
max_tested = dose_amount


tmz_treat_dosetimes = spaced_list(end_treat,1.0,0.0,0.0).*hours
tmz_adjuv_dosetimes = spaced_list(end_time,5.0,23.0,end_treat+28.0).*hours
rg_dosetimes = spaced_list(end_time-1.0,18.0,10.0,0.0).*hours #its the opposite18 days off and 10days on
const num_dose_times = length(rg_dosetimes)

doses = ones(length(rg_dosetimes)).*dose_amount

inject_times = sort(unique([rg_dosetimes;tmz_treat_dosetimes;tmz_adjuv_dosetimes]));

function adjust_dose(x)
    return min(dose_amount,relu(x))
end

# Convert arrays to sets for O(1) lookup
tmz_treat_dosetimes_set = Set(tmz_treat_dosetimes)
tmz_adjuv_dosetimes_set = Set(tmz_adjuv_dosetimes)
rg_dosetimes_set = Set(rg_dosetimes)

# Create a dictionary for rg_dosetimes and doses
# rg_dose_dict = Dict(zip(rg_dosetimes, doses))
map = Dict(zip(rg_dosetimes, Int64.(1:length(rg_dosetimes))))
function affect_dose!(integrator)
    SciMLBase.set_proposed_dt!(integrator, 0.1)
    
    current_time = integrator.t

    if current_time in tmz_treat_dosetimes_set
        integrator.u[states["AbsTMZ"]] += tmz_treat_dose
    elseif current_time in tmz_adjuv_dosetimes_set
        integrator.u[states["AbsTMZ"]] += tmz_adjuv_dose
    end
    if current_time in rg_dosetimes_set
        p = integrator.p
        # map = p[length(ode_params)]
        # doses = p[p_num+nb_scaling_params+1:end]
        doses = p[p_num+1:end]
        current_dose = doses[map[current_time]]
        integrator.u[states["AbsRG"]] += adjust_dose(current_dose)
    end
end
hit = PresetTimeCallback(inject_times, affect_dose!)

min_drug_dosage = 20.0*avg_human_surface_area*length(rg_dosetimes); #minimum drug dosage
min_tested = min_drug_dosage/length(rg_dosetimes)


function dose_heaviside(t, dosetimes, dose)
    if t in dosetimes
        return dose
    else
        return 0.0
    end
end