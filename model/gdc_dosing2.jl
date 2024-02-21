using Random, DifferentialEquations

include("../utilities/utils.jl")
include("gdc_params.jl")

hours = 24.0
end_time = 28.0*5.0
end_treat = 42.0

avg_human_surface_area = 1.7 #m^2
tmz_treat_dose = 75.0*avg_human_surface_area #360*mg_correcttion #0.075*1.7
tmz_adjuv_dose = 150.0*avg_human_surface_area #0.150*1.7
dose_amount = 600.0
const max_tested = dose_amount

tmz_treat_dosetimes = spaced_list(end_treat,1.0,0.0,0.0).*hours
tmz_adjuv_dosetimes = spaced_list(end_time,5.0,23.0,end_treat+28.0).*hours
gdc_dosetimes = spaced_list(end_time-1.0,18.0,10.0,0.0).*hours
const num_dose_times = length(rg_dosetimes)

doses = ones(length(gdc_dosetimes)).*dose_amount

inject_times = sort(unique([gdc_dosetimes;tmz_treat_dosetimes;tmz_adjuv_dosetimes]));

function adjust_dose(x)
    return min(dose_amount,relu(x))
end

# Convert arrays to sets for O(1) lookup
tmz_treat_dosetimes_set = Set(tmz_treat_dosetimes)
tmz_adjuv_dosetimes_set = Set(tmz_adjuv_dosetimes)
gdc_dosetimes_set = Set(gdc_dosetimes)

map = Dict(zip(gdc_dosetimes, Int64.(1:length(gdc_dosetimes))))
function affect_dose!(integrator)
    SciMLBase.set_proposed_dt!(integrator, 0.01)
    
    current_time = integrator.t

    if current_time in tmz_treat_dosetimes_set
        integrator.u[states["AbsTMZ"]] += tmz_treat_dose
    elseif current_time in tmz_adjuv_dosetimes_set
        integrator.u[states["AbsTMZ"]] += tmz_adjuv_dose
    end
    if current_time in gdc_dosetimes_set
        p = integrator.p
        # map = p[length(ode_params)]
        # doses = p[length(ode_params)+1:end]
        doses = p[p_num+1:end]
        current_dose = doses[map[current_time]]
        integrator.u[states["AbsGDC"]] += adjust_dose(current_dose)
    end
end
hit = PresetTimeCallback(inject_times, affect_dose!);

min_drug_dosage = 40.0*length(gdc_dosetimes); #minimum drug dosage
min_tested = min_drug_dosage/length(gdc_dosetimes); #minimum drug dose
# tmz_treat_dosetimes = spaced_list(end_treat,1,0,0).*hours
# function tmz_treat!(integrator)
#     SciMLBase.set_proposed_dt!( integrator, 0.01)
#     integrator.u[3] += tmz_treat_dose
# end
# tmz_treat_hit = PresetTimeCallback(tmz_treat_dosetimes, tmz_treat!);

# tmz_adjuv_dosetimes = spaced_list(end_time,5,23,end_treat+28).*hours
# function tmz_adjuv!(integrator)
#     SciMLBase.set_proposed_dt!( integrator, 0.01)
#     integrator.u[3] += tmz_adjuv_dose
# end
# tmz_adjuv_hit = PresetTimeCallback(tmz_adjuv_dosetimes, tmz_adjuv!);

# gdc_dosetimes = spaced_list(end_time-1,21,7,0).*hours
# doses = ones(length(gdc_dosetimes)).*dose_amount
# function gdc_dose!(integrator)
#     SciMLBase.set_proposed_dt!( integrator, 0.01)
#     hit_gdc = integrator.p[length(ode_params)+1:end][findall(x->x==integrator.t,gdc_dosetimes)][1]
#     integrator.u[6] += relu(hit_gdc)
# end
# hit_gdc = PresetTimeCallback(gdc_dosetimes, gdc_dose!);

# hit = CallbackSet(tmz_treat_hit, tmz_adjuv_hit, hit_gdc);

# tmz_treat_dosetimes = spaced_list(end_treat,1,0,0).*hours
# tmz_adjuv_dosetimes = spaced_list(end_time,5,23,end_treat+28).*hours
# gdc_dosetimes = spaced_list(end_time-1,21,7,0).*hours
# # println(gdc_dosetimes)
# inject_times = sort(unique([gdc_dosetimes;tmz_treat_dosetimes;tmz_adjuv_dosetimes;]));
# rng = Random.default_rng()
# Random.seed!(rng, 42)
# avg_huma_surface_area = 1.7 #m^2
# tmz_treat_dose = 75.0*avg_huma_surface_area #360*mg_correcttion #0.075*1.7
# tmz_adjuv_dose = 150.0*avg_huma_surface_area #0.150*1.7
# dose_amount = 400.0
# doses = ones(length(gdc_dosetimes)).*dose_amount

# function affect_dose!(integrator)
#     SciMLBase.set_proposed_dt!(integrator, 0.01)
#     ev = 0
#     if integrator.t in tmz_treat_dosetimes
#         integrator.u[3] += tmz_treat_dose
#         ev += 1
#     else 
#         nothing
#     end
#     if integrator.t in tmz_adjuv_dosetimes
#         integrator.u[3] += tmz_adjuv_dose
#         ev += 1
#     else 
#         nothing
#     end
#     if integrator.t in gdc_dosetimes
#         hit_gdc = integrator.p[length(ode_params[1,:])+1:end][findall(x->x==integrator.t,gdc_dosetimes)][1]
#         integrator.u[6] += relu(hit_gdc)
#         ev += 1
#     else 
#         nothing
#     end
#     if ev == 0
#         println("this should not get here!")
#     else 
#         nothing
#     end
# end

# hit = PresetTimeCallback(inject_times, affect_dose!);