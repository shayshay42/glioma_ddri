include("../../model/$(drug)_pkpd2.jl")
include("../../model/$(drug)_dosing2.jl")
include("../../model/$(drug)_params.jl")

u0 = zeros(9)
u0[1] = 17.7
tspan = (0,end_time+7).*hours
p = [ode_params; 1.0; 1.0;doses]
prob = ODEProblem(pk_pd!,u0,tspan,p)
sys = modelingtoolkitize(prob)
sys = structural_simplify(sys)
prob_jac = ODEProblem(sys, u0, tspan, p)