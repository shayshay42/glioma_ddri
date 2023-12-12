using DifferentialEquations, ModelingToolkit

u0 = zeros(9)
u0[1] = 17.7
tspan = (0,end_time+7).*hours
p = [ode_params; 1.0; 1.0;doses]
prob = ODEProblem(pk_pd!,u0,tspan,p)
sys = modelingtoolkitize(prob)
sys = structural_simplify(sys)
func = ODEFunction(sys, jac=true)
prob_jac = ODEProblem(func, u0, tspan, p)