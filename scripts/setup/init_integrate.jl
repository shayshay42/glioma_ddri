using DifferentialEquations, ModelingToolkit, ComponentArrays

u0 = zeros(length(states))
u0[1] = C0
tspan = (0,end_time+7).*hours
# p = [ode_params..., default_scaling..., doses...]
# p = ComponentArray(θ=ode_params, doses=doses)
p = [ode_params..., doses...]
# p = [ode_params; default_scaling; doses]
prob = ODEProblem(pk_pd!,u0,tspan,p)
sys = modelingtoolkitize(prob)
sys = structural_simplify(sys)
func = ODEFunction(sys, jac=true)
prob = ODEProblem(func, u0, tspan, p)