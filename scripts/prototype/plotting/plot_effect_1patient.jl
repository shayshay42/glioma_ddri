drug = "rg"

include("../../../model/$(drug)_pkpd2.jl")
include("../../../model/$(drug)_dosing2.jl")
include("../../../model/$(drug)_params.jl")

println("this is the drug: ", drug)
println("this is the IC50 for TMZ: ", IC50_1, " converted to amount: ", IC50_1 * 140)
# Determine IC50 values
IC50_value = IC50_2 * 1000 * Vpla
println("this the IC50 for $drug: ", IC50_2, " converted to amount: ",IC50_value)

# include("../../../scripts/setup/init_integrate.jl")

function dose_affect!(integrator)
    # SciMLBase.set_proposed_dt!(integrator, 0.01)
    integrator.u[2] += integrator.u[end]
end
# cb = PeriodicCallback(dose_affect!, 24.0, initial_affect=true)

event_times = collect(1:18).*24.0
cb = PresetTimeCallback(event_times, dose_affect!)
tspan = (0.0, event_times[end]+(10.0*24.0))
    
p = ode_params
u0 = [17.7, 0.0, 0.0, 0.0, 5000.0]
prob = ODEProblem(simple_pkpd!, u0, tspan, p)
sol = solve(prob, Rodas4P2(; linsolve = nothing), callback=cb)

# auxillary functions computing for plotting
cPlaRG_sol = sol[3,:] ./ (1000 * Vpla)
exp2_sol = (cPlaRG_sol ./(psi*IC50_2)).^gamma_2
E_sol = (exp2_sol.*Imax_2) ./ (exp2_sol .+ 1)
t_sol = (-log.(log.(sol[1,:]./K) ./ log(C0/K)) ./ r) .+ 72
fun_sol = K .* (C0/K).^exp.(-r .* t_sol)
delta_sol = (E_sol .* fun_sol) ./ (72 .* sol[1,:])

# Get max values for each y-axis
max_y1 = 40 # Replace with actual max value of data associated with primary y-axis
# Set the same relative position for IC50 on both y-axes
relative_position_ic50 = IC50_value / max_y1
# Apply the relative position to secondary y-axis
ic50_secondary_axis = 0.5/relative_position_ic50
max_y2 = ic50_secondary_axis
p = plot(sol.t, sol[3,:], label="PlaRG", xlabel="Time", ylabel="Value", legend=:topleft, ylim=[0.0, max_y1])
plot!(sol.t, ones(length(sol.t)).*(IC50_2 * 1000 * Vpla), label="IC50", ls=:dash, alpha=0.5, color="blue")
p2 = twinx(p)
plot!(p2, sol.t, E_sol, label="Effect", color="red", ylabel="Effect", ylim=[0.0, max_y2],ls=:dash, legend=:topright, alpha=0.5)
plot!(p2, sol.t, ones(length(sol.t)).*(0.5), label="IC50", color="red", ls=:dash, alpha=0.5)
