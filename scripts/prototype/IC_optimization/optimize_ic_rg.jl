drug = "rg"

include("../../model/$(drug)_params.jl")

println("this is the drug: ", drug)
println("this is the IC50 for TMZ: ", IC50_1)
println("this the IC50 for $drug: ", IC50_2)

using DifferentialEquations, ModelingToolkit

function rg_pk(du, u, p, t)
    Cl2, ka2, Vpla, Q, Vtis = p[1:end-1]
    AbsRG, PlaRG, TisRG = u
    dAbsRG = -ka2 * AbsRG
    dPlaRG = ka2 * AbsRG - (Cl2 / Vpla) * PlaRG + (Q / Vtis) * TisRG - (Q / Vpla) * PlaRG
    dTisRG = -(Q / Vtis) * TisRG + (Q / Vpla) * PlaRG
    du .= [dAbsRG, dPlaRG, dTisRG]
end

dose_amount = 3000.0
function dose!(integrator)
    # SciMLBase.set_proposed_dt!(integrator, 0.01)
    integrator.u[1] += integrator.p[end]
end
cb = PeriodicCallback(dose!, 24.0, initial_affect=true)

p = [Cl2, ka2, Vpla, Q, Vtis, dose_amount]
u0 = [0.0, 0.0, 0.0]
tspan = (0.0, 18.0).*24.0
prob = ODEProblem(rg_pk, u0, tspan, p)
# sys = modelingtoolkitize(prob)
# sys = structural_simplify(sys)
sol = solve(prob, callback=cb)

using Plots
# plot(sol.t, sol[2,:])
# plot!(sol[2,:])
# # plot(sol, vars=(0, 1), label="AbsRG")
# plot(sol, vars=(0, 2), label="PlaRG")
# plot!(sol.t, IC50_2*(Vpla*1000).*ones(length(sol.t)), label="IC50_2")
# plot(sol, vars=(0, 3), label="TisRG")
using Optim

mults = [0.01, 0.1, 0.5, 1.0, 1.5, 2.0]
function objective(dose, mult)
    # final_cond = [guess[1], 0.009, guess[2]]  # AbsRG and TisRG are the guesses
    # prob = ODEProblem(rg_pk, final_cond, tspan, params)
    # sol = solve(prob, Tsit5(), reverse=true)
    tmp_prob = remake(prob, p=[Cl2, ka2, Vpla, Q, Vtis, dose[1]])
    tmp_sol = solve(tmp_prob, callback=cb)
    # Analyze the solution to find the peak of PlaRG
    peak_PlaRG = maximum(tmp_sol[2, :])  # Assuming PlaRG is the second variable
    return abs2(peak_PlaRG - mult*(IC50_2*(Vpla*1000)))  # Objective is to minimize this difference
end
avg_human_surface_area=1.7
# Initial guess for AbsRG and TisRG
initial_guess = [20.0, 100.0, 500.0, 1000.0, 2000.0, 3000.0].*avg_human_surface_area # Adjust this as needed
i=1
println("Initial guess for multiplier $(mults[i]): $(initial_guess[i])")
objective([initial_guess[i]], mults[i])

result = optimize(x-> objective(x, mults[i]), [initial_guess[i]])
optimal_dose = Optim.minimizer(result)
println("Optimal dose for multiplier $(mults[i]): $optimal_dose")

# Optimization for each multiplier
for (i,mult) in enumerate(mults)
    result = optimize(x -> objective(x, mult), [initial_guess[i]], LBFGS())
    optimal_dose = Optim.minimizer(result)
    println("Optimal dose for multiplier $mult: $optimal_dose, reached loss: $(result.minimum)")
end

# Perform optimization
opt_result = optimize(objective_function, initial_guess)

# Use the result of optimization to get the final condition
optimized_guess = Optim.minimizer(opt_result)
final_cond = [optimized_guess[1], 0.009, optimized_guess[2]]

# Now solve the ODE with this condition
prob = ODEProblem(rg_pk, final_cond, tspan, params)
sol = solve(prob, Tsit5(), reverse=true)

# Initial condition at the start of the process (end of the backward solution)
initial_cond = sol[end]
