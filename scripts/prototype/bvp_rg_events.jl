using Pkg
Pkg.activate(".")

using DifferentialEquations, Plots
using BoundaryValueDiffEq

function drug_pk!(du, u, p, t)
    Cl2, ka2, Vpla, Q, Vtis = p
    AbsRG, PlaRG, TisRG, dose = u
    dAbsRG = -ka2 * AbsRG
    dPlaRG = ka2 * AbsRG - (Cl2 / Vpla) * PlaRG + (Q / Vtis) * TisRG - (Q / Vpla) * PlaRG
    dTisRG = -(Q / Vtis) * TisRG + (Q / Vpla) * PlaRG
    du .= [dAbsRG, dPlaRG, dTisRG, 0.0]
end

include("../../model/rg_params.jl")
p = [Cl2, ka2, Vpla, Q, Vtis]

# Initial guess for the solution
global u0_guess = [0.0, 0.0, 0.0, 5000.0] # Adjust as needed

# Boundary condition function
function bc_drug_pk_1!(residual, sol, p, t)
    residual[1] = maximum(sol[2,:]) - (IC50_2 * (1000 * Vpla)) # Ensure the second state matches ic50 amount conversion at end time
    residual[2] = sol[1][2] # Ensure the seocnd state is 0 at the beginning
    residual[3] = sol[1][3] # Ensure the third state is 0 at the beginning
    residual[4] = sol[1][1]
end

function dose_affect!(integrator)
    # SciMLBase.set_proposed_dt!(integrator, 0.01)
    integrator.u[1] += integrator.u[4]
end
event_times = collect(1:18).*24.0
cb = PresetTimeCallback(event_times, dose_affect!)
tspan = (0.0, event_times[end]+(10.0*24.0))

bvp_drug_pk = BVProblem(drug_pk!, bc_drug_pk!, u0_guess, tspan, p)
sol_drug_pk = solve(bvp_drug_pk, Shooting(KenCarp4()), callback=cb)

# Define and solve the boundary value problem
prob = ODEProblem(drug_pk!, u0_guess, tspan, p)
sol = solve(prob, callback=cb)
plot(sol.t, sol[1,:], label="AbsRG")
plot(sol.t, sol[2,:], label="PlaRG")
println(count)


# Plotting the solution (optional)
# plot(sol_drug_pk, title="Solution of the Pharmacokinetics BVP", xlabel="Time", ylabel="States")
plot(sol_drug_pk.t, sol_drug_pk[2,:], label="PlaRG")
plot!(sol_drug_pk.t, IC50_2*(Vpla*1000).*ones(length(sol_drug_pk.t)), label="IC50_2")

# plot(sol_drug_pk, vars=(0, 3), label="TisRG")
println("The peak of PlaRG is: ", maximum(sol_drug_pk[2, :]), " compared to IC50_2*(Vpla*1000) = ", IC50_2*(Vpla*1000))
println("The dose amount is: ", sol_drug_pk[4, 1])
# println("The initial condition for PlaRG is: ", sol_drug_pk[2, 1])
# println("The initial condition for TisRG is: ", sol_drug_pk[3, 1])
#actually need a two point BVP problem