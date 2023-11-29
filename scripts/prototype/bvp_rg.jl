using Pkg
Pkg.activate(".")

using DifferentialEquations, Plots
using BoundaryValueDiffEq

function drug_pk!(du, u, p, t)
    Cl2, ka2, Vpla, Q, Vtis = p
    AbsRG, PlaRG, TisRG = u
    dAbsRG = -ka2 * AbsRG
    dPlaRG = ka2 * AbsRG - (Cl2 / Vpla) * PlaRG + (Q / Vtis) * TisRG - (Q / Vpla) * PlaRG
    dTisRG = -(Q / Vtis) * TisRG + (Q / Vpla) * PlaRG
    du .= [dAbsRG, dPlaRG, dTisRG]
end

# u0 = [x, 0.0, 0.0]
# ut = [?, ic50, ?] where t=?

include("../../model/rg_params.jl")
p = [Cl2, ka2, Vpla, Q, Vtis]
tspan = (0.0, 24.0)

# Initial guess for the solution
u0_guess = [5000.0, 0.0, 0.0] # Adjust as needed

# Boundary condition function
function bc_drug_pk!(residual, sol, p, t)
    residual[1] = maximum(sol[2,:]) - (IC50_2 * (1000 * Vpla)) # Ensure the second state matches ic50 amount conversion at end time
    residual[2] = sol[1][2] # Ensure the seocnd state is 0 at the beginning
    residual[3] = sol[1][3] # Ensure the third state is 0 at the beginning
end

# Define and solve the boundary value problem
bvp_drug_pk = BVProblem(drug_pk!, bc_drug_pk!, u0_guess, tspan, p)
sol_drug_pk = solve(bvp_drug_pk, Shooting(KenCarp4()))

# Plotting the solution (optional)
plot(sol_drug_pk, title="Solution of the Pharmacokinetics BVP", xlabel="Time", ylabel="States")
plot(sol_drug_pk, vars=(0, 2), label="PlaRG")
plot!(sol_drug_pk.t, IC50_2*(Vpla*1000).*ones(length(sol_drug_pk.t)), label="IC50_2")
plot(sol_drug_pk, vars=(0, 3), label="TisRG")
println("The peak of PlaRG is: ", maximum(sol_drug_pk[2, :]), " compared to IC50_2*(Vpla*1000) = ", IC50_2*(Vpla*1000))
println("The initial condition for AbsRG is: ", sol_drug_pk[1, 1])
println("The initial condition for PlaRG is: ", sol_drug_pk[2, 1])
println("The initial condition for TisRG is: ", sol_drug_pk[3, 1])
#actually need a two point BVP problem