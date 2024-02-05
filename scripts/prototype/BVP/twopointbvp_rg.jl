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
tspan = (0.0, 3*24.0)

# Initial guess for the solution
u_guess = [5000.0, 0.0, 0.0]

# Boundary condition at the start of the time span
function bc_start!(resid_start, u_start, p)
    resid_start[1] = u_start[2] # State 2 should be 0 at the start
    resid_start[2] = u_start[3] # State 3 should be 0 at the start
end

# Boundary condition at the end of the time span
# Here, you can specify conditions for the final state or leave it free
function bc_end!(resid_end, u_end, p)
    resid_end[1] = u_end[2] - (IC50_2 * (1000 * Vpla)) # Ensure the second state matches ic50 amount conversion at end time
end

# Define the structure of the boundary residuals
bcresid_prototype = (zeros(3), zeros(3)) # Tuple of two zero arrays each of length 3

# Define the BVP
bvp_drug_pk = TwoPointBVProblem(drug_pk!, (bc_start!, bc_end!), u_guess, tspan, p, bcresid_prototype=bcresid_prototype)

# Solve the BVP
sol_drug_pk = solve(bvp_drug_pk, Shooting(KenCarp4()))#, dt=0.05)

# Plotting the solution (optional)
plot(sol_drug_pk.t, sol_drug_pk[1,:], title="Solution of the Pharmacokinetics BVP", xlabel="Time", ylabel="States")
plot(sol_drug_pk.t, sol_drug_pk[2,:], label="PlaRG")
