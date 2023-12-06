using DifferentialEquations, Plots

# Define the ODE system
function odesys!(du, u, p, t)
    a, b, c, d = p
    du[1] = -a * u[1]
    du[2] = a * u[1] - (b + d) * u[2] + c * u[3]
    du[3] = -c * u[3] + d * u[2]
end

# Parameters
params = [0.0168, 6.8, 7.363636363636362, 0.0216]
# Initial conditions
u0 = [3000.0, 0.0, 0.0]
# Time span
tspan = (0.0, 48.0)
# Solve the system numerically
prob = ODEProblem(odesys!, u0, tspan, params)
sol = solve(prob, Tsit5())

# Plot the solution for y(t)
plot(sol.t,sol[2,:], label="numerical", lw=2)
xlabel!("Time")
ylabel!("Values")
title!("Solution of the ODE System")


solution_string = read("scripts/prototype/analytic_sol.txt", String)
solution_string = replace(solution_string, "{" => "")
solution_string = replace(solution_string, "}" => "")
(x_rhs, y_rhs, z_rhs) = split(solution_string, ",")

function parse_math(solution_string)
    analyticSol = replace(solution_string, "E^" => "exp")
    analyticSol = replace(analyticSol, "Sqrt" => "sqrt")
    analyticSol = replace(analyticSol, "*^" => "*10^")
    analyticSol = replace(analyticSol, ".*" => "*")
    analyticSol = replace(analyticSol, "[t] -> " => "=")
    analyticSol = replace(analyticSol, "[" => "(")
    analyticSol = replace(analyticSol, "]" => ")")
    right_hand_side = split(analyticSol, "=")[2]
    function_body = "((t, p, u0) -> begin a, b, c, d = p; x0, y0, z0 = u0; $right_hand_side end)"
    rhs = eval(Meta.parse(function_body))
    return rhs
end

y_state = parse_math(y_rhs)

plot!(sol.t, [y_state(_t, params, u0) for _t in sol.t], label="analytical", lw=2, ls=:dash)


# Define your ODE system
function drug_pk!(du, u, p, t)
    Cl2, ka2, Vpla, Q, Vtis = p
    AbsRG, PlaRG, TisRG = u
    dAbsRG = -ka2 * AbsRG
    dPlaRG = ka2 * AbsRG - (Cl2 / Vpla) * PlaRG + (Q / Vtis) * TisRG - (Q / Vpla) * PlaRG
    dTisRG = -(Q / Vtis) * TisRG + (Q / Vpla) * PlaRG
    du .= [dAbsRG, dPlaRG, dTisRG]
end
# Initial conditions and parameters
u0 = [3000.0, 0.0, 0.0]  # [x0, y0, z0]
include("../../assets/pk_params.jl")
drug = "rg"
include("../../model/$(drug)_params.jl")
drug_params = eval(Symbol(uppercase(drug) * "_params"))
pk_param_names = collect(keys(drug_params))
pk_indices = indexin(pk_param_names, param_order)
pk_param_values= ode_params[pk_indices]
p_numeric = pk_param_values  # [a, b, c, d]
tspan = (0.0, 48.0)
# Assuming the same setup as before for the ODE problem
prob = ODEProblem(drug_pk!, u0, tspan, p_numeric)
sol = solve(prob)

# Plotting
plot!(sol.t, sol[2,:], label = "Simulation PK")
