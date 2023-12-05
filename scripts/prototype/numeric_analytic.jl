using Plots, DifferentialEquations

#(*Define the ODEs*)
# eqns = {x'[t] == -a x[t], y'[t] == a x[t] - (b + d) y[t] + c z[t], 
# z'[t] == -c z[t] + d y[t]};

# (*Define initial conditions*)
# initConds = {x[0] == x0, y[0] == y0, z[0] == z0};

# (*Solve the system of ODEs*)
# sol = DSolve[{eqns, initConds}, {x[t], y[t], z[t]}, t];

# (*Output the solution*)
# sol

# (*Export the solution to a \
# file*)Export["/Users/shayanhajhashemi/Documents/glioma_ddri/scripts/\
# prototype/analytic_sol.txt", sol]

function parse_mathematica_solution(sol_str)
    # Replace Mathematica sqrt and exp
    julia_str = replace(sol_str, "Sqrt" => "sqrt")
    julia_str = replace(julia_str, "E^" => "exp")  # Adding an opening parenthesis for 'exp'

    # Replace curly braces and other Mathematica-specific syntax
    julia_str = replace(julia_str, "{" => "(")
    julia_str = replace(julia_str, "}" => ")")
    julia_str = replace(julia_str, "[" => "(")
    julia_str = replace(julia_str, "]" => ")")
    julia_str = replace(julia_str, "->" => "=")
    
    return julia_str
end

function create_solution_functions(sol_str)
    x_solution_str = parse_mathematica_solution(match(r"x\[t\] -> ([^,]+),", sol_str).captures[1])
    y_solution_str = parse_mathematica_solution(match(r"y\[t\] -> ([^,]+),", sol_str).captures[1])
    z_solution_str = parse_mathematica_solution(match(r"z\[t\] -> ([^\}]+)", sol_str).captures[1])

    # Create Julia functions with extracted parameters and initial conditions
    x_func = eval(Meta.parse("((t, p, u0) -> begin a, b, c, d = p; x0, y0, z0 = u0; $(x_solution_str) end)"))
    y_func = eval(Meta.parse("((t, p, u0) -> begin a, b, c, d = p; x0, y0, z0 = u0; $(y_solution_str) end)"))
    z_func = eval(Meta.parse("((t, p, u0) -> begin a, b, c, d = p; x0, y0, z0 = u0; $(z_solution_str) end)"))

    return x_func, y_func, z_func
end

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

# Evaluate the analytical solutions
solution_string = read("scripts/prototype/analytic_sol.txt", String)
x_analytic, y_analytic, z_analytic = create_solution_functions(solution_string)

function convert_params(p)
    Cl2, ka2, Vpla, Q, Vtis = p
    return [ka2, Cl2/Vpla, Q/Vtis, Q/Vpla]
end
p_analytic = convert_params(p_numeric)
t_vals = range(tspan[1], tspan[2], length = length(sol.t))
y_analytic_vals = [y_analytic(t, p_analytic, u0) for t in sol.t]

# Plotting
plot(sol.t, sol[2,:], label = "Simulation (2nd State)", lw=3)
plot!(sol.t, y_analytic_vals, label = "Analytical (2nd State)", ls=:dash, lw=3)
