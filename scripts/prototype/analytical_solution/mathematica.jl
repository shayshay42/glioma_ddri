# using MathLink

# # Open a link to a Mathematica session
# wmlink = W""

# # Define the ODE system as a string in Mathematica syntax
# ode_system = """
#     DSolve[{{
#         x'[t] == -a x[t],
#         y'[t] == a x[t] - b y[t] + c z[t],
#         z'[t] == -c z[t] + d y[t]},
#         {x[0] == x0, y[0] == y0, z[0] == z0}},
#         {x, y, z}, t]
# """

# # Send the command to Mathematica and solve the ODE system
# weval(ode_system)

# # Retrieve the solution from Mathematica
# solution = MathLink.get(wmlink)

# # Close the link
# MathLink.close(wmlink)

# # Process the solution as needed in Julia
# # Depending on the format of the solution, you may need to parse or convert it to a usable form

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

    # Create Julia functions
    x_func = eval(Meta.parse("((t, a, b, c, d, x0, y0, z0) -> $x_solution_str)"))
    y_func = eval(Meta.parse("((t, a, b, c, d, x0, y0, z0) -> $y_solution_str)"))
    z_func = eval(Meta.parse("((t, a, b, c, d, x0, y0, z0) -> $z_solution_str)"))

    return x_func, y_func, z_func
end

# Replace this string with your actual solution string
solution_string = read("scripts/prototype/analytic_sol.txt", String)
x, y, z = create_solution_functions(solution_string)

# Example of using the created functions
println(x(1, 1, 1, 1, 1, 1, 1, 1))
println(y(1, 1, 1, 1, 1, 1, 1, 1))
println(z(1, 1, 1, 1, 1, 1, 1, 1))
