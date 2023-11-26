# using Symbolics
# using Latexify

# # Define the symbols
# @variables t a b c d x(t) y(t) z(t)
# D = Differential(t)

# # System of equations
# eqs = [
#     D(x) ~ -a*x,
#     D(y) ~ a*x - b*y - c*y + d*z,
#     D(z) ~ c*y - d*z
# ]

# # Attempt to solve the system symbolically
# sol = Symbolics.solve_for(eqs, [x, y, z])

# println(sol)
# latexify(sol)

using SymPy

# Define symbols
@vars t a b c d x0 y0 z0
x = SymFunction("x")
y = SymFunction("y")
z = SymFunction("z")

# Define the equations
eq1 = Eq(diff(x(t), t), -a*x(t))
eq2 = Eq(diff(y(t), t), a*x(t) - b*y(t) - c*y(t) + d*z(t))
eq3 = Eq(diff(z(t), t), c*y(t) - d*z(t))

# Solve the equations with initial conditions
sol_x = dsolve(eq1, x(t), ics={x(0): x0})
sol_y = dsolve(eq2, y(t), ics={y(0): y0})
sol_z = dsolve(eq3, z(t), ics={z(0): z0})

# Display the solutions
(sol_x, sol_y, sol_z)
