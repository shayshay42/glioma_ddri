from sympy import symbols, Eq, Function, dsolve
# from sympy.abc import a, b, c, d
k, c, v1, v2, q = symbols('k c v1 v2 q', real=True, positive=True)
a, b, c, d = symbols('ka b c d', real=True, positive=True)
x0, y0, z0 = symbols('x0 y0 z0', real=True, positive=True)

# Define the symbols and functions
t = symbols('t')
x = Function('x')(t)
y = Function('y')(t)
z = Function('z')(t)

# Define the equations
eq1 = Eq(x.diff(t), -k * x)
eq2 = Eq(y.diff(t), k * x - (c/v1) * y - (q/v1) * y + (q/v2) * z)
eq3 = Eq(z.diff(t), (q/v1) * y - (q/v2) * z)

# Solve the equations symbolically
sol_x = dsolve(eq1, ics={x.subs(t, 0): x0})
sol_z = dsolve(eq3, ics={z.subs(t, 0): z0})
sol_y = dsolve(eq2, ics={y.subs(t, 0): y0})

sol_x, sol_y, sol_z

# (Eq(x(t), C1*exp(-a*t)),
#  Eq(y(t), (C1 + Integral((a*x(t) + d*z(t))*exp(b*t)*exp(c*t), t))*exp(-t*(b + c))),
#  Eq(z(t), (C1 + c*Integral(y(t)*exp(d*t), t))*exp(-d*t)))

from sympy import latex

# Convert solutions to LaTeX
latex_x = latex(sol_x)
latex_y = latex(sol_y)
latex_z = latex(sol_z)

latex_x, latex_y, latex_z

