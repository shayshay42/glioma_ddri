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

analyticSol = "{x[t] -> (1.26318548017794*^-21*E^(6.6272300439191065*t) - 4.031722685188401*^-20*E^(7.591606319717258*t) + 3000.*E^(14.185236363636363*t))/E^(14.202036363636363*t), y[t] -> (-1.460188411332483*(1.*E^(6.6272300439191065*t) + 4.088504690657046*E^(7.591606319717258*t) - 5.088504690657046*E^(14.185236363636363*t)))/E^(14.202036363636363*t), z[t] -> (0.14935869794233111*(1.*E^(6.6272300439191065*t) - 1.1462587784535356*E^(7.591606319717258*t) + 0.14625877845353558*E^(14.185236363636363*t)))/E^(14.202036363636363*t)}"
#parse the solution string into 3 function x , y, z
analyticSol = replace(analyticSol, "{" => "")
analyticSol = replace(analyticSol, "}" => "")
analyticSol = replace(analyticSol, "E^" => "exp")
analyticSol = replace(analyticSol, "*^" => "*10^")
analyticSol = replace(analyticSol, ".*" => "*")
analyticSol = replace(analyticSol, "[t] -> " => "=")

(x_rhs, y_rhs, z_rhs) = split(analyticSol, ",")


x_state = eval(Meta.parse("t -> $(split(x_rhs, '=')[2])"))
y_state = eval(Meta.parse("t -> $(split(y_rhs, '=')[2])"))
z_state = eval(Meta.parse("t -> $(split(z_rhs, '=')[2])"))

plot!(sol.t, y_state.(sol.t), label="analytical", lw=2, ls=:dash)