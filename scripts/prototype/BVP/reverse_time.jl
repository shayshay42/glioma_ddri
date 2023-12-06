using DifferentialEquations

# Define the differential equation
function differential_eq!(du, u, p, t)
    du[1] = -p[1] * u[1]  # Example of a simple first-order linear differential equation
end

# Parameters and final condition
p = [0.1]  # Parameter of the differential equation
final_condition = [1.0]  # Final condition

# Time span (going backwards in time)
tspan = (10.0, 0.0)  # Solving from t = 10 to t = 0

# Create a differential equation problem
prob = ODEProblem(differential_eq!, final_condition, tspan, p)

# Solve the problem
sol = solve(prob)

# Print the value at the beginning of the time span (which is the end in normal time direction)
println("The value at t = ", tspan[2], " is ", sol(tspan[2]))
