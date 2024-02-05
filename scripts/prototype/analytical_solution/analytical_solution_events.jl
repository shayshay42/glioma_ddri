solution_string = read("scripts/prototype/analytical_solution/analytic_sol.txt", String)
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
    println(function_body)
    return rhs
end

x_state = parse_math(x_rhs)
y_state = parse_math(y_rhs)
z_state = parse_math(z_rhs)

event_times = collect(1:18).*24.0
tspan = (0.0, event_times[end]+(10.0*24.0))
u0 = [3000.0, 0.0, 0.0]  # [x0, y0, z0]
include("../../../assets/pk_params.jl")
drug = "rg"
include("../../../model/$(drug)_params.jl")
drug_params = eval(Symbol(uppercase(drug) * "_params"))
pk_param_names = collect(keys(drug_params))
pk_indices = indexin(pk_param_names, param_order)
pk_param_values= ode_params[pk_indices]
p_numeric = pk_param_values
function convert_params(p)
    Cl2, ka2, Vpla, Q, Vtis = p
    return [ka2, Cl2/Vpla, Q/Vtis, Q/Vpla]
end
p_analytic = convert_params(p_numeric)

doses = ones(length(event_times)).*3000.0
# start_state = u0
# for (i, (t, dose)) in enumerate(zip(event_times, doses))
#     state_before_dose = [x_state(t, p_analytic, start_state), y_state(t, p_analytic, start_state), z_state(t, p_analytic, start_state)]
#     start_state = [state_before_dose[1]+dose, state_before_dose[2], state_before_dose[3]]
# end

using Plots

# Assuming functions and parameters are already defined

# Time span and time step
total_time = tspan[2]
time_step = 0.1  # Adjust time step as needed for desired resolution

# Create a continuous time grid
time_grid = 0:time_step:total_time

# Initialize arrays for storing time points and y_state values
times = []
y_values = []

# Start with the initial state
start_state = u0

# Time since last reset (initial condition or event)
elapsed_time_since_reset = 0.0
last_event_time = 0.0  # Time of the last event

# Iterate over the time grid
for t in time_grid
    # Update elapsed time since last reset
    elapsed_time_since_reset = t - last_event_time

    # Check if the current time is an event time and apply dose if it is
    if t in event_times
        dose = doses[findfirst(isequal(t), event_times)]  # Get the corresponding dose
        state_before_dose = [x_state(elapsed_time_since_reset, p_analytic, start_state), y_state(elapsed_time_since_reset, p_analytic, start_state), z_state(elapsed_time_since_reset, p_analytic, start_state)]
        start_state = [state_before_dose[1] + dose, state_before_dose[2], state_before_dose[3]]
        last_event_time = t
        elapsed_time_since_reset = 0.0  # Reset the elapsed time
    end

    # Calculate state without dose
    current_state = [x_state(elapsed_time_since_reset, p_analytic, start_state), y_state(elapsed_time_since_reset, p_analytic, start_state), z_state(elapsed_time_since_reset, p_analytic, start_state)]

    # Store the time and y_state value
    push!(times, t)
    push!(y_values, current_state[2])
end

# Plot the y_state values over time
plot(times./24, y_values, title="RG PK analytical solution Plasma Compartment", xlabel="Time (days)", ylabel="Amount (mg)", label="RG Plasma Compartment")
plot!(times./24, ones(length(times)).*(1.0*(IC50_2*(Vpla*1000))), label="IC50 = $(1.0*(IC50_2*(Vpla*1000)))")