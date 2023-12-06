using BoundaryValueDiffEq, DifferentialEquations
using Plots, Makie
const g = 9.81
L = 1.0
tspan = (0.0, pi / 2)
function simplependulum!(du, u, p, t)
    θ = u[1]
    dθ = u[2]
    du[1] = dθ
    du[2] = -(g / L) * sin(θ)
    du
end

function bc1!(residual, u, p, t)
    residual[1] = u[end][1] - pi / 2 # the solution at the end of the time span should be pi/2
    residual[2] = u[end ÷ 2][1] + pi / 2 # the solution at the middle of the time span should be -pi/2
end
u0 = [pi/2, pi/2]
bvp1 = BVProblem(simplependulum!, bc1!, u0, tspan)
# sol1 = solve(bvp1, Shooting(Tsit5()), dt = 0.05)
sol2 = solve(bvp1, MIRK4(), dt = 0.05)


# Define finer ranges for θ and dθ
θ_range = range(-2π, stop=2π, length=40)
dθ_range = range(-10, stop=10, length=40)

# Initialize arrays for the vector field components
u_vals = zeros(length(θ_range), length(dθ_range))
v_vals = zeros(length(θ_range), length(dθ_range))

scale = 1/10
# Create the vector field plot
plot(size=(1500, 900))
for (θ, dθ) in Iterators.product(θ_range, dθ_range)
    du = zeros(2)
    simplependulum!(du, [θ, dθ], [], 0.0)
    norm = sqrt(du[1]^2 + du[2]^2)

    # Normalize and scale the vector
    u = du[1] / norm * scale
    v = du[2] / norm * scale

    # Add to the quiver plot
    # quiver!([θ], [dθ], quiver=([u], [v]), 
    #         linewidth=norm*scale, arrow=(:closed, norm*scale), color=:gray)
    arrows!([θ], [dθ], [u], [v], 
            linewidth=norm*scale, arrow=(:closed, norm*scale), color=:gray)
end

xlabel!("θ (rad)")
ylabel!("dθ/dt (rad/s)")
title!("Vector Field of a Simple Pendulum")

plot!(sol2[1,:], sol2[2,:], lw=5, label="Solution2")

sol1 = solve(bvp1, Shooting(Vern9()))#, dt = 0.05)

plot!(sol1[1,:], sol1[2,:], lw=5, label="Solution1")

prob = ODEProblem(simplependulum!, u0, tspan)
sol = solve(prob, Vern9())
plot!(sol[1,:], sol[2,:], lw=5, label="Solution")

function bc2a!(resid_a, u_a, p) # u_a is at the beginning of the time span
    resid_a[1] = u_a[1] + pi / 2 # the solution at the beginning of the time span should be -pi/2
end
function bc2b!(resid_b, u_b, p) # u_b is at the ending of the time span
    resid_b[1] = u_b[1] - pi / 2 # the solution at the end of the time span should be pi/2
end
bvp2 = TwoPointBVProblem(simplependulum!, (bc2a!, bc2b!), [pi / 2, pi / 2], tspan;
                         bcresid_prototype = (zeros(1), zeros(1)))
sol3 = solve(bvp2, MIRK4(), dt = 0.05)
plot!(sol3[1,:], sol3[2,:], lw=3, label="Solution3")

savefig("simple_pendulum.png")

plot(sol1.t, sol1[1,:], label="θ(t)")
plot!(sol1.t, sol1[2,:], label="dθ/dt(t)")
plot!(sol1.t, pi/2 .* ones(length(sol1.t)), label="pi/2")
plot!(sol1.t, -pi/2 .* ones(length(sol1.t)), label="-pi/2")

#plot the vector field
plot!([0, 0], [-pi/2, pi/2], label="")
plot!
