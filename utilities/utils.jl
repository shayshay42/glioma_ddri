using Plots, LinearAlgebra, Statistics, ProgressMeter

"""
    trapezoidal_rule(t, y)

Calculate the area under a curve defined by `x` and `y` coordinates using the trapezoidal rule. This numerical integration method works by approximating the region under the curve as a series of trapezoids and then calculating the sum of their areas.

# Arguments
- `t::Vector`: A vector of `x` values (time points) in ascending order.
- `y::Vector`: A vector of `y` values (concentration levels) corresponding to each `x` value.

# Returns
- `Float`: The approximate area under the curve, as calculated by the trapezoidal rule.

# Example
```julia
julia> t = 0:0.1:1;
julia> y = t.^2;
julia> trapezoidal_rule(t, y)
0.33000000000000007
"""
function trapezoidal_rule(t, y)
  return 0.5*sum((t[2:end] - t[1:end-1]) .* (y[1:end-1] + y[2:end]))
end

"""
    spaced_list(p, n, m, b=1)

Create a list of integers where `n` numbers are sequentially appended 
followed by a jump of `m` numbers. This pattern is repeated until 
the end number `p` is reached or surpassed. 

Optionally, the start of the sequence can be adjusted from the 
default of `1` with the `b` parameter.

# Arguments
- `p::Integer`: The final number of the sequence. The function will stop 
  adding numbers to the list once this number is reached or surpassed.
- `n::Integer`: The number of sequential integers to append to the list 
  at a time.
- `m::Integer`: The number of integers to skip in the sequence after 
  each set of `n` integers is added.
- `b::Integer`: (optional) The beginning number of the sequence. Default 
  is `1`.

# Returns
- `Array{Integer}`: An array of integers that follows the specified 
  sequential pattern.

# Example
```julia
julia> spaced_list(20, 2, 3)
[1, 2, 6, 7, 11, 12, 16, 17]
"""
function spaced_list(p, n, m, b=1)
  # Initialize an empty list
  spaced_list = []
  # Initialize a counter variable
  counter = b
  # Use a while loop to iterate through the range of integers
  while counter <= p
    # Append `n` integers spaced by 1 to the list
    for i in 1:n
      spaced_list = [spaced_list; counter]
      counter += 1
      # Check if the counter has reached the end of the range
      if counter > p
        break
      end
    end
    # Add `m` to the counter to create the jump
    counter += m
  end
  return spaced_list
end

"""
    create_callback()

Create and return a callback function for an optimization routine and a list to store loss values.

The returned callback function prints the epoch number and loss value at each step of the optimization,
and also stores each loss value in the returned list. The callback is intended to be used with an 
optimization routine from the Optimization.jl package.

# Returns
- `track_progress::Function`: the callback function that takes `iter` (the current epoch number) and 
  `x` (the current parameter values) as arguments, computes the loss, and stores it in `loss_values`.
- `loss_values::Vector`: a list to hold the loss values computed during the optimization.

# Example
```julia
callback, loss_values = create_callback()
options = Optimization.Options(callback = callback, maxiters = 50)
res = Optimization.solve(optprob, opt, options)
@show res.u
@show loss_values

"""
struct StopException <: Exception end
function create_callback(total;verbose=true, animate=true, progress_bar=true, saving=true, early_stopping=true)
  loss_values = []
  iter = Ref(1)
  # current_res = []
  if animate
    anim = Animation()
  end
  # p = Progress(total,1)
  # pt = ProgressThresh(1.0, "Starting...") # create a ProgressThresh object
  function track_progress(x, current_loss)
      if animate || early_stopping
        # current_res = x
        push!(loss_values, current_loss)
      end

      if verbose
        println("Epoch: $iter, Loss: $current_loss")
      end
      #   if length(loss_values) > 2
      #     if abs(loss_values[end]-current_loss) < 1e-5
      #       throw(StopException())
      #     end
      #   end
      # end

      if progress_bar
        @info "Epoch: $iter, Loss: $current_loss"
        # println(x)
        # ProgressMeter.next!(p)
        # ProgressMeter.update!(pt, "Epoch: $iter, Loss: $current_loss")
      end

      if saving
        #deprecated
        iter_results[iter] = copy(res.u) # Store a copy of the result in the SharedArray
        serialize(open("./results/iteration_$(iter).jls", "w"), iter_results[iter]) # Serialize and save the results after each iteration
      end

      iter[] += 1
      #save loss in local list
      if animate
        # push!(loss_values, current_loss)
        frame(anim, plotter_gif(sim(x)))
      end

      if early_stopping
        window = 5
        if length(loss_values) > window+1
          avg_change = mean([abs(loss_values[end-i] - loss_values[end-(i+1)]) for i in 0:window])
          if avg_change < 1e-5
            println("Early stopping triggered!")
            return true
              # throw(StopException())
          end
        end
      end
      #continue the optimization
      return false  
  end
  if animate
    return track_progress, loss_values, iter, anim
  else
    return track_progress, loss_values, iter
  end
end

#short functions
logit(x::Real) = log(x / (one(x) - x))
logistic(x::Real) = inv(exp(-x) + one(x))

relu(x::Real) = max(zero(x), x)
erelu(x::Real) = max(eps(), x)

# function relu(x)
#   return max(0,x)
# end

# function erelu(x)
#   return max(eps(),x)
# end


function loss(θ, ode_params, scaling)
  drug_min_scaling = scaling[4]
  drug_max_scaling = scaling[3]
  prob_temp = remake(prob, p=[parameters..., scaling..., doses...])
  sol_temp = solve(prob_temp, Rodas4P2(), callback=hit, sensealg=nothing)
  cell = (sol_temp[end,end] - scaling[1])/(scaling[2]-scaling[1])
  drug = (sum(abs,doses)-drug_min_scaling)/(drug_max_scaling-drug_min_scaling)
  loss = cell + drug
  return loss
end