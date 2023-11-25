struct ConditionalOutput
    loss::Float64
    ftv::Float64
    drug_auc::Float64
    tumor_auc::Float64
    trajectory::Array{Float64,2}
end

struct Patient
    idx::Int
    ode_parameters::Vector{Float64}
    minimal_tumor::Float64
    maximal_tumor::Float64
    optimal_doses::Vector{Float64}
    output_measures::Dict{String, ConditionalOutput}
end