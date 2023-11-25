using Distributions
using Plots

function confidence_interval(conf, param, mu, sigma, multiplier)
    
    # Initialize arrays for lower and upper bounds
    lower_bounds = zeros(length(mu))
    upper_bounds = zeros(length(mu))

    # Compute the quantiles for the lower and upper bounds of the confidence interval
    alpha = (1.0 - conf) / 2.0
    lower_quantile = alpha
    upper_quantile = 1.0 - alpha

    # Loop over all mu, sigma, and param values and calculate bounds
    for i in 1:length(mu)
        if param[i] == "Q"
            # For the logit distribution
            logit_dist = LogitNormal(logit(mu[i]), sigma[i])
            lower_bounds[i] = quantile(logit_dist, lower_quantile) * multiplier[i]
            upper_bounds[i] = quantile(logit_dist, upper_quantile) * multiplier[i]
        else
            # For the log-normal distribution
            log_dist = LogNormal(log(mu[i])+log(multiplier[i]), sigma[i])
            lower_bounds[i] = quantile(log_dist, lower_quantile) #* multiplier[i]
            upper_bounds[i] = quantile(log_dist, upper_quantile) #* multiplier[i]
        end
    end

    return lower_bounds, upper_bounds
end

include("../../assets/pk_params.jl")
# # Define PK parameters
# GDC_params = ["ka2","V2" ,"kel","k12","k21"]
# GDC_mu =     [0.0039, 0.077, 0.029, 1.12, 2.86]
# GDC_sigma =  [0.017, 0.036, 0.3, 0.064, 0.073]
# GDC_mult =   [60, 1000, 60, 60, 60]
 
# RG_params =  ["ka2", "Cl2", "Vpla", "Q", "Vtis"]
# RG_mu =      [0.00028, 0.00017, 0.0015, 0.00000054, 0.0000044]
# RG_sigma =   [0.25, 0.73, 1.76, 2.25, 5.08]
# RG_mult =    [60, 1000*60, 1000, 1000*60, 1000]

# TMZ_params = ["Cl1", "k23", "ka1"]
# TMZ_mu =     [10.0, 7.2 * 10^-4, 5.8]
# TMZ_sigma =  [0.047, 0.1649, 0.8961]
# TMZ_mult =   [1.0, 1.0, 1.0]

# Define confidence level
confidence_level = 0.9999 # Change this to 0.68 for one standard deviation

# Calculate the confidence intervals for GDC parameters
GDC_lower, GDC_upper = confidence_interval(confidence_level, GDC_params, GDC_mu, GDC_sigma, GDC_mult)

RG_lower, RG_upper = confidence_interval(confidence_level, RG_params, RG_mu, RG_sigma, RG_mult)

TMZ_lower, TMZ_upper = confidence_interval(confidence_level, TMZ_params, TMZ_mu, TMZ_sigma, TMZ_mult)

# using Printf

# function format_number(num)
#     if abs(num) < 0.001 || abs(num) > 9999
#         return @sprintf("%.1e", num)  # Use scientific notation
#     else
#         return @sprintf("%.3f", num)  # Use fixed decimal notation
#     end
# end

# function print_table()
#     datasets = [
#         ("GDC", GDC_params, GDC_mu, GDC_sigma, GDC_mult, GDC_lower, GDC_upper),
#         ("RG", RG_params, RG_mu, RG_sigma, RG_mult, RG_lower, RG_upper),
#         ("TMZ", TMZ_params, TMZ_mu, TMZ_sigma, TMZ_mult, TMZ_lower, TMZ_upper)
#     ]

#     println("+---------+-----------+---------+--------+------+----------------+---------------+")
#     println("| Dataset | Parameter |   Mu    | Sigma  | Mult | Lower Bound    | Upper Bound   |")
#     println("+---------+-----------+---------+--------+------+----------------+---------------+")

#     for dataset in datasets
#         name, params, mus, sigmas, mults, lowers, uppers = dataset
#         for i in 1:length(params)
#             println("| $(rpad(name, 7)) | $(rpad(params[i], 9)) | $(rpad(format_number(mus[i]), 7)) | $(rpad(format_number(sigmas[i]), 6)) | $(rpad(mults[i], 4)) | $(rpad(format_number(lowers[i]), 14)) | $(rpad(format_number(uppers[i]), 13)) |")
#         end
#     end

#     println("+---------+-----------+---------+--------+------+----------------+---------------+")
# end

# print_table()

# +---------+-----------+---------+--------+------+----------------+---------------+
# | Dataset | Parameter |   Mu    | Sigma  | Mult | Lower Bound    | Upper Bound   |
# +---------+-----------+---------+--------+------+----------------+---------------+
# | GDC     | ka2       | 0.004   | 0.017  | 60   | 0.219          | 0.250         |
# | GDC     | V2        | 0.077   | 0.036  | 1000 | 66.936         | 88.577        |
# | GDC     | kel       | 0.029   | 0.300  | 60   | 0.542          | 5.590         |
# | GDC     | k12       | 1.120   | 0.064  | 60   | 52.388         | 86.200        |
# | GDC     | k21       | 2.860   | 0.073  | 60   | 129.173        | 227.962       |
# | RG      | ka2       | 2.8e-04 | 0.250  | 60   | 0.006          | 0.044         |
# | RG      | Cl2       | 1.7e-04 | 0.730  | 60000| 0.596          | 174.604       |
# | RG      | Vpla      | 0.002   | 1.760  | 1000 | 0.002          | 1412.204      |
# | RG      | Q         | 5.4e-07 | 2.250  | 60000| 5.1e-06        | 204.552       |
# | RG      | Vtis      | 4.4e-06 | 5.080  | 1000 | 1.1e-11        | 1.7e+06       |
# | TMZ     | Cl1       | 10.000  | 0.047  | 1.0  | 8.329          | 12.006        |
# | TMZ     | k23       | 7.2e-04 | 0.165  | 1.0  | 3.8e-04        | 0.001         |
# | TMZ     | ka1       | 5.800   | 0.896  | 1.0  | 0.178          | 189.467       |
# +---------+-----------+---------+--------+------+----------------+---------------+
