using OrderedCollections
# Define PK parameters
# GDC Parameters
GDC_params = OrderedDict(
    "ka2" => (mu = 0.0039, sigma = 0.017, mult = 60),
    "V2" => (mu = 0.077, sigma = 0.036, mult = 1000),
    "kel" => (mu = 0.029, sigma = 0.3, mult = 60),
    "k12" => (mu = 1.12, sigma = 0.064, mult = 60),
    "k21" => (mu = 2.86, sigma = 0.073, mult = 60)
)

# RG Parameters
RG_params = OrderedDict(
    "Cl2" => (mu = 0.00017, sigma = 0.73, mult = 1000*60),
    "ka2" => (mu = 0.00028, sigma = 0.25, mult = 60),
    "Vpla" => (mu = 0.0015, sigma = 1.76, mult = 1000),
    "Q" => (mu = 0.00000054, sigma = 2.25, mult = 1000*60),
    "Vtis" => (mu = 0.0000044, sigma = 5.08, mult = 1000)
)

# TMZ Parameters
TMZ_params = OrderedDict(
    "Cl1" => (mu = 10.0, sigma = 0.047, mult = 1.0),
    "k23" => (mu = 7.2 * 10^-4, sigma = 0.1649, mult = 1.0),
    "ka1" => (mu = 5.8, sigma = 0.8961, mult = 1.0)
)
# For RG_params_dict
# mu_ka2_RG = RG_params_dict["ka2"].mu