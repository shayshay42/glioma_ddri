function pk_pd!(du,u,p,t)#(du::Vector{Float64}, u::Vector{Float64}, p::Vector{Float64}, t::Float64)
    # gamma_1, psi, C0, D0, r, K, BW, IC50_1, Imax_1, IC50_2, gamma_2, Imax_2, xi, VD1, Cl1, k23, ka1, k32, Cl2, ka2, Vpla, Q, Vtis = p[1:length(ode_params)]
    gamma_1,psi,C0,D0,r,K,BW,IC50_1,Imax_1,IC50_2,gamma_2,Imax_2,xi,VD1,Cl1,k23,ka1,k32,Cl2,ka2,Vpla,Q,Vtis = p[1:length(ode_params)]
    # for name in keys(p)
    #     eval(:($name = nt.$name))
    # end
    xi = IC50_1/IC50_2
    C, D, AbsTMZ, PlaTMZ, CSFTMZ, AbsRG, PlaRG, TisRG, cAUC = u

    dAbsTMZ = -ka1 * AbsTMZ
    dPlaTMZ = ka1 * AbsTMZ - (Cl1 / VD1) * PlaTMZ - k23 * PlaTMZ + k32 * CSFTMZ
    dCSFTMZ = k23 * PlaTMZ - k32 * CSFTMZ

    dAbsRG = -ka2 * AbsRG
    dPlaRG = ka2 * AbsRG - (Cl2 / Vpla) * PlaRG + (Q / Vtis) * TisRG - (Q / Vpla) * PlaRG
    dTisRG = -(Q / Vtis) * TisRG + (Q / Vpla) * PlaRG
   
    cCSFTMZ = CSFTMZ / 140
    cPlaRG = PlaRG / (1000 * Vpla)

    #domain error
    cCSFTMZ = erelu(cCSFTMZ)
    cPlaRG = erelu(cPlaRG)
    C = erelu(C)

    pi1 = psi * IC50_1
    exp1 = (cCSFTMZ/pi1)^gamma_1
    exp2 = ((xi * cPlaRG)/pi1)^gamma_2
    exp12 = exp1 * exp2

    Imax_sum = Imax_1 + Imax_2
    E_num = Imax_1 * exp1 + Imax_2 * exp2 + Imax_sum * exp12 - Imax_1 * Imax_2 * exp12
    E_den = exp1 + exp2 + exp12 + 1
    E = E_num / E_den

    log_C0K = log(C0/K)
    t1 = -log(log(C/K) / log_C0K) / r
    t2 = t1 + 72
    fun = K * exp(log_C0K * exp(-r * t2))
    delta = (E * fun) / (72 * C)

    dC = C * r * log(K / C) - delta * C
    dD = delta * C

    dcAUC = C

    du .= [dC, dD, dAbsTMZ, dPlaTMZ, dCSFTMZ, dAbsRG, dPlaRG, dTisRG, dcAUC]#, dpRGauc]
end

function drug_pk!(du, u, p, t)
    Cl2, ka2, Vpla, Q, Vtis = p[1:end-1]
    AbsRG, PlaRG, TisRG = u
    dAbsRG = -ka2 * AbsRG
    dPlaRG = ka2 * AbsRG - (Cl2 / Vpla) * PlaRG + (Q / Vtis) * TisRG - (Q / Vpla) * PlaRG
    dTisRG = -(Q / Vtis) * TisRG + (Q / Vpla) * PlaRG
    du .= [dAbsRG, dPlaRG, dTisRG]
end

function drug_pk_dose!(du, u, p, t)
    Cl2, ka2, Vpla, Q, Vtis = p
    AbsRG, PlaRG, TisRG, dose = u
    dAbsRG = -ka2 * AbsRG
    dPlaRG = ka2 * AbsRG - (Cl2 / Vpla) * PlaRG + (Q / Vtis) * TisRG - (Q / Vpla) * PlaRG
    dTisRG = -(Q / Vtis) * TisRG + (Q / Vpla) * PlaRG
    du .= [dAbsRG, dPlaRG, dTisRG, 0.0]
end
