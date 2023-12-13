function pk_pd!(du,u,p,t)#(du::Vector{Float64}, u::Vector{Float64}, p::Vector{Float64}, t::Float64)
    # gamma_1, psi, C0, D0, r, K, BW, IC50_1, Imax_1, IC50_2, gamma_2, Imax_2, xi, VD1, Cl1, k23, ka1, k32, Cl2, ka2, Vpla, Q, Vtis = p[1:length(ode_params)]
    gamma_1,psi,C0,D0,r,K,BW,IC50_1,Imax_1,IC50_2,gamma_2,Imax_2,xi,VD1,Cl1,k23,ka1,k32,Cl2,ka2,Vpla,Q,Vtis = p[1:p_num]#length(ode_params)]
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

    # pi1 = psi * IC50_1
    exp1 = (cCSFTMZ /(psi*IC50_1))^gamma_1
    exp2 = (cPlaRG  /(psi*IC50_2))^gamma_2
    exp12 = exp1 * exp2

    Imax_sum = Imax_1 + Imax_2
    E_num = Imax_1 * exp1 + Imax_2 * exp2 + Imax_sum * exp12 - Imax_1 * Imax_2 * exp12
    E_den = exp1 + exp2 + exp12 + 1
    E = E_num / E_den

    t = (-log(log(C/K) / log(C0/K)) / r) + 72
    fun = K * (C0/K) ^ exp(-r * t)
    # delta = (E * fun) / (72 * C)

    dC = C * r * log(K / C) - (E * fun / 72)#delta * C
    dD = (E * fun / 72)#delta * C

    dcAUC = C
    du .= [dC, dD, dAbsTMZ, dPlaTMZ, dCSFTMZ, dAbsRG, dPlaRG, dTisRG, dcAUC]#, dpRGauc]
end
# create a dictionary of the state and their index
states = OrderedDict(zip(["C", "D", "AbsTMZ", "PlaTMZ", "CSFTMZ", "AbsRG", "PlaRG", "TisRG", "cAUC"], 1:9))

function pk_pd_alt!(du,u,p,t)#(du::Vector{Float64}, u::Vector{Float64}, p::Vector{Float64}, t::Float64)
    gamma_1,psi,C0,D0,r,K,BW,IC50_1,Imax_1,IC50_2,gamma_2,Imax_2,xi,VD1,Cl1,k23,ka1,k32,Cl2,ka2,Vpla,Q,Vtis = p[1:p_num]
    drug_doses = p[p_num+1:end]
    C, D, AbsTMZ, PlaTMZ, CSFTMZ, AbsRG, PlaRG, TisRG, cAUC = u

    # treatment_dose = dose_heaviside(t, tmz_treat_dosetimes_set, tmz_treat_dose)
    # adjuvant_dose = dose_heaviside(t, tmz_adjuv_dosetimes_set, tmz_adjuv_dose)
    # AbsTMZ = treatment_dose > 0 ? AbsTMZ + treatment_dose : (adjuvant_dose > 0 ? AbsTMZ + adjuvant_dose : AbsTMZ)    
    AbsTMZ = t in tmz_treat_dosetimes_set ? AbsTMZ + tmz_treat_dose : (t in tmz_adjuv_dosetimes_set ? AbsTMZ + tmz_adjuv_dose : AbsTMZ)
    dAbsTMZ = -ka1 * AbsTMZ
    dPlaTMZ = ka1 * AbsTMZ - (Cl1 / VD1) * PlaTMZ - k23 * PlaTMZ + k32 * CSFTMZ
    dCSFTMZ = k23 * PlaTMZ - k32 * CSFTMZ
    print("time is $t and the set its supposed to be in is $tmz_treat_dosetimes_set\n")
    AbsRG = t in rg_dosetimes_set ? AbsRG + drug_doses[map[t]] : AbsRG
    dAbsRG = -ka2 * AbsRG
    dPlaRG = ka2 * AbsRG - (Cl2 / Vpla) * PlaRG + (Q / Vtis) * TisRG - (Q / Vpla) * PlaRG
    dTisRG = -(Q / Vtis) * TisRG + (Q / Vpla) * PlaRG
   
    cCSFTMZ = CSFTMZ / 140
    cPlaRG = PlaRG / (1000 * Vpla)
    exp1 = (cCSFTMZ /(psi*IC50_1))^gamma_1
    exp2 = (cPlaRG  /(psi*IC50_2))^gamma_2
    exp12 = exp1 * exp2

    Imax_sum = Imax_1 + Imax_2
    E_num = Imax_1 * exp1 + Imax_2 * exp2 + Imax_sum * exp12 - Imax_1 * Imax_2 * exp12
    E_den = exp1 + exp2 + exp12 + 1
    E = E_num / E_den

    t = (-log(log(C/K) / log(C0/K)) / r) + 72
    fun = K * (C0/K) ^ exp(-r * t)
    dC = C * r * log(K / C) - (E * fun / 72)#delta * C
    dD = (E * fun / 72)#delta * C

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

function simple_pkpd!(du, u, p, t)
    gamma_1,psi,C0,D0,r,K,BW,IC50_1,Imax_1,IC50_2,gamma_2,Imax_2,xi,VD1,Cl1,k23,ka1,k32,Cl2,ka2,Vpla,Q,Vtis = p
    C, AbsRG, PlaRG, TisRG, dose = u
    C = erelu(C)

    cPlaRG = PlaRG / (1000 * Vpla)
    cPlaRG = erelu(cPlaRG)
    exp2 = (cPlaRG /(psi*IC50_2))^gamma_2
    E = Imax_2 * exp2 / (exp2 + 1)

    t = (-log(log(C/K) / log(C0/K)) / r) + 72
    fun = K * (C0/K)^exp(-r * t)
    # delta = (E * fun) / (72 * C)

    # E = 1.0
    # fun = 1.0

    dC = C*r*log(K/C) - (E*fun)/72
    dAbsRG = -ka2 * AbsRG
    dPlaRG = ka2 * AbsRG - (Cl2 / Vpla) * PlaRG + (Q / Vtis) * TisRG - (Q / Vpla) * PlaRG
    dTisRG = -(Q / Vtis) * TisRG + (Q / Vpla) * PlaRG

    du .= [dC, dAbsRG, dPlaRG, dTisRG, 0.0]
end

