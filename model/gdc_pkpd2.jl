function pk_pd!(du, u, p, t)
    gamma_1, psi, C0, D0, r, K, BW, IC50_1, Imax_1, IC50_2, gamma_2, Imax_2, xi, VD1, Cl1, k23, ka1, k32, ka2, V2, kel, k12, k21 = p[1:length(ode_params)]
    xi = IC50_1 / IC50_2
    C, D, AbsTMZ, PlaTMZ, CSFTMZ, AbsGDC, PlaGDC, PeriphGDC, cAUC = u
    
    cCSFTMZ = CSFTMZ / 140
    cPlaGDC = PlaGDC / (1000 * V2)

    #inplace if negative make it = postiive epsilon
    cCSFTMZ = erelu(cCSFTMZ)
    cPlaGDC = erelu(cPlaGDC)
    C = erelu(C)

    pi1 = psi*IC50_1
    exp1 = (cCSFTMZ/pi1)^gamma_1
    exp2 = ((xi*cPlaGDC)/pi1)^gamma_2
    exp12 = exp1*exp2

    Imax_sum = Imax_1+Imax_2
    E_num =(Imax_1*exp1)+(Imax_2*exp2)+(Imax_sum*exp12)-(Imax_1*Imax_2*exp12)
    E_den = exp1+exp2+exp12+1
    E = E_num/E_den

    lck = log(C0/K)
    t1=-log(log(C/K)/lck)/r
    t2=t1+72
    fun = K*exp(lck*exp(-r*t2))
    delta=(E/72)*(fun/C)

    dC = (C*r*log(K/C))-delta*C
    dD = delta*C

    dAbsTMZ = -ka1 * AbsTMZ
    dPlaTMZ = ka1 * AbsTMZ - Cl1 / VD1 * PlaTMZ - k23 * PlaTMZ + k32 * CSFTMZ
    dCSFTMZ = k23 * PlaTMZ - k32 * CSFTMZ

    El = kel * PlaGDC
    dAbsGDC = -ka2 * AbsGDC
    dPlaGDC = ka2 * AbsGDC - El + k21 * PeriphGDC - k12 * PlaGDC
    dPeriphGDC = -k21 * PeriphGDC + k12 * PlaGDC

    dcAUC = C

    du .= [dC, dD, dAbsTMZ, dPlaTMZ, dCSFTMZ, dAbsGDC, dPlaGDC, dPeriphGDC, dcAUC]
end

function drug_pk!(du, u, p, t)
    kel, ka2, V2, k12, k21 = p[1:end-1]
    AbsGDC, PlaGDC, PeriphGDC = u
    dAbsGDC = -ka2 * AbsGDC
    dPlaGDC = ka2 * AbsGDC - kel * PlaGDC + k21 * PeriphGDC - k12 * PlaGDC
    dPeriphGDC = -k21 * PeriphGDC + k12 * PlaGDC
    du .= [dAbsGDC, dPlaGDC, dPeriphGDC]
end

