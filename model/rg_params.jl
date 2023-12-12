include("../assets/pk_params.jl")

gamma_1 = 0.5689#0.42676#0.5689
psi     = 1#0.8
C0      = 17.7
D0      = 0.0
r       = 0.007545/24
K       = 158.04
BW      = 70.0
IC50_1  = 5.807*0.000001*194.151#15.5936*0.000001*194.151#5.807*0.000001*194.151
Imax_1  = 0.905#1.1026#0.905
IC50_2  = 0.009238
gamma_2 = 2.712
Imax_2  = 1.0
xi      = IC50_1/IC50_2
VD1     = 30.3 
Cl1     = TMZ_params["Cl1"].mu * TMZ_params["Cl1"].mult#10.5705229706946 #V
k23     = TMZ_params["k23"].mu * TMZ_params["k23"].mult#0.000823753578728557 #V
ka1     = TMZ_params["ka1"].mu * TMZ_params["ka1"].mult#9.75543575220511 #V
k32     = 0.76
Cl2     = RG_params["Cl2"].mu  * RG_params["Cl2"].mult#32.6682 #V
ka2     = RG_params["ka2"].mu  * RG_params["ka2"].mult#0.0385233 #V
Vpla    = RG_params["Vpla"].mu * RG_params["Vpla"].mult#0.934662 #V
Q       = RG_params["Q"].mu    * RG_params["Q"].mult#0.0302696 #V
Vtis    = RG_params["Vtis"].mu * RG_params["Vtis"].mult#0.00299745 #V

ode_params = [gamma_1,psi,C0,D0,r,K,BW,IC50_1,Imax_1,IC50_2,gamma_2,Imax_2,xi,VD1,Cl1,k23,ka1,k32,Cl2,ka2,Vpla,Q,Vtis];
param_order = ["gamma_1", "psi", "C0", "D0", "r", "K", "BW", "IC50_1", "Imax_1", "IC50_2", "gamma_2", "Imax_2", "xi", "VD1", "Cl1", "k23", "ka1", "k32", "Cl2", "ka2", "Vpla", "Q", "Vtis"]
const p_num = length(ode_params)
# Create the named tuple
# ode_p = NamedTuple{Tuple(Symbol.(param_order))}(ode_params)


# ics = (C=17.0, D=0.0, AbsTMZ=0.0, PlaTMZ=0.0, CSFTMZ=0.0, AbsRG=0.0, PlaRG=0.0, TisRG=0.0, cAUC=0.0)
# u0 = [ics.C, ics.D, ics.AbsTMZ, ics.PlaTMZ, ics.CSFTMZ, ics.AbsRG, ics.PlaRG, ics.TisRG, ics.cAUC]