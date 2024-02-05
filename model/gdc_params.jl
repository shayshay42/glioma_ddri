include("../assets/pk_params.jl")

gamma_1 =0.5689#0.6689
psi     =1
C0      =17.7
D0      =0.0
r       =0.007545/24
K       =158.04
BW      =70
IC50_1  =5.807*0.000001*194.151
Imax_1  =0.905#0.605
IC50_2  =6.702*0.000001*458#0.0005364#6.702#2.48#20#0.0005364
gamma_2 =0.916#0.8992#0.8502#0.1969#0.26#0.76
Imax_2  =0.8824#0.9357#0.9834#1.306#0.86#0.96#0.26
xi      =IC50_1/IC50_2
VD1     =30.3#(FIXED)
Cl1     =TMZ_params["Cl1"].mu * TMZ_params["Cl1"].mult#9.64999551636561#10.5705229706946(VARIED)
k23     =TMZ_params["k23"].mu * TMZ_params["k23"].mult#0.000791962919186130#0.000823753578728557(VARIED)
ka1     =TMZ_params["ka1"].mu * TMZ_params["ka1"].mult#2.68455859105533#9.75543575220511(VARIED)
k32     =0.76#(FIXED)
ka2     =GDC_params["ka2"].mu * GDC_params["ka2"].mult#0.889860150615071#0.939691715960798(VARIED)
V2      =GDC_params["V2"].mu  * GDC_params["V2"].mult#120.826298714951#123.004969002618(VARIED)
kel     =GDC_params["kel"].mu * GDC_params["kel"].mult#2.27911938377711#1.25132101717989(VARIED)
k12     =GDC_params["k12"].mu * GDC_params["k12"].mult#74.4##V
k21     =GDC_params["k21"].mu * GDC_params["k21"].mult#26.4##V

ode_params = [gamma_1,psi,C0,D0,r,K,BW,IC50_1,Imax_1,IC50_2,gamma_2,Imax_2,xi,VD1,Cl1,k23,ka1,k32,ka2,V2,kel,k12,k21];

param_order = ["gamma_1", "psi", "C0", "D0", "r", "K", "BW", "IC50_1", "Imax_1", "IC50_2", "gamma_2", "Imax_2", "xi", "VD1", "Cl1", "k23", "ka1", "k32", "ka2", "V2", "kel", "k12", "k21"]
const p_num = length(ode_params)
default_scaling = [1.0,1.0]
const nb_scaling_params = length(default_scaling)