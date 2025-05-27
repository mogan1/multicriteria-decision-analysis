* Cost-minimizing LP for generation expansion: Multi-zone MCDA implementation

*$INCLUDE mcda_test_data_2024.gms
$INCLUDE mcda_nordic_data_2024.gms

Scalar starttime;
starttime = jnow;

Scalars
DN Set to '1' to turn off all adoption /0/
NAC Set to '1' to turn off all availability costs /0/
NIC Set to '1' to turn off all investment costs /0/
NX Set to '1' to turn off all net imports /0/
NUC Indicate derating of nuclear capacity /0/
;

Variables
OF_Cost        Objective function value for cost minimisation
OF_Emissions   Objective function value for emission minimisation
q              Objective function value for minimax implementation
f(L,T)         Power flow on line l in period t
;

Positive Variables
a_gen(N,U)    Generation capacity made available of thermal unit u at node n
a_vre(N,E)    Generation capacity made available of VRE unit e at node n
a_hyd(N,W)    Generation capacity made available of hydro unit w at node n
b_vre(N,E)    Generation-capacity adopted of VRE unit e at node n
dd(N,T)       DSR activated by node n in period t
gg_gen(N,T,U) Generation by thermal unit u at node n in period t
gg_vre(N,T,E) Generation by VRE unit e at node n in period t
r_in(N,T,W)   Volume of water pumped by hydro unit w at node n in period t
r_out(N,T,W)  Volume of water turbined by hydro unit w at node n in period t
r_sto(N,T,W)  Volume of water stored by hydro unit w at node n in period t
z(N,T,W)      Volume of water spilled by hydro unit w at node n in period t
;

* Added to prevent capacity availability where it does not exist
a_gen.fx(N,U)$(NOT ntou(N,U))=0;
a_vre.fx(N,E)$(NOT ntoe(N,E))=0;
a_hyd.fx(N,W)$(NOT ntow(N,W))=0;

* Added to turn off investment option
b_vre.fx(N,E)$(DN eq 1)=0;

* Added to turn off availability costs
C_ava_gen(U)$(NAC eq 1)=0;
C_ava_vre(E)$(NAC eq 1)=0;
C_ava_hyd(W)$(NAC eq 1)=0;

* Added to turn off investment costs
C_inv_vre(E)$(NIC eq 1)=0;

* Added to turn off net imports
X(T,N)$(NX eq 1)=0;

* Added to derate nuclear capacity
G_gen(N,'nuclear')=(1-NUC)*G_gen(N,'nuclear');

Display G_gen;

Equations
CostEq              Objective function equation for cost minimisation
EmissionEq          Objective function equation for emission minimisation
EnergyBalance       Energy balance by node and by period
FlowPosLimit        Power-flow limit in positive direction by line and by period
FlowNegLimit        Power-flow limit in negative direction by line and by period
LoadShedLimit       Maximum load shedding by node and by period
GenCapacityLimit    Thermal generation limit for each technology by period
GenAvailLimit       Thermal capacity limit for each technology
VRECapacityLimit    VRE generation limit for each technology by period
VREAvailLimit       VRE capacity limit for each technology
VREInvLimit         VRE adoption limit for each technology
GenUpRampLimit      Thermal generation up ramping by technology and by period
GenDownRampLimit    Thermal generation down ramping by technology and by period
StorageBalance      Storage balance by hydro unit at each node and by period
StorageMaxLimit     Maximum storage volume by technology and by period
StorageMinLimit     Minimum storage volume by technology and by period
HydroPumpLimit      Maximum pumped-hydro charging by technology and by period
HydroCapacityLimit  Hydro generation limit for each technology by period
HydroAvailLimit     Hydro capacity limit for each technology
CostGoal            Cost goal in minimax implementation
EmissionGoal        Emission goal in minimax implementation
CostAux             Auxiliary cost constraint to obtain maximum emissions in cost-minimization mode
CostAux2            Auxiliary cost constraint in PGP emission-minimization mode
EmissionsAux        Auxiliary emission constraint to obtain maximum cost in emission-minimization mode
;

CostEq.. OF_Cost=E=SUM((N,T,U), (C_opr(U)+S*P(U))*gg_gen(N,T,U))+SUM((N,T), C_dsr(T,N)*dd(N,T))+SUM((N,E), C_ava_vre(E)*a_vre(N,E) + C_inv_vre(E)*b_vre(N,E))+SUM((N,U), C_ava_gen(U)*a_gen(N,U))+SUM((N,W), C_ava_hyd(W)*a_hyd(N,W));
EmissionEq.. OF_Emissions=E=SUM((N,T,U), P(U)*gg_gen(N,T,U));
EnergyBalance(N,T).. dd(N,T)+SUM(U, gg_gen(N,T,U))+SUM(E, gg_vre(N,T,E))+SUM(W, Q_hyd(N,W)*r_out(N,T,W)-F_hyd(N,W)*r_in(N,T,W))+X(T,N)+TT(T)*V*SUM(L$NMinus(L,N), f(L,T))-TT(T)*V*SUM(L$NPlus(L,N), f(L,T))-D(T,N)=E=0;
FlowPosLimit(L,T).. TT(T)*K_pos(L)-V*TT(T)*f(L,T)=G=0;
FlowNegLimit(L,T).. TT(T)*K_neg(L)+V*TT(T)*f(L,T)=G=0;
LoadShedLimit(N,T).. TT(T)*D_dsr(T,N)-dd(N,T)=G=0;
GenCapacityLimit(N,T,U).. TT(T)*a_gen(N,U)-gg_gen(N,T,U)=G=0;
GenAvailLimit(N,U).. G_gen(N,U)-a_gen(N,U)=G=0;
VRECapacityLimit(N,T,E).. A(T,E,N)*TT(T)*a_vre(N,E)-gg_vre(N,T,E)=G=0;
VREAvailLimit(N,E).. G_vre(N,E)+b_vre(N,E)-a_vre(N,E)=G=0;
VREInvLimit(N,E).. M(N,E)*G_vre(N,E)-b_vre(N,E)=G=0;
GenUpRampLimit(N,T,U)$( ORD(T) ge 2 ).. -gg_gen(N,T,U)+gg_gen(N,T-1,U)+TT(T)*R_up(U)*a_gen(N,U)=G=0;
GenDownRampLimit(N,T,U)$( ORD(T) ge 2 ).. gg_gen(N,T,U)-gg_gen(N,T-1,U)+TT(T)*R_down(U)*a_gen(N,U)=G=0;
StorageBalance(N,T,W).. -r_sto(N,T,W)+r_sto(N,T-1,W)$(ORD(T) > 1)+RR_ini(N,W)$( ORD(T) eq 1 )+r_in(N,T,W)-r_out(N,T,W)-z(N,T,W)+II_hourly_ror(N,T)$( ORD(W) eq 1)+II_hourly_res(N,T)$( ORD(W) eq 2)=E=0;
StorageMaxLimit(N,T,W).. RR_max(N,W)-r_sto(N,T,W)=G=0;
StorageMinLimit(N,T,W)$(ORD(T) eq CARD(T) ).. -RR_min(N,W)+r_sto(N,T,W)=G=0;
HydroPumpLimit(N,T,W).. TT(T)*RR_in(N,W)*RR_max(N,W)-r_in(N,T,W)=G=0;
HydroCapacityLimit(N,T,W).. TT(T)*a_hyd(N,W)-Q_hyd(N,W)*r_out(N,T,W)=G=0;
HydroAvailLimit(N,W).. Y_hyd(N,W)-a_hyd(N,W)=G=0;
CostGoal.. -W_C*(OF_Cost-Cmin)/Cmin+q=G=0;
EmissionGoal.. -W_E*(OF_Emissions-Emin)/Emin+q=G=0;
CostAux.. -OF_Cost+Cmin=G=0;
CostAux2.. -OF_Cost+Chat=G=0;
EmissionsAux.. -OF_Emissions+Emin=G=0;

MODEL GEP_LP /CostEq, EmissionEq, EnergyBalance, FlowPosLimit, FlowNegLimit, LoadShedLimit, GenCapacityLimit, GenAvailLimit, VRECapacityLimit, VREAvailLimit, VREInvLimit, GenUpRampLimit, GenDownRampLimit, StorageBalance, StorageMaxLimit, StorageMinLimit, HydroPumpLimit, HydroCapacityLimit, HydroAvailLimit/;
option LP = CPLEX;
*option LP = GUROBI;
option optcr = 0.0000;

SOLVE GEP_LP USING LP MINIMIZING OF_Cost;

$ontext
Cmin = OF_Cost.L;
Emax = OF_Emissions.L;

MODEL GEP_LP_Cost_Aux /CostEq, EmissionEq, EnergyBalance, FlowPosLimit, FlowNegLimit, LoadShedLimit, GenCapacityLimit, GenAvailLimit, VRECapacityLimit, VREAvailLimit, VREInvLimit, GenUpRampLimit, GenDownRampLimit, StorageBalance, StorageMaxLimit, StorageMinLimit, HydroPumpLimit, HydroCapacityLimit, HydroAvailLimit, CostAux/;

SOLVE GEP_LP_Cost_Aux USING LP MINIMIZING OF_Emissions;
$offtext

Parameters Tot_Cons_C, Tot_Gen_C, Tot_Th_Gen_C, Tot_VRE_Gen_C, Tot_Hyd_Gen_C, Flow_C(L,T), Price_NT_C(N,T), Ave_Price_N_C(N);

Tot_Cons_C=SUM((T,N), D(T,N));
Tot_Th_Gen_C=SUM((N,T,U), gg_gen.L(N,T,U));
Tot_VRE_Gen_C=SUM((N,T,E), gg_vre.L(N,T,E));
Tot_Hyd_Gen_C=SUM((N,T,W), Q_hyd(N,W)*r_out.L(N,T,W));
Tot_Gen_C=Tot_Th_Gen_C+Tot_VRE_Gen_C+Tot_Hyd_Gen_C;
Flow_C(L,T)=f.L(L,T);
Price_NT_C(N,T)=EnergyBalance.M(N,T);
Ave_Price_N_C(N)=SUM(T, Price_NT_C(N,T))/CARD(T);

Display Tot_Cons_C, Tot_Gen_C, Tot_Th_Gen_C, Tot_VRE_Gen_C, Tot_Hyd_Gen_C, Flow_C, Price_NT_C, Ave_Price_N_C;

Cmin = OF_Cost.L;
Emax = OF_Emissions.L;

SOLVE GEP_LP USING LP MINIMIZING OF_Emissions;

Emin = OF_Emissions.L;
Cmax = OF_Cost.L;

MODEL GEP_LP_Emissions_Aux /CostEq, EmissionEq, EnergyBalance, FlowPosLimit, FlowNegLimit, LoadShedLimit, GenCapacityLimit, GenAvailLimit, VRECapacityLimit, VREAvailLimit, VREInvLimit, GenUpRampLimit, GenDownRampLimit, StorageBalance, StorageMaxLimit, StorageMinLimit, HydroPumpLimit, HydroCapacityLimit, HydroAvailLimit, EmissionsAux/;

SOLVE GEP_LP_Emissions_Aux USING LP MINIMIZING OF_Cost;

Parameters Tot_Cons_E, Tot_Gen_E, Tot_Th_Gen_E, Tot_VRE_Gen_E, Tot_Hyd_Gen_E, Flow_E(L,T), Price_NT_E(N,T), Ave_Price_N_E(N);

Tot_Cons_E=SUM((T,N), D(T,N));
Tot_Th_Gen_E=SUM((N,T,U), gg_gen.L(N,T,U));
Tot_VRE_Gen_E=SUM((N,T,E), gg_vre.L(N,T,E));
Tot_Hyd_Gen_E=SUM((N,T,W), Q_hyd(N,W)*r_out.L(N,T,W));
Tot_Gen_E=Tot_Th_Gen_E+Tot_VRE_Gen_E+Tot_Hyd_Gen_E;
Flow_E(L,T)=f.L(L,T);
Price_NT_E(N,T)=EnergyBalance.M(N,T);
Ave_Price_N_E(N)=SUM(T, Price_NT_E(N,T))/CARD(T);

Display Tot_Cons_E, Tot_Gen_E, Tot_Th_Gen_E, Tot_VRE_Gen_E, Tot_Hyd_Gen_E, Flow_E, Price_NT_E, Ave_Price_N_E;

Emin = OF_Emissions.L;
Cmax = OF_Cost.L;

Display Cmin, Emax, Cmax, Emin;


*$ontext

*MODEL GEP_Minimax /CostEq, EmissionEq, EnergyBalance, FlowPosLimit, FlowNegLimit, LoadShedLimit, GenCapacityLimit, GenAvailLimit, VRECapacityLimit, VREAvailLimit, VREInvLimit, GenUpRampLimit, GenDownRampLimit, StorageBalance, StorageMaxLimit, StorageMinLimit, HydroPumpLimit, HydroCapacityLimit, HydroAvailLimit, CostGoal, EmissionGoal/;

MODEL GEP_LP_Cost_Aux /CostEq, EmissionEq, EnergyBalance, FlowPosLimit, FlowNegLimit, LoadShedLimit, GenCapacityLimit, GenAvailLimit, VRECapacityLimit, VREAvailLimit, VREInvLimit, GenUpRampLimit, GenDownRampLimit, StorageBalance, StorageMaxLimit, StorageMinLimit, HydroPumpLimit, HydroCapacityLimit, HydroAvailLimit, CostAux2/;

*SOLVE GEP_Minimax USING LP MINIMIZING q;

*Display q.L;

Parameters TotalCost_PGP(it), TotalEmissions_PGP(it), a_vre_PGP(it,N,E), b_vre_PGP(it,N,E), a_gen_PGP(it,N,U), a_hyd_PGP(it,N,W);

Loop (it,
Chat=Cmin+(Cmax-Cmin)*(ord(it)-1)*(1/(CARD(it)-1));
Display Chat;
*W_C = (ord(it)-1)*(1/(CARD(it)-1));
*W_E = 1-W_C;
option optcr = 0.0000;
SOLVE GEP_LP_Cost_Aux USING LP MINIMIZING OF_Emissions;

TotalCost_PGP(it)=OF_Cost.L;
TotalEmissions_PGP(it)=OF_Emissions.L;
a_vre_PGP(it,N,E)=a_vre.L(N,E);
b_vre_PGP(it,N,E)=b_vre.L(N,E);
a_gen_PGP(it,N,U)=a_gen.L(N,U);
a_hyd_PGP(it,N,W)=a_hyd.L(N,W);
);


execute_unload 'res_MCDA_2024.gdx', TotalCost_PGP, TotalEmissions_PGP, a_vre_PGP, b_vre_PGP, a_gen_PGP, a_hyd_PGP;
execute 'gdxxrw.exe res_MCDA_2024.gdx par=TotalCost_PGP rng=TotalCost!a1';
execute 'gdxxrw.exe res_MCDA_2024.gdx par=TotalEmissions_PGP rng=TotalEmissions!a1';
execute 'gdxxrw.exe res_MCDA_2024.gdx par=a_vre_PGP rng=a_vre!a1';
execute 'gdxxrw.exe res_MCDA_2024.gdx par=b_vre_PGP rng=b_vre!a1';
execute 'gdxxrw.exe res_MCDA_2024.gdx par=a_gen_PGP rng=a_gen!a1';
execute 'gdxxrw.exe res_MCDA_2024.gdx par=a_hyd_PGP rng=a_hyd!a1';
*$offtext

scalar elapsed; elapsed = (jnow - starttime)*24*3600;

display elapsed