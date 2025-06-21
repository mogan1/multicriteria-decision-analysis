* Cost-minimizing LP for energy portfolio in case of SE multiple bidding zones - An MCDA implementation
* Author: Pushan Deb and Mousam Ganguly 

$TITLE Combined LP Model with Cost, Emission, and CVaR Minimization

$GDXIN in.gdx
$LOAD *
$GDXIN

Scalars
* Emin/Emax and Cmin/Cmax are defined in INCLUDE
    Qmin                    'Scalar variable: Minimized CVaR'
    Qmax                    'Scalar variable: Maximised CVaR'
    C2                      'Scalar variable: Intermediary cost during emission minimization'
    C3                      'Scalar variable: Intermediary cost during CVaR minimization'
    E1                      'Scalar variable: Intermediary emission during cost minimization'
    E3                      'Scalar variable: Intermediary emission during CVaR minimization'
    Q1                      'Scalar variable: Intermediary CVaR during cost minimization'
    Q2                      'Scalar variable: Intermediary CVaR during emission minimization'
    alpha                   'Scalar variable: Confidence level'
;


Variables
    of_cost                 'Variable: Objective function value for cost minimization'
    of_emission             'Variable: Objective function value for emission minimization'
    of_cvar                 'Variable: Objective function value for CVaR of residual load minimization'
    f(L,T)                  'Variable: Power flow on transmission line l in period t [MW]'
    h(N)                    'Variable: Minimum threshold residual load at VaR at node n [MWh]'
    ;
    
Positive Variables
    g(N,T,U)                'Variable: Generation of thermal unit u at node n in period t [MWh]'
    dsr(N,T)                'Variable: Demand-side response activated at node n at time T [MWh]'
    b_vre(N,E)              'Variable: Adopted capacity of VRE unit e at node n [MW]'
    a_gen(N,U)              'Variable: Available capacity of thermal unit u at node n [MW]'
    a_vre(N,E)              'Variable: Available capacity of VRE unit e at node n [MW]'
    a_hyd(N,W)              'Variable: Available capacity of hydro unit w at node n [MW]'
    z(N,T)                  'Variable: Amount of residual load that falls above VaR threshold at node n in period t [MWh]'
    uv(N,T,E)                'Variable: VRE generation level by unit e at node n in period T [MWh]'
    r_in(N,T,W)             'Variable: Volume of water pumped into hydro unit w at node n in period t [m3]'
    r_out(N,T,W)            'Variable: Volume of water turbined out from hydro unit w at node n in period t [m3]'
    r_sto(N,T,W)            'Variable: Volume of water stored by hydro unit w at node n in period t [m3]'
    sp(N,T,W)                'Variable: Volume of water spilled from hydro unit w at node n in period t [m3]' 
;

Equations
    cost_eq                 'Objective function equation: Minimize total cost'
    emission_eq             'Objective function equation: Minimize total CO₂ emissions'
    cvar_eq                 'Objective function equation: Minimize CVaR of residual load'
    
    energy_balance          'Constraint: Energy balance (supply must meet demand)'
    flow_pos_limit          'Constraint: Power-flow limit in positive direction (inward) by line and by period'
    flow_neg_limit          'Constraint: Power-flow limit in negative direction (outward) by line and by period'
    load_shed_limit         'Constraint: Maximum allowable load shedding by node and time period'
    gen_capacity_limit      'Constraint: Thermal generation capped by available capacity for each technology and time'
    gen_avail_limit         'Constraint: Thermal generation must not exceed installed capacity per technology'
    vre_capacity_limit      'Constraint: VRE generation limitation on available capacity'
    vre_avail_limit         'Constraint: VRE generation availability limitation on installed capacity'
    vre_inv_limit           'Constraint: Limitation on additional VRE generation capacity'
    gen_up_ramp_limit       'Constraint: Thermal generation ramp-up limit between consecutive time periods'
    gen_down_ramp_limit     'Constraint: Thermal generation ramp-down limit between consecutive time periods'
    storage_balance         'Constraint: Storage state-of-charge balance for hydro units across time'
    storage_max_limit       'Constraint: Upper bound on storage capacity for hydro technologies by period'
    storage_min_limit       'Constraint: Lower bound on storage level for hydro technologies by period'
    hydro_pump_limit        'Constraint: Pumped-hydro charging limit per technology and period'
    hydro_capacity_limit    'Constraint: Hydro generation limit based on available water and capacity'
    hydro_avail_limit       'Constraint: Hydro generation restricted by installed capacity'
    residual_load_balance   'Constraint: Residual load deviation'
    emission_cap            'Limit for emission'
    cvar_cap                'Limit for emission'


cost_eq..                                       of_cost =e= sum((N,T,U), (C_opr(U)+S*P(U))*g(N,T,U)) + sum((N,T), C_dsr(T,N)*dsr(N,T)) + sum((N,E), C_inv_vre(E)*b_vre(N,E)) + sum((N,U), C_ava_gen(U)*a_gen(N,U)) + sum((N,E), C_ava_vre(E)*a_vre(N,E)) + sum((N,W), C_ava_hyd(W)*a_hyd(N,W));
emission_eq..                                   of_emission =e= sum((N,T,U), P(U)*g(N,T,U));
cvar_eq..                                       of_cvar =e= sum(N, (h(N) + (1/(1-alpha) * card(T)) * sum(T, z(N,T))));
energy_balance(N,T)..                           sum(U, g(N,T,U)) + sum(E, uv(N,T,E)) + sum(W, Q_hyd(N,W)*r_out(N,T,W) - F_hyd(N,W)*r_in(N,T,W)) + sum(L $ NMinus(L,N), V*TT(T)*f(L,T)) - sum(L $ NPlus(L,N), V*TT(T)*f(L,T)) + dsr(N,T) + X(T,N) - D(T,N) =e= 0;
flow_pos_limit(L,T)..                           TT(T)*K_pos(L)-V*TT(T)*f(L,T) =g= 0;
flow_neg_limit(L,T)..                           V*TT(T)*f(L,T) + TT(T)*K_neg(L) =g= 0;
load_shed_limit(N,T)..                          TT(T)*D_dsr(T,N) - dsr(N,T) =g= 0;
gen_capacity_limit(N,T,U)..                     TT(T)*a_gen(N,U) - g(N,T,U) =g= 0;
gen_avail_limit(N,U)..                          G_gen(N,U) - a_gen(N,U) =g= 0;
gen_up_ramp_limit(N,T,U)$( ORD(T) ge 2 )..      TT(T)*R_up(U)*a_gen(N,U) - g(N,T,U) + g(N,T-1,U) =g= 0;
gen_down_ramp_limit(N,T,U)$( ORD(T) ge 2 )..    g(N,T,U) - g(N,T-1,U) + TT(T)*R_down(U)*a_gen(N,U) =g= 0;
vre_capacity_limit(N,T,E)..                     TT(T)*A(T,E,N)*a_vre(N,E) - uv(N,T,E) =g= 0;
vre_avail_limit(N,E)..                          G_vre(N,E) + b_vre(N,E) - a_vre(N,E) =g= 0;
vre_inv_limit(N,E)..                            M(N,E)*G_vre(N,E) - b_vre(N,E) =g= 0;
storage_balance(N,T,W)..                        - r_sto(N,T,W) + r_sto(N,T-1,W)$(ORD(T) > 1) + RR_ini(N,W)$( ORD(T) eq 1 ) + r_in(N,T,W) - r_out(N,T,W) - sp(N,T,W) + II_hourly_ror(N,T)$( ORD(W) eq 1) + II_hourly_res(N,T)$( ORD(W) eq 2)  =e= 0;
storage_max_limit(N,T,W)..                      RR_max(N,W) - r_sto(N,T,W) =g= 0;
storage_min_limit(N,T,W)$(ORD(T) eq CARD(T))..  r_sto(N,T,W) - RR_min(N,W) =g= 0;
hydro_pump_limit(N,T,W)..                       TT(T)*RR_in(N,W)*RR_max(N,W) - r_in(N,T,W) =g= 0;
hydro_capacity_limit(N,T,W)..                   TT(T)*a_hyd(N,W) - Q_hyd(N,W)*r_out(N,T,W) =g= 0;
hydro_avail_limit(N,T,W)..                      Y_hyd(N,W) - a_hyd(N,W) =g= 0;
residual_load_balance(N,T)..                    z(N,T) + h(N) - sum(U, g(N,T,U)) =g= 0;

option LP = CPLEX;
option optcr = 0.0000;
alpha = 0.95;


*Model 1: Cost Minimization ===
model cost_lp / cost_eq, energy_balance, flow_pos_limit, flow_neg_limit, load_shed_limit, gen_capacity_limit, gen_avail_limit, gen_up_ramp_limit, gen_down_ramp_limit, vre_capacity_limit, vre_avail_limit, vre_inv_limit, storage_balance, storage_max_limit, storage_min_limit, hydro_pump_limit, hydro_capacity_limit, hydro_avail_limit, residual_load_balance /;
cost_lp.optfile = 1;

Solve cost_lp using lp minimizing of_cost;

Cmin = of_cost.L;
E1 = sum((N,T,U), P(U)*g.L(N,T,U));
Q1 = sum(N, (h.L(N) + (1/(1-alpha) * card(T)) * sum(T, z.L(N,T))));


Scalar run_mode;
run_mode = 1;

* Report Summary
Parameters
    Tot_Cons              'Total consumption across all nodes and times [MWh]'
    Tot_Th_Gen            'Total thermal generation [MWh]'
    Tot_VRE_Gen           'Total VRE generation [MWh]'
    Tot_Hyd_Gen           'Total hydro generation [MWh]'
    Tot_Gen               'Total generation [MWh]'
    Flow(L,T)             'Power flow [MW]'
    Price_NT(N,T)         'Nodal marginal price [€/MWh]'
    Ave_Price_N(N)        'Average nodal price [€/MWh]'
    Gen_Node_Tech(N,U)    'Node-level thermal generation by tech [MWh]'
    VRE_Node_Tech(N,E)    'Node-level VRE generation by tech [MWh]'
    Hyd_Node_Tech(N,W)    'Node-level hydro generation by tech [MWh]'
    Flow_Tot(L)           'Total flow across all time periods per line [MWh]';

Price_NT_Final(N,T) = EnergyBalance.M(N,T);
Ave_Price_Final(N)  = SUM(T, Price_NT_Final(N,T)) / CARD(T);
Tot_Th_Gen_Final    = SUM((N,T,U), gg_gen.L(N,T,U));
Tot_VRE_Gen_Final   = SUM((N,T,E), gg_vre.L(N,T,E));
Tot_Hyd_Gen_Final   = SUM((N,T,W), Q_hyd(N,W)*r_out.L(N,T,W));
Tot_Cons_Final      = SUM((T,N), D(T,N));
Tot_Gen_Final       = Tot_Th_Gen_Final + Tot_VRE_Gen_Final + Tot_Hyd_Gen_Final;
Flow_Avg_Final(L)   = SUM(T, f.L(L,T)) / CARD(T);

* Export all key results to GDX
execute_unload 'output_cost_summary.gdx',
    Tot_Cons, Tot_Th_Gen, Tot_VRE_Gen, Tot_Hyd_Gen, Tot_Gen,
    Flow, Price_NT, Ave_Price_N,
    Gen_Node_Tech, VRE_Node_Tech, Hyd_Node_Tech,Flow_Tot;


file fsum /'results.txt'/;
put fsum;
put 'Run Type,Tot_Cons,Tot_Th_Gen,Tot_VRE_Gen,Tot_Hyd_Gen,Tot_Gen',/;
put 'Cost', ',', Tot_Cons_Final:12:2, ',', Tot_Th_Gen_Final:12:2, ',', Tot_VRE_Gen_Final:12:2, ',',
    Tot_Hyd_Gen_Final:12:2, ',', Tot_Gen_Final:12:2, /;

*Model 2: Emission Minimization ===
model emission_lp / emission_eq, energy_balance, flow_pos_limit, flow_neg_limit, load_shed_limit, gen_capacity_limit, gen_avail_limit, gen_up_ramp_limit, gen_down_ramp_limit, vre_capacity_limit, vre_avail_limit, vre_inv_limit, storage_balance, storage_max_limit, storage_min_limit, hydro_pump_limit, hydro_capacity_limit, hydro_avail_limit, residual_load_balance /;

emission_lp.optfile = 1;
Solve emission_lp using lp minimizing of_emission;
execute_unload 'output_emission_summary.gdx',
    Tot_Cons, Tot_Th_Gen, Tot_VRE_Gen, Tot_Hyd_Gen, Tot_Gen,
    Flow, Price_NT, Ave_Price_N,
    Gen_Node_Tech, VRE_Node_Tech, Hyd_Node_Tech,Flow_Tot;

run_mode = 2;

Price_NT_Final(N,T) = EnergyBalance.M(N,T);
Ave_Price_Final(N)  = SUM(T, Price_NT_Final(N,T)) / CARD(T);
Tot_Th_Gen_Final    = SUM((N,T,U), gg_gen.L(N,T,U));
Tot_VRE_Gen_Final   = SUM((N,T,E), gg_vre.L(N,T,E));
Tot_Hyd_Gen_Final   = SUM((N,T,W), Q_hyd(N,W)*r_out.L(N,T,W));
Tot_Cons_Final      = SUM((T,N), D(T,N));
Tot_Gen_Final       = Tot_Th_Gen_Final + Tot_VRE_Gen_Final + Tot_Hyd_Gen_Final;
Flow_Avg_Final(L)   = SUM(T, f.L(L,T)) / CARD(T);

put 'Emission', ',', Tot_Cons_Final:12:2, ',', Tot_Th_Gen_Final:12:2, ',', Tot_VRE_Gen_Final:12:2, ',',
    Tot_Hyd_Gen_Final:12:2, ',', Tot_Gen_Final:12:2, /;

* === Model 3: CVaR Minimization ===
Model cvar_lp / all /;
cvar_lp.optfile = 1;
Solve cvar_lp using lp minimizing of_cvar;

run_mode = 3;

Price_NT_Final(N,T) = EnergyBalance.M(N,T);
Ave_Price_Final(N)  = SUM(T, Price_NT_Final(N,T)) / CARD(T);
Tot_Th_Gen_Final    = SUM((N,T,U), gg_gen.L(N,T,U));
Tot_VRE_Gen_Final   = SUM((N,T,E), gg_vre.L(N,T,E));
Tot_Hyd_Gen_Final   = SUM((N,T,W), Q_hyd(N,W)*r_out.L(N,T,W));
Tot_Cons_Final      = SUM((T,N), D(T,N));
Tot_Gen_Final       = Tot_Th_Gen_Final + Tot_VRE_Gen_Final + Tot_Hyd_Gen_Final;
Flow_Avg_Final(L)   = SUM(T, f.L(L,T)) / CARD(T);

execute_unload 'output_cvar_summary.gdx',
    Tot_Cons, Tot_Th_Gen, Tot_VRE_Gen, Tot_Hyd_Gen, Tot_Gen,
    Flow, Price_NT, Ave_Price_N,
    Gen_Node_Tech, VRE_Node_Tech, Hyd_Node_Tech,Flow_Tot;


put 'CVaR', ',', Tot_Cons_Final:12:2, ',', Tot_Th_Gen_Final:12:2, ',', Tot_VRE_Gen_Final:12:2, ',',
    Tot_Hyd_Gen_Final:12:2, ',', Tot_Gen_Final:12:2, /;

putclose fsum;
