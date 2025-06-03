* Cost-minimizing LP for energy portfolio in case of SE multiple bidding zones - An MCDA implementation
* Author: Pushan Deb and Mousam Ganguly 

$INCLUDE mcda_nordic_data_2024.gms

Variables
    of_cost                 'Variable: Objective function value for cost minimization'
    of_emission             'Variable: Objective function value for emission minimization'
    of_cvar                 'Variable: Objective function value for CVaR of residual load minimization'
    f(L,T)                  'Variable: Power flow on transmission line l in period t [MW]';
    
Positive Variables
    g(N,T,U)                'Variable: Generation of thermal unit u at node n in period t [MWh]'
    dsr(N,T)                'Variable: Demand-side response activated at node n at time T [MWh]'
    b_vre(N,E)              'Variable: Adopted capacity of VRE unit e at node n [MW]'
    a_gen(N,U)              'Variable: Available capacity of thermal unit u at node n [MW]'
    a_vre(N,E)              'Variable: Available capacity of VRE unit e at node n [MW]'
    a_hyd(N,W)              'Variable: Available capacity of hydro unit w at node n [MW]'
    h(N)                    'Variable: Minimum threshold residual load at VaR at node n'
    z(N,T)                  'Variable: Amount of residual load that falls above VaR threshold at node n in period t [MWh]'
    u(N,T,E)                'Variable: VRE generation level by unit e at node n in period T [MWh]'
    r_in(N,T,W)             'Variable: Volume of water pumped into hydro unit w at node n in period t [m3]'
    r_out(N,T,W)            'Variable: Volume of water turbined out from hydro unit w at node n in period t [m3]'
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
    vre_capacity_limit      'Constraint: VRE generation capped by available VRE output in each period'
    vre_avail_limit         'Constraint: VRE generation limited by installed capacity per technology'
    vre_inv_limit           'Constraint: Limit on new VRE adoption per technology'
    gen_up_ramp_limit       'Constraint: Thermal generation ramp-up limit between consecutive time periods'
    gen_down_ramp_limit     'Constraint: Thermal generation ramp-down limit between consecutive time periods'
    storage_balance         'Constraint: Storage state-of-charge balance for hydro units across time'
    storage_max_limit       'Constraint: Upper bound on storage capacity for hydro technologies by period'
    storage_min_limit       'Constraint: Lower bound on storage level for hydro technologies by period'
    hydro_pump_limit        'Constraint: Pumped-hydro charging limit per technology and period'
    hydro_capacity_limit    'Constraint: Hydro generation limit based on available water and capacity'
    hydro_avail_limit       'Constraint: Hydro generation restricted by installed capacity'
    
*    CostGoal                'Auxiliary constraint: Cost goal (used in minimax or benchmarking)'
*    EmissionGoal            'Auxiliary constraint: Emission goal (used in minimax or benchmarking)'
*    CostAux                 'Auxiliary constraint: Upper bound on cost when minimizing emissions'
*    CostAux2                'Auxiliary constraint: Cost control in policy-guided planning (PGP) emission minimization'
*    EmissionsAux            'Auxiliary constraint: Upper bound on emissions when minimizing cost'
;

cost_eq..                   of_cost =e= sum((N,T,U), (C_opr(U)+S*P(U))*g(N,T,U)) + sum((N,T), C_dsr(T,N)*dsr(N,T)) + sum((N,E), C_inv_vre(E)*b_vre(N,E)) + sum((N,U), C_ava_gen(U)*a_gen(N,U)) + sum((N,E), C_ava_vre(E)*a_vre(N,E)) + sum((N,W), C_ava_hyd(W)*a_hyd(N,W));
emission_eq..               of_emission =e= sum((N,T,U), P(U)*g(N,T,U));
cvar_eq..                   of_cvar =e= sum(N, (h(N) + (1/(1-alpha)*T)*sum(T, z(N,T))));
energy_balance(N,T)..       sum(U, g(N,T,U)) + sum(E, u(N,T,E)) + sum(W, Q_hyd(N,W)*r_out(N,T,W) - F_hyd(N,W)*r_in(N,T,W)) + sum(L $ NMinus(L,N), V*TT(T)*f(L,T)) - sum(L $ NPlus(L,N), V*TT(T)*f(L,T)) + dsr(N,T) + X(T,N) - D(T,N) =e= 0;
flow_pos_limit(L,T)..       TT(T)*K_pos(L)-V*TT(T)*f(L,T) =g= 0;
flow_neg_limit(L,T)..       V*TT(T)*f(L,T) + TT(T)*K_neg(L) =g= 0;
load_shed_limit(N,T)..      TT(T)*D_dsr(T,N) - dsr(N,T) =g= 0;
gen_capacity_limit(N,T,U).. TT(T)*a_gen(N,U) - g(N,T,U) =g= 0;
gen_avail_limit(N,U)..      G_gen(N,U) - a_gen(N,U) =g= 0;
vre_capacity_limit(N,T,E).. TT(T)*A(T,E,N)*a_vre(N,E) - u(N,T,E) =g= 0;
       


    


