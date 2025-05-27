* Cost-minimizing LP for energy portfolio in case of SE multiple bidding zones - An MCDA implementation
* Author: Pushan Deb and Mousam Ganguly 

$INCLUDE mcda_nordic_data_2024.gms

Variables
    of_cost               'Objective function: Value for cost minimization'
    of_emission           'Objective function: Value for emission minimization'
    of_cvar               'Objective function: Value for CVaR of residual load minimization';
    
Positive Variables
    g(N,T,U)            'Variable: Generation of thermal unit u at node n at time t [MWh]'
    d(N,T)              'Variable: Demand-side response activated at node N at time T [MWh]'
    b_vre(N,E)          'Variable: Adopted capacity of VRE unit e at node n [MW]'
    ;

Equations
    cost_eq               'Objective function equation: Minimize total cost'
    emission_eq           'Objective function equation: Minimize total CO₂ emissions'
    cvar_eq               'Objective function equation: Minimize CVaR of residual load'
    
    energy_balance_eq     'Constraint: Energy balance (supply must meet demand)'
    flow_pos_limit        'Constraint: Power-flow limit in positive direction (inward) by line and by period'
    flow_neg_limit        'Constraint: Power-flow limit in negative direction (outward) by line and by period'
    load_shed_limit       'Constraint: Maximum allowable load shedding by node and time period'
    gen_capacity_limit    'Constraint: Thermal generation capped by available capacity for each technology and time'
    gen_avail_limit       'Constraint: Thermal generation must not exceed installed capacity per technology'
    vre_capacity_limit    'Constraint: VRE generation capped by available VRE output in each period'
    vre_avail_limit       'Constraint: VRE generation limited by installed capacity per technology'
    vre_inv_limit           'Constraint: Limit on new VRE adoption per technology'
    gen_up_ramp_limit        'Constraint: Thermal generation ramp-up limit between consecutive time periods'
    gen_down_ramp_limit      'Constraint: Thermal generation ramp-down limit between consecutive time periods'
    storage_balance        'Constraint: Storage state-of-charge balance for hydro units across time'
    storage_max_limit       'Constraint: Upper bound on storage capacity for hydro technologies by period'
    storage_min_limit       'Constraint: Lower bound on storage level for hydro technologies by period'
    hydro_pump_limit        'Constraint: Pumped-hydro charging limit per technology and period'
    hydro_capacity_limit    'Constraint: Hydro generation limit based on available water and capacity'
    hydro_avail_limit       'Constraint: Hydro generation restricted by installed capacity'
    
    CostGoal              'Auxiliary constraint: Cost goal (used in minimax or benchmarking)'
    EmissionGoal          'Auxiliary constraint: Emission goal (used in minimax or benchmarking)'
    CostAux               'Auxiliary constraint: Upper bound on cost when minimizing emissions'
    CostAux2              'Auxiliary constraint: Cost control in policy-guided planning (PGP) emission minimization'
    EmissionsAux          'Auxiliary constraint: Upper bound on emissions when minimizing cost';

cost_eq.. of_cost =e= sum((N,T,U), (C_opr(U)+S*P(U))*g(N,T,U)) + sum((N,T), C_dsr(T,N)*d(N,T)) + sum((N,E), C_inv_vre(E)*b_vre(N,E)) + sum((N,U), C_ava_gen(U)*a_gen(N,U))


    


