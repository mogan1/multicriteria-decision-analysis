* Cost-minimizing LP for energy portfolio in case of SE multiple bidding zones - An MCDA implementation
* Author: Pushan Deb and Mousam Ganguly
*
* Solve for 36 points to prove Cost vs Emission effect on generation switch for a fixed CVar


$INCLUDE mcda_nordic_data_2024.gms
* ---- Scenario targets (Cost/Emission/CVaR tuples)
$INCLUDE pareto_targets.inc

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

* Targets for each run (assigned inside loop over scen)
Scalar
    C_target                'Fixed Cost target'
    E_target                'Fixed Emission target'
    Q_target                'Fixed CVaR target';

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
    uv(N,T,E)               'Variable: VRE generation level by unit e at node n in period T [MWh]'
    r_in(N,T,W)             'Variable: Volume of water pumped into hydro unit w at node n in period t [m3]'
    r_out(N,T,W)            'Variable: Volume of water turbined out from hydro unit w at node n in period t [m3]'
    r_sto(N,T,W)            'Variable: Volume of water stored by hydro unit w at node n in period t [m3]'
    sp(N,T,W)               'Variable: Volume of water spilled from hydro unit w at node n in period t [m3]'
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
;


cost_eq..     of_cost =e=
    sum((N,T,U), (C_opr(U)+S*P(U))*g(N,T,U))
  + sum((N,T),   C_dsr(T,N)*dsr(N,T))
  + sum((N,E),   C_inv_vre(E)*b_vre(N,E))
  + sum((N,U),   C_ava_gen(U)*a_gen(N,U))
  + sum((N,E),   C_ava_vre(E)*a_vre(N,E))
  + sum((N,W),   C_ava_hyd(W)*a_hyd(N,W));

emission_eq.. of_emission =e= sum((N,T,U), P(U)*g(N,T,U));

cvar_eq..     of_cvar =e= sum(N, (h(N) + (1/(1-alpha) * card(T)) * sum(T, z(N,T))));

energy_balance(N,T).. sum(U, g(N,T,U))
                    + sum(E, uv(N,T,E))
                    + sum(W, Q_hyd(N,W)*r_out(N,T,W) - F_hyd(N,W)*r_in(N,T,W))
                    + sum(L $ NMinus(L,N), V*TT(T)*f(L,T))
                    - sum(L $ NPlus(L,N),  V*TT(T)*f(L,T))
                    + dsr(N,T) + X(T,N) - D(T,N) =e= 0;

flow_pos_limit(L,T).. TT(T)*K_pos(L) - V*TT(T)*f(L,T) =g= 0;
flow_neg_limit(L,T).. V*TT(T)*f(L,T) + TT(T)*K_neg(L) =g= 0;

load_shed_limit(N,T).. TT(T)*D_dsr(T,N) - dsr(N,T) =g= 0;

gen_capacity_limit(N,T,U).. TT(T)*a_gen(N,U) - g(N,T,U) =g= 0;
gen_avail_limit(N,U)..      G_gen(N,U) - a_gen(N,U) =g= 0;

gen_up_ramp_limit(N,T,U)$( ORD(T) ge 2 )..
    TT(T)*R_up(U)*a_gen(N,U) - g(N,T,U) + g(N,T-1,U) =g= 0;

gen_down_ramp_limit(N,T,U)$( ORD(T) ge 2 )..
    g(N,T,U) - g(N,T-1,U) + TT(T)*R_down(U)*a_gen(N,U) =g= 0;

vre_capacity_limit(N,T,E).. TT(T)*A(T,E,N)*a_vre(N,E) - uv(N,T,E) =g= 0;
vre_avail_limit(N,E)..      G_vre(N,E) + b_vre(N,E) - a_vre(N,E) =g= 0;
vre_inv_limit(N,E)..        M(N,E)*G_vre(N,E) - b_vre(N,E) =g= 0;

storage_balance(N,T,W)..
    - r_sto(N,T,W)
    + r_sto(N,T-1,W)$(ORD(T) > 1)
    + RR_ini(N,W)$( ORD(T) eq 1 )
    + r_in(N,T,W) - r_out(N,T,W) - sp(N,T,W)
    + II_hourly_ror(N,T)$( ORD(W) eq 1)
    + II_hourly_res(N,T)$( ORD(W) eq 2)
    =e= 0;

storage_max_limit(N,T,W).. RR_max(N,W) - r_sto(N,T,W) =g= 0;
storage_min_limit(N,T,W)$(ORD(T) eq CARD(T)).. r_sto(N,T,W) - RR_min(N,W) =g= 0;

hydro_pump_limit(N,T,W).. TT(T)*RR_in(N,W)*RR_max(N,W) - r_in(N,T,W) =g= 0;
hydro_capacity_limit(N,T,W).. TT(T)*a_hyd(N,W) - Q_hyd(N,W)*r_out(N,T,W) =g= 0;
hydro_avail_limit(N,T,W).. Y_hyd(N,W) - a_hyd(N,W) =g= 0;

residual_load_balance(N,T).. z(N,T) + h(N) - sum(U, g(N,T,U)) =g= 0;

* Fixed-point equalities for a given scenario tuple
Equation
   cost_expr     'Total cost at fixed cost target'
   emission_expr 'Total emission at fixed emission target'
   cvar_expr     'Total CVaR at fixed CVaR target'
;

cost_expr..
    sum((N,T,U), (C_opr(U)+S*P(U))*g(N,T,U))
  + sum((N,T),   C_dsr(T,N)*dsr(N,T))
  + sum((N,E),   C_inv_vre(E)*b_vre(N,E))
  + sum((N,U),   C_ava_gen(U)*a_gen(N,U))
  + sum((N,E),   C_ava_vre(E)*a_vre(N,E))
  + sum((N,W),   C_ava_hyd(W)*a_hyd(N,W))
  =e= C_target;

emission_expr.. sum((N,T,U), P(U)*g(N,T,U)) =e= E_target;

cvar_expr.. sum(N, h(N) + (1/(1-alpha) * (1/card(T)) * sum(T, z(N,T)))) =e= Q_target;

* Dummy objective (fixed-point feasibility solve)
Variable dummy;
Equation dummy_obj;
dummy_obj.. dummy =e= 0;

* ---- Fixed point model (same as your report file)
model fixed_point /
    dummy_obj,
    energy_balance, flow_pos_limit, flow_neg_limit,
    load_shed_limit, gen_capacity_limit, gen_avail_limit,
    gen_up_ramp_limit, gen_down_ramp_limit,
    vre_capacity_limit, vre_avail_limit, vre_inv_limit,
    storage_balance, storage_max_limit, storage_min_limit,
    hydro_pump_limit, hydro_capacity_limit, hydro_avail_limit,
    residual_load_balance,
    cost_expr, emission_expr, cvar_expr
/;

option LP = CPLEX;
option optcr = 0.0000;
alpha = 0.95;

* ============================================================
* Scenario-indexed reporting (pattern from mcda_SE_2025_generation_from_pareto.gms L146+)
* ============================================================

Parameter
    CostAct(scen)               'Cost achieved (recomputed)'
    EmAct(scen)                 'Emissions achieved (recomputed)'
    CVaRAct(scen)               'CVaR achieved (recomputed)'
    SolveStat(scen)
    ModelStat(scen)

    Thermal_Node_Tech_run(scen,N,U) 'Thermal gen by node x thermal-tech [MWh]'
    VRE_Tech_run(scen,E)            'VRE gen by tech (all nodes, all time) [MWh]'
    Hyd_Node_Tech_run(scen,N,W)     'Hydro gen by node x hydro-tech [MWh]'
;

loop(scen,

    C_target = C_target_s(scen);
    E_target = E_target_s(scen);
    Q_target = Q_target_s(scen);

    solve fixed_point using lp minimizing dummy;

    SolveStat(scen) = fixed_point.solvestat;
    ModelStat(scen) = fixed_point.modelstat;

* IMPORTANT GUARD: only read .L when a solution exists
    if (fixed_point.solvestat = 1 and fixed_point.modelstat <= 2,

        CostAct(scen) = sum((N,T,U), (C_opr(U)+S*P(U))*g.L(N,T,U))
                      + sum((N,T),   C_dsr(T,N)*dsr.L(N,T))
                      + sum((N,E),   C_ava_vre(E)*a_vre.L(N,E))
                      + sum((N,W),   C_ava_hyd(W)*a_hyd.L(N,W));

        EmAct(scen)   = sum((N,T,U), P(U)*g.L(N,T,U));

        CVaRAct(scen) = sum(N, h.L(N) + (1/(1-alpha) * (1/card(T)) * sum(T, z.L(N,T))));

* Requested outputs
        Thermal_Node_Tech_run(scen,N,U) = sum(T, g.L(N,T,U));
        VRE_Tech_run(scen,E)            = sum((N,T), uv.L(N,T,E));
        Hyd_Node_Tech_run(scen,N,W)     = SUM(T, r_out.L(N,T,W));

    else
        CostAct(scen) = 0;
        EmAct(scen)   = 0;
        CVaRAct(scen) = 0;
        Thermal_Node_Tech_run(scen,N,U) = 0;
        VRE_Tech_run(scen,E)            = 0;
        Hyd_Node_Tech_run(scen,N,W)     = 0;
    );

);



* ============================================================
* GDX Output
* ============================================================
execute_unload  'output_generation_summary_36_points_2026_01_22.gdx'
    C_target_s, E_target_s, Q_target_s,
    CostAct, EmAct, CVaRAct,
    SolveStat, ModelStat,
    Thermal_Node_Tech_run,
    VRE_Tech_run,
    Hyd_Node_Tech_run;

* ============================================================
* CSV outputs (Run encoded as "C=<cost>|E=<em>|Q=<cvar>" no commas)
* ============================================================

File fth /'gen_node_tech_runs_36_points_2026_01_20.csv'/;;
put fth;
put 'Run,Node,Thermal_Tech,Generation(MWh)' /;
loop((scen,N,U),
    put 'C=', C_target_s(scen):0:0, '|E=', E_target_s(scen):0:2, '|Q=', Q_target_s(scen):0:2,
        ',', N.tl:0, ',', U.tl:0, ',', Thermal_Node_Tech_run(scen,N,U):18:6 /;
);
putclose fth;


File  fvre  /'gen_vre_tech_runs_36_points_points_2026_01_20.csv'/;
put fvre;
put 'Run,VRE_Tech,Generation(MWh)' /;
loop((scen,E),
    put 'C=', C_target_s(scen):0:0, '|E=', E_target_s(scen):0:2, '|Q=', Q_target_s(scen):0:2,
        ',', E.tl:0, ',', VRE_Tech_run(scen,E):18:6 /;
);
putclose fvre;

* diagnostics file
File fmap  /run_points_map_36_points_points_2026_01_20.csv/;
put fmap;
put 'scen,CostTarget,EmissionTarget,CVaRTarget,CostAct,EmAct,CVaRAct,SolveStat,ModelStat' /;
loop(scen,
    put scen.tl:0, ',',
        C_target_s(scen):0:0, ',', E_target_s(scen):0:2, ',', Q_target_s(scen):0:2, ',',
        CostAct(scen):0:2, ',', EmAct(scen):0:2, ',', CVaRAct(scen):0:2, ',',
        SolveStat(scen):0:0, ',', ModelStat(scen):0:0 /;
);
putclose fmap;
