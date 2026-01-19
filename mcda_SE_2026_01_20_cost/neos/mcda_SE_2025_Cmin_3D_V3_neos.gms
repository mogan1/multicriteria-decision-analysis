* Cost-minimizing LP for energy portfolio in case of SE multiple bidding zones - An MCDA implementation
* Author: Pushan Deb and Mousam Ganguly

$offlisting
option limrow=0, limcol=0, solprint=off;
option LP = CPLEX;
option optcr = 0.0000;

$gdxin in

* =============================================================================
* STEP 1 — DECLARE EVERYTHING BEFORE USING $GDXIN OR $LOAD
* =============================================================================

Sets
    ii      "Firms"
    e       "VRE technologies"
    it      "Iterator"
    l       "Transmission lines"
    n       "Nodes"
    t       "Periods"
    u       "Thermal units"
    w       "Hydro units"
    dy      "Days of the year"

    ntou(n,u)   "Node-to-thermal"
    ntoe(n,e)   "Node-to-vre"
    ntow(n,w)   "Node-to-hydro"
;

Alias (e,ee), (n,nn), (ii,iii), (u,uu), (w,ww);

Parameters
    HY, HD, GWMW

    TT(t), TD(t,dy)

    K_pos(l), K_neg(l)
    NPlus(l,n), NMinus(l,n)

    D(t,n), C_dsr(t,n), D_max, D_dsr(t,n)
    X(t,n)

    P(u), R_up(u), R_down(u), C_opr(u)

    K_gen(ii,n,u), G_gen(n,u)
    K_vre(n,ii,e), G_vre(n,e)

    Y_hyd(n,w), Q_hyd(n,w)

    C_inv_vre(e), C_ava_vre(e), C_ava_hyd(w), C_ava_gen(u)

    A(t,e,n)

    F_hyd(n,w)

    R_max_table(ii,n,w), RR_max(n,w), RR_min(n,w)
    E_sto(n,w), RR_in(n,w), RR_ini(n,w)

    II_daily_ror(dy,n), II_hourly_ror(n,t), Tot_IN_ror(n)
    II_daily_res(dy,n), II_hourly_res(n,t), Tot_IN_res(n)

    I(n,t,w)

    M(n,e)

    V, S, W_C, W_E, Chat
    Cmin, Cmax, Emin, Emax
;


* 6x6 grid for cost, emissions and CVaR thresholds
Set
    egrid /e1*e6/
    qgrid /q1*q6/
    cgrid /c1*c6/
    ;


* Override known bounds
Cmax = 3.96093E+09;
Cmin = 3.14058E+09;
Emin = 98006.33;
Emax = 440714.77;


$load ii, e, it, l, n, t, u, w, dy, ntou, ntoe, ntow
$load HY, HD, GWMW, TT, TD
$load K_pos, K_neg, NPlus, NMinus
$load D, C_dsr, D_max, D_dsr, X
$load P, R_up, R_down, C_opr
$load K_gen, G_gen
$load K_vre, G_vre
$load Y_hyd, Q_hyd
$load C_inv_vre, C_ava_vre, C_ava_hyd, C_ava_gen
$load A, F_hyd
$load R_max_table, RR_max, RR_min, E_sto, RR_in, RR_ini
$load II_daily_ror, II_hourly_ror, Tot_IN_ror
$load II_daily_res, II_hourly_res, Tot_IN_res
$load I
$load V, S, W_C, W_E, Chat
$load M

$gdxin

Scalars
* Emin/Emax and Cmin/Cmax are defined in INCLUDE
    Qmin /7431.91/
    Qmax /10407.77/
    Qmin                    'Scalar variable: Minimized CVaR'
    Qmax                    'Scalar variable: Maximised CVaR'
    C2                      'Scalar variable: Intermediary cost during emission minimization'
    C3                      'Scalar variable: Intermediary cost during CVaR minimization'
    E1                      'Scalar variable: Intermediary emission during cost minimization'
    E3                      'Scalar variable: Intermediary emission during CVaR minimization'
    Q1                      'Scalar variable: Intermediary CVaR during cost minimization'
    Q2                      'Scalar variable: Intermediary CVaR during emission minimization'
    alpha                   /0.95/
    eps_em                  'emission cap used in each iteration'
    eps_q                   'cvar cap used in each iteration'
    eps_c                   'cost cap used in each iteration'
    Emin                    'Scalar variable: Minimized Emission'
    Emax                    'Scalar variable: Maximised Emission'
    Cmin                    'Scalar variable: Minimized Cost'
    Cmax                    'Scalar variable: Maximized Cost'
;

Parameters
    epsilon_em(egrid)      'Grid of emission caps'
    epsilon_cvar(qgrid)    'Grid of CVaR caps'
    epsilon_cost(cgrid)    'Grid of cost caps'
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
    cost_cap                'Limit for cost'

*    CostGoal                'Auxiliary constraint: Cost goal (used in minimax or benchmarking)'
*    EmissionGoal            'Auxiliary constraint: Emission goal (used in minimax or benchmarking)'
*    CostAux                 'Auxiliary constraint: Upper bound on cost when minimizing emissions'
*    CostAux2                'Auxiliary constraint: Cost control in policy-guided planning (PGP) emission minimization'
*    EmissionsAux            'Auxiliary constraint: Upper bound on emissions when minimizing cost'
;

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
emission_cap..                                  of_emission =l= eps_em;
cvar_cap..                                      of_cvar     =l= eps_q;
cost_cap..                                      of_cost     =l= eps_c;


model cost_lp / cost_eq, energy_balance, flow_pos_limit, flow_neg_limit, load_shed_limit, gen_capacity_limit, gen_avail_limit, gen_up_ramp_limit, gen_down_ramp_limit, vre_capacity_limit, vre_avail_limit, vre_inv_limit, storage_balance, storage_max_limit, storage_min_limit, hydro_pump_limit, hydro_capacity_limit, hydro_avail_limit, residual_load_balance /;
model emission_lp / emission_eq, energy_balance, flow_pos_limit, flow_neg_limit, load_shed_limit, gen_capacity_limit, gen_avail_limit, gen_up_ramp_limit, gen_down_ramp_limit, vre_capacity_limit, vre_avail_limit, vre_inv_limit, storage_balance, storage_max_limit, storage_min_limit, hydro_pump_limit, hydro_capacity_limit, hydro_avail_limit, residual_load_balance /;
model cvar_lp / cvar_eq, energy_balance, flow_pos_limit, flow_neg_limit, load_shed_limit, gen_capacity_limit, gen_avail_limit, gen_up_ramp_limit, gen_down_ramp_limit, vre_capacity_limit, vre_avail_limit, vre_inv_limit, storage_balance, storage_max_limit, storage_min_limit, hydro_pump_limit, hydro_capacity_limit, hydro_avail_limit, residual_load_balance /;
Model pareto3D_cost / cost_eq, emission_eq, cvar_eq,
                 energy_balance, flow_pos_limit, flow_neg_limit, load_shed_limit,
                 gen_capacity_limit, gen_avail_limit, gen_up_ramp_limit, gen_down_ramp_limit,
                 vre_capacity_limit, vre_avail_limit, vre_inv_limit,
                 storage_balance, storage_max_limit, storage_min_limit,
                 hydro_pump_limit, hydro_capacity_limit, hydro_avail_limit,
                 residual_load_balance,
                 emission_cap, cvar_cap /;
Model pareto3D_emission / cost_eq, emission_eq, cvar_eq,
                 energy_balance, flow_pos_limit, flow_neg_limit, load_shed_limit,
                 gen_capacity_limit, gen_avail_limit, gen_up_ramp_limit, gen_down_ramp_limit,
                 vre_capacity_limit, vre_avail_limit, vre_inv_limit,
                 storage_balance, storage_max_limit, storage_min_limit,
                 hydro_pump_limit, hydro_capacity_limit, hydro_avail_limit,
                 residual_load_balance,
                 cost_cap, cvar_cap /;
Model pareto3D_cvar / cost_eq, emission_eq, cvar_eq,
                 energy_balance, flow_pos_limit, flow_neg_limit, load_shed_limit,
                 gen_capacity_limit, gen_avail_limit, gen_up_ramp_limit, gen_down_ramp_limit,
                 vre_capacity_limit, vre_avail_limit, vre_inv_limit,
                 storage_balance, storage_max_limit, storage_min_limit,
                 hydro_pump_limit, hydro_capacity_limit, hydro_avail_limit,
                 residual_load_balance,
                 emission_cap, cost_cap /;



* Run cost minimization independently, record emission and cvar
*solve cost_lp using lp minimizing of_cost;

*Cmin = of_cost.L;
*E1 = sum((N,T,U), P(U)*g.L(N,T,U));
*Q1 = sum(N, (h.L(N) + (1/(1-alpha) * card(T)) * sum(T, z.L(N,T))));

* Run emission minimization independently, record cost and cvar
*solve emission_lp using lp minimizing of_emission;

*C2   = sum((N,T,U), (C_opr(U)+S*P(U))*g.L(N,T,U)) + sum((N,T), C_dsr(T,N)*dsr.L(N,T)) + sum((N,E), C_inv_vre(E)*b_vre.L(N,E)) + sum((N,U), C_ava_gen(U)*a_gen.L(N,U)) + sum((N,E), C_ava_vre(E)*a_vre.L(N,E)) + sum((N,W), C_ava_hyd(W)*a_hyd.L(N,W));
*Emin = of_emission.L;
*Q2 = sum(N, (h.L(N) + (1/(1-alpha) * card(T)) * sum(T, z.L(N,T))));

* Run CVaR minimization independently, record emission and cost
*solve cvar_lp using lp minimizing of_cvar;

*C3 = sum((N,T,U), (C_opr(U)+S*P(U))*g.L(N,T,U)) + sum((N,T), C_dsr(T,N)*dsr.L(N,T)) + sum((N,E), C_inv_vre(E)*b_vre.L(N,E)) + sum((N,U), C_ava_gen(U)*a_gen.L(N,U)) + sum((N,E), C_ava_vre(E)*a_vre.L(N,E)) + sum((N,W), C_ava_hyd(W)*a_hyd.L(N,W));;
*E3 = sum((N,T,U), P(U)*g.L(N,T,U));
*Qmin = of_cvar.L;

* Determine Max tuple (Cmax, Emax, Qmax)
*Cmax = max(C2, C3);
*Emax = max(E1, E3);
*Qmax = max(Q1, Q2);


** discretize the range of iteration between max and min of emission and CVar
**pre populate using loop
loop(egrid,
  epsilon_em(egrid) = Emin + (ord(egrid) - 1) * (Emax - Emin) / (card(egrid) - 1);
);
loop(qgrid,
  epsilon_cvar(qgrid) = Qmin + (ord(qgrid) - 1) * (Qmax - Qmin) / (card(qgrid) - 1);
);
loop(cgrid,
  epsilon_cost(cgrid) = Cmin + (ord(cgrid) - 1) * (Cmax - Cmin) / (card(cgrid) - 1);
);
*
*

file res /'results.txt'/;
put res;
Put '"epsilon_em","epsilon_cvar","cost","emission","cvar"' /;

loop(egrid,
  loop(qgrid,
    eps_q = epsilon_cvar(qgrid);
    eps_em  = epsilon_em(egrid);

    g.L(N,T,U)     = 0;
    dsr.L(N,T)     = 0;
    b_vre.L(N,E)   = 0;
    a_gen.L(N,U)   = 0;
    a_vre.L(N,E)   = 0;
    a_hyd.L(N,W)   = 0;
    z.L(N,T)       = 0;
    h.L(N)         = 0;
    uv.L(N,T,E)    = 0;
    r_in.L(N,T,W)  = 0;
    r_out.L(N,T,W) = 0;
    r_sto.L(N,T,W) = 0;
    sp.L(N,T,W)    = 0;
    f.L(L,T)       = 0;

    solve pareto3D_cost using lp minimizing of_cost;
*    
    put epsilon_em(egrid):12:2, ',', epsilon_cvar(qgrid):12:2, ',', of_cost.l:12:2, ',', of_emission.l:12:2, ',', of_cvar.l:12:2 /;
  );
);
putclose res;
