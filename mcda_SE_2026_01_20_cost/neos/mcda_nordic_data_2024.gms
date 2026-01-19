* mcda_nordic_data_2024.gms: Nordic data for multi-zone MCDA problem instances

* Set definitions

Sets
ii Firms                          /vattenfall, eon, okg, fortum, statkraft, sydkraft, skellefteakraft, fringe_se/
e VRE technologies               /wind, solar/
it Loops through minimax         /it1*it5/
l Transmission lines             /l1*l3/
n Nodes                          /se1, se2, se3, se4/
t Periods                        /t1*t8760/
u Thermal units                  /gas, oil, nuclear, chp_coal, chp_waste, chp_gas, chp_peat, chp_biomass/
w Hydro units                    /ror,  res, phs/
dy days of the year                /d1*d365/
*dh hours of a day                 /h1*h24/
*dhh(dd,dh)  each day has 24 h      / #dd.#dh    /
*tdd(t,dd) map hours of year to days of year / #t:#dd   /
ntou(N,U) Node-to-thermal relation  /se4.gas, se3.oil, se4.oil, se3.nuclear, se3.chp_coal, se3.chp_waste, se3.chp_gas, se3.chp_peat, se1.chp_biomass, se3.chp_biomass/
ntoe(N,E) Node-to-VRE relation  /se1.wind, se2.wind, se3.wind, se4.wind, se1.solar, se2.solar, se3.solar, se4.solar/
ntow(N,W) Node-to-hydro relation  /se1.res, se2.res, se3.res, se4.res, se3.ror, se4.ror/
;

Alias(e,ee);
Alias(n,nn);
Alias(ii,iii);
Alias(u,uu);
Alias(w,ww);

Scalars
HY number of hours in a year /8760/
HD number of hours in a day /24/
GWMW convert GW to MW /1000/
;

* Time parameters

Parameter TT(T) Duration of each period in hours;
TT(T) = 1;

Parameter
TD(T,DY)   Mapping from hours to days;

TD(T,DY)= 1$(ORD(T)<ORD(DY)*HD+1 AND ORD(T)>ORD(DY)*HD-HD);

Display TD;

* Transmission-line parameters

Parameters
K_pos(L)      Maximum rated capacity of each line in positive direction
/ l1      3300
  l2      7300
  l3      5400
/

K_neg(L)      Maximum rated capacity of each line in negative direction
/ l1      3300
  l2      7300
  l3      2000
/
;

Parameters
NPlus(L,N)
/
l1.se1 1
l2.se2 1
l3.se3 1
/

NMinus(L,N)
/
l1.se2 1
l2.se3 1
l3.se4 1
/
;

* Demand parameters
*$CALL GDXXRW.EXE consumption_se_2018.xlsx par=D rdim=1 cdim=1 rng=A1:E8761 trace=3

Parameter D(T,N) Nodal demand by period ;
$GDXIN consumption_se_2018.gdx
$LOAD D
$GDXIN

Parameter C_dsr(T,N)  Cost of nodal demand-side response by period;
C_dsr(T,N)=4000;

Parameter D_max Maximum proportion of DSR allowed;
D_max = 0.1;

Parameter D_dsr(T,N) Maximum DSR allowed by node and by period;
D_dsr(T,N) = D(T,N)*D_max;

display D, C_dsr, D_max, D_dsr;

* Import parameters

* $CALL GDXXRW.EXE netimport_se_2018.xlsx par=X rdim=1 cdim=1 rng=A1:E8761 trace=3

Parameter X(T,N)  Nodal net import by period ; ;
$GDXIN netimport_se_2018
$LOAD X
$GDXIN

Display X;

***********************************************************
***************************** Power generation parameters *
***********************************************************


Parameter P(U) Emission rate of thermal unit u
/
gas       0.5
oil       0.72
nuclear   0
chp_coal  0.83
chp_waste 0.94
chp_gas   0.5
chp_peat  1.09
chp_biomass 0
/
;

Parameter R_up(U)  Ramp-up limit of thermal unit u as a share of total installed capacity
/
gas      0.5
oil      0.7
nuclear  0.1
chp_coal    0.1
chp_waste   0.1
chp_gas   0.5
chp_peat   0.1
chp_biomass 0.1
/
;

Parameter R_down(U)  Ramp-down limit of thermal unit u as a share of total installed capacity
/
gas         0.5
oil         0.7
nuclear     0.1
chp_coal   0.1
chp_waste  0.1
chp_gas    0.1
chp_peat   0.1
chp_biomass   0.1
/
;

Parameter C_opr(U) Marginal cost of production of thermal unit u
/
gas       65
oil       67
nuclear   21
chp_coal  37.3
chp_waste 22
chp_gas   56.8
chp_peat  21.8
chp_biomass 26.7
/
;

Table K_gen(II,N,U) Installed conventional generation capacity of generation type u at node n owned by firm ii
                              gas            oil      nuclear       chp_coal       chp_waste       chp_gas        chp_peat         chp_biomass
vattenfall.se3                                        4939.8                       130                                             58.7
eon.se1                                                                                                                            65.4
eon.se3                                      763      656.4
okg.se3                                               826.5
fortum.se3                                            1354.33                      148.3                           165             128
fringe_se.se3                                                       130                             260                            608
fringe_se.se4                  448           680
;

Parameter G_gen(N,U) Installed generation capacity of thermal unit u at node n;
G_gen(N,U)=SUM(II, K_gen(II,N,U));

Display K_gen, G_gen;

Table K_vre(N,II,E)  VRE generation capacity type e of firm ii at node n
                      wind         solar
se1.eon               237.8
se1.fringe_se         428.2        1.391
se2.vattenfall        94.5
se2.fortum            75
se2.statkraft         547
se2.fringe_se         134.5        7.606
se3.fringe_se         2170         99.589
se4.vattenfall        166
se4.fringe_se         1396         41.428
;

Parameter G_vre(N,E) Installed generation capacity of VRE unit e at node n;
G_vre(N,E)=SUM(II, K_vre(N,II,E));

Display K_vre, G_vre;

Table Y_hyd(N,W) Installed hydro generation capacity of generation type w at node n

        ror       res        
se1               5370.57
se2               8165
se3     400       2190
se4     100       244
;

Display Y_hyd;

Parameter Q_hyd(N,W) Generation efficiency of hydro unit w at node n;

Q_hyd(N,W) = 1;

Parameter C_inv_vre(E) Amortized per-MW investment cost of VRE unit e per annum
/
wind      81848
solar     85176
/
;

Parameter C_ava_vre(E) Amortized per-MW O&M cost of VRE unit e per annum
/
wind      26260
solar     15184
/
;

Parameter C_ava_hyd(W) Amortized per-MW O&M cost of hydro unit w per annum
/
ror                      0
res                      29796
/
;

Parameter C_ava_gen(U) Amortized per-MW O&M cost of thermal unit u per annum
/
gas            6968
oil           14092
nuclear      121316
chp_coal      40456
chp_waste     14092
chp_gas        6968
chp_peat      14092
chp_biomass   14092
/
;

Display C_inv_vre, C_ava_vre, C_ava_hyd, C_ava_gen;

* $CALL GDXXRW.EXE vre_se_2018.xlsx par=A rdim=2 cdim=1 rng=A1:F17521 trace=3

Parameter A(T,E,N) Nodal VRE availability factor by period ;
$GDXIN vre_se_2018.gdx
$LOAD A
$GDXIN

Display A;

Parameter F_hyd(N,W) Pumped-hydro efficiency of hydro unit w at node n;

F_hyd(N,W)$(ntow(N,W) AND ORD(W)=3) = 1.36;

Table R_max_table(II,N,W) Table for storage capacity
                    ror   res
vattenfall.se1            12210.36513
vattenfall.se2            3997.632655
vattenfall.se3            670.1296548
fortum.se2                4605.050604
fortum.se3                1347.155496
statkraft.se1             665.2928729
statkraft.se2             1517.553091
statkraft.se4             350.2718
sydkraft.se2              3905.438404
sydkraft.se4              199.7282
skellefteakraft.se1       1625.505084
fringe_se.se2             2074.325245
fringe_se.se3             382.7148494
;

Parameter RR_max(N,W) Maximum reservoir volume of hydro unit w at node n;
RR_max(N,W)=SUM(II, R_max_table(II,N,W))*GWMW;

Parameter RR_min(N,W) Minimum reservoir volume of hydro unit w at node n;

RR_min(N,W)$(ntow(N,W) AND ORD(W)=2) = 0.53*RR_max(N,W);
RR_min(N,W)$(ntow(N,W) AND ORD(W)=3) = 0.53*RR_max(N,W);

Parameter E_sto(N,W) Self-discharge rate of hydro unit w at node n;

E_sto(N,W)$(ntow(N,W) AND ORD(W)=2) = 0;
E_sto(N,W)$(ntow(N,W) AND ORD(W)=3) = 0;

Parameter RR_in(N,W) Maximum charging rate of hydro unit w at node n;

RR_in(N,W)$(ntow(N,W) AND ORD(W)=3) = 0.20;

Parameter RR_ini(N,W)  Initial storage level of hydro unit w at node n;

RR_ini(N,W)$(ntow(N,W) AND ORD(W)=2) = 0.6*RR_max(N,W);
RR_ini(N,W)$(ntow(N,W) AND ORD(W)=3) = 0.6*RR_max(N,W);

Display Q_hyd, F_hyd, RR_max, RR_min, E_sto, RR_in, RR_ini;


* $CALL GDXXRW.EXE inflow_ror_se_2018.xlsx par=II_daily_ror rdim=1 cdim=1 rng=A1:E366 trace=3

Parameter II_daily_ror(DY,N) Natural daily inflow of ror plants;
$GDXIN inflow_ror_se_2018.gdx
$LOAD II_daily_ror
$GDXIN

Parameter II_hourly_ror(N,T) Natural hourly inflow of ror plants;

II_hourly_ror(N,T)=SUM((DY)$TD(T,DY), II_daily_ror(DY,N))*(1/HD)*GWMW;
II_hourly_ror(N,T)$(II_hourly_ror(N,T) <0 )=-II_hourly_ror(N,T);

Display II_daily_ror, II_hourly_ror;

Parameter Tot_IN_ror(N) Total annual natural inflow of ror plants;

Tot_IN_ror(N)=SUM(T, II_hourly_ror(N,T));

Display Tot_IN_ror;

* $CALL GDXXRW.EXE inflow_res_se_2018.xlsx par=II_daily_res rdim=1 cdim=1 rng=A1:E366 trace=3

Parameter II_daily_res(DY,N) Natural daily inflow of res plants;
$GDXIN inflow_res_se_2018.gdx
$LOAD II_daily_res
$GDXIN

Parameter II_hourly_res(N,T) Natural hourly inflow of res plants;

II_hourly_res(N,T)=SUM((DY)$TD(T,DY), II_daily_res(DY,N))*(1/HD)*GWMW;
II_hourly_res(N,T)$(II_hourly_res(N,T) <0 )=-II_hourly_res(N,T);

Display II_daily_res, II_hourly_res;

Parameter Tot_IN_res(N) Total annual natural inflow of res plants;

Tot_IN_res(N)=SUM(T, II_hourly_res(N,T));

Display Tot_IN_res;

Parameter I(N,T,W)   Natural inflow to hydro unit w at node n in period t;

I(N,T,'ror')=II_hourly_ror(N,T);
I(N,T,'res')=II_hourly_res(N,T);

Display I;

* Scaling factor for power flow

Scalar V /1/;

* Price of emissions

Scalar S /15/;

* Weight on cost minimization

Scalar W_C /0/;

* Minimized cost

Scalar Cmin /0/;

* Emissions associated with minimized cost

Scalar Emax /0/;

* Weight on emission minimization

Scalar W_E /1/;

* Minimized emissions

Scalar Emin /0/;

* Cost associated with minimized emissions

Scalar Cmax /0/;

* Cost steps used in PGP

Scalar Chat /0/;

* Maximum adoption of VRE capacity
Table
M(N,E)   Maximum adoption factor for VRE unit e at node n
         wind   solar
se1       1      1        
se2       1      1
se3       1      1     
se4       1      1          
;

execute_unload 'in.gdx';