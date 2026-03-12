# Energy Portfolio Optimization: Cost, Emission, and Risk (CVaR)

![GAMS](https://img.shields.io/badge/Language-GAMS-blue)
![Solver](https://img.shields.io/badge/Solver-CPLEX-orange)

This repository contains the GAMS implementation of a **Multi-Criteria Decision Analysis (MCDA)** framework for energy portfolio optimization. The model is designed for a multi-bidding zone system (specifically the Nordic/Swedish market) and utilizes Linear Programming (LP) to balance economic, environmental, and risk-based objectives.

## 📖 Overview

The model optimizes an energy portfolio by simultaneously considering three competing objective functions:
1.  **Cost Minimization:** Minimizes operational costs, emission taxes, investment in VRE, and capacity availability.
2.  **Emission Minimization:** Minimizes total $CO_2$ footprint from thermal generation.
3.  **Risk Minimization:** Minimizes the **Conditional Value at Risk (CVaR)** of the residual load at a $95\%$ confidence level.

The solution space is explored using the **Epsilon-Constraint Method**, generating a Pareto frontier to visualize trade-offs between cost, carbon, and system risk.

## 🛠 Features

* **Multi-Zone Modeling:** Handles power flows across transmission lines ($L$) and nodes ($N$) with thermal, hydro, and VRE technologies.
* **Hydro Management:** Includes reservoir balance, pumped-hydro logic, and spilled water variables.
* **Demand-Side Response (DSR):** Integrated DSR activation and utilization tracking.
* **Advanced Risk Metrics:** Implements CVaR to handle "tail risk," ensuring grid reliability under extreme scenarios.
* **Boundary Correction:** Includes a **CVaR Sweep** script to handle spurious/degenerate values by applying tighter physical bounds ($C_{ub}$ and $E_{lb}$).

## 📂 Project Structure

* `mcda_nordic_data_2024.gms`: Core data file containing sets (nodes, time, technology) and parameters (demand, VRE profiles, costs).
* `portfolio_optimization.gms`: Primary optimization script containing the LP formulation and epsilon-constraint loop.
* `cvar_sweep.gms`: Specialized script for refining the Pareto boundary when CVaR minimization yields degenerate solutions.

## 🚀 Getting Started

### Prerequisites
* **GAMS Distribution:** Version 30.x or higher recommended.
* **Solver:** An LP solver (the code is configured for **CPLEX**).
* **GDX Viewer:** Recommended to inspect the `output_dsr_summary.gdx` results.
* **NEOS server:** Internet-based service that provides access to optimization solvers.

### Execution
1.  Ensure `mcda_nordic_data_2024.gms` is in the working directory.
2.  Run the primary optimization:
    ```bash
    gams mcda_SE_2025_Cmin_3D.gms lp=cplex
    ```
3.  To address spurious values in risk minimization, run the sweep:
    ```bash
    gams mcda_SE_2025_Cvar_sweep.gms lp=cplex
    ```
4.  To get the 36 points for Pareto frtontier, run the following in NEOS server:
    ```bash
    gams mcda_SE_2025_Cmin_3D_V3_neos.gms lp=cplex
    ```
5.  To get the generation mix for the 36 points, run the following:
    ```bash
    gams mcda_SE_2025_report_summary_36points.gms lp=cplex
    ```
6. Execute the *optimization.py* file to plot the images for the analysis

## 📝 Authors
* **Pushan Deb**
* **Mousam Ganguly**

---
*This work was developed as part of a thesis/research project focusing on Swedish energy market integration and MCDA implementation. It is based out of data compiled by Hasanzadeh et al. (2023) on the Nordic market*