# Urban Waste Collection and Routing Optimization System

This repository provides an implementation of optimization algorithms for urban waste collection logistics within an Operations Research (OR) framework. The project explores a variety of solvers—ranging from classical heuristics to the Location-Routing Problem (LRP)—utilizing **Google OR-Tools**, **Matplotlib**, and **Python**.

## Project Overview

This implementation addresses logistical challenges in urban waste management by optimizing vehicle trajectories and facility selection. The project is structured into three progressive modules:

### 1. Capacitated Vehicle Routing Problem (CVRP)
This module addresses the standard CVRP for 30 collection points.
* **Algorithmic Sequence**: Implementation of the Clarke-Wright Savings heuristic for initial solutions, followed by 2-opt local search and the Google OR-Tools CP-SAT solver for global optimization.
* **Performance Analysis**: The implementation demonstrates a reduction in total travel distance from a heuristic baseline of **1157.19 km** to a calculated optimum of **1108.18 km**.
* **Spatial Modeling**: Incorporates Manhattan distance metrics to approximate grid-based urban block structures, resulting in a calculated path of 1442 km.

### 2. Heterogeneous Constraints and Multi-Fleet Management
This module addresses the distribution of four distinct waste categories (Food, Recyclable, Other, and Hazardous) through a hybrid solving strategy.
* **High-Volume Waste**: Utilizes the OR-Tools Routing Model under dual constraints of vehicle capacity and maximum operational distance.
* **Hazardous Waste**: Modeled as a Pure Traveling Salesman Problem (TSP) and addressed using a CP-SAT model incorporating Miller-Tucker-Zemlin (MTZ) subtour elimination.

### 3. Location-Routing Problem (LRP) and Asymmetric Networks
The third module incorporates facility selection and non-uniform road conditions.
* **Facility Location**: Employs CP-SAT to identify an optimal subset of transfer stations based on fixed construction costs and transportation distances.
* **Asymmetric Routing**: Accounts for directed movement in urban environments where $D(i,j) \neq D(j,i)$, simulating one-way streets and traffic constraints.
* **Visualization**: Implementation of a directed-graph plotter to illustrate vehicle flow within asymmetric networks.

## Technical Characteristics

* **Hybrid Optimization**: Integrates meta-heuristics (Guided Local Search) with Constraint Programming to balance computational efficiency and precision.
* **Asymmetric Matrix Support**: Capable of modeling directed traffic restrictions.
* **Data Scaling**: Implements automated float-to-integer scaling to facilitate high-precision solving within the CP-SAT environment.
* **Visualization Implementation**: Custom scripts for route animation, multi-fleet color coding, and directed path representation.

## 📁 Repository Structure

```text
.
├── data/                 # Logistics data including collection point coordinates
├── docs/                 # Technical reports and mathematical derivations (PDF)
├── src/                  
│   ├── q1_cvrp_solver.py           # CVRP and heuristic implementation
│   ├── q2_multitype_solver.py      # Multi-fleet and hybrid strategy logic
│   └── q3_lrp_asymmetric_solver.py # LRP and directed graph optimization
└── README.md
