# Urban Waste Collection & Routing Optimizer

A multi-stage Operations Research (OR) project designed to optimize urban waste collection logistics. This repository implements a series of solvers—ranging from classic heuristics to advanced Location-Routing Problems (LRP)—using **Google OR-Tools**, **Matplotlib**, and **Python**.

---

## 🚀 Project Overview

This project tackles the increasing pressure on urban waste management by optimizing vehicle routes and facility selection. It is structured into three progressive modules of increasing complexity:

### Phase 1: Basic CVRP Optimization
Solves the standard Capacitated Vehicle Routing Problem (CVRP) for 30 collection points.
* **Algorithm Chain**: Clarke-Wright Savings (Initial Solution) ➔ 2-opt Heuristic (Local Optimization) ➔ OR-Tools CP-SAT (Exact Global Search).
* **Key Metric**: Successfully reduced total daily travel distance from an initial heuristic result of **1157.19 km** to a global optimum of **1108.18 km**.
* **Real-world Robustness**: Includes an evaluator for **Manhattan Distance** to simulate grid-based urban blocks, resulting in a 1442 km path.

### Phase 2: Multi-Fleet & Heterogeneous Constraints
Addresses the complexity of four distinct waste types (Food, Recyclable, Other, and Hazardous) with a hybrid solving strategy.
* **High-Volume Waste**: Implements the OR-Tools Routing Model with dual constraints: **Vehicle Capacity** and **Maximum Travel Distance**.
* **Hazardous Waste**: Since volume is low, this is treated as a Pure TSP and solved using an exact **CP-SAT model with MTZ (Miller-Tucker-Zemlin)** subtour elimination.

### Phase 3: Location-Routing Problem (LRP) & Asymmetric Networks
The most advanced module, handling facility selection and non-standard road conditions.
* **Facility Location**: Uses CP-SAT to select the optimal subset of transfer stations based on construction costs and distance.
* **Asymmetric Routing**: Accounts for one-way streets and time-dependent congestion where $D(i,j) \neq D(j,i)$.
* **Visualization**: Features a directed-graph plotter to show the flow of garbage trucks through one-way urban networks.

---

## 🛠️ Technical Highlights

* **Hybrid Solving**: Combines the speed of meta-heuristics (Guided Local Search) with the precision of Constraint Programming.
* **Asymmetric Matrix Support**: Capable of modeling real-world traffic restrictions.
* **Dynamic Scaling**: Automatically handles float-to-integer scaling for high-precision CP-SAT solving.
* **Advanced Visualization**: Custom Matplotlib scripts for route animation, multi-fleet color coding, and directed pathing.

---

## 📂 Repository Structure

```text
.
├── data/                 # Raw Excel data from urban logistics sensors
├── docs/                 # Technical reports and mathematical derivations (PDF)
├── src/                  
│   ├── q1_cvrp_solver.py           # Basic CVRP & Heuristic chain
│   ├── q2_multitype_solver.py      # Multi-fleet & Hybrid strategy
│   └── q3_lrp_asymmetric_solver.py # LRP & Directed graph optimization
└── README.md