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
├── data/                  # Logistics data including collection point coordinates
├── docs/                  # Technical reports and mathematical derivations (PDF)
├── src/                   
│   ├── q1_cvrp_solver.py            # CVRP and heuristic implementation
│   ├── q2_multitype_solver.py       # Multi-fleet and hybrid strategy logic
│   └── q3_lrp_asymmetric_solver.py # LRP and directed graph optimization
└── README.md
```

---

# 城市垃圾收集与路径优化系统

本项目在运筹学（OR）框架下实现了城市垃圾收集物流的优化算法。项目利用 **Google OR-Tools**、**Matplotlib** 和 **Python**，探索了从经典启发式算法到选址-路径问题（LRP）的多种求解器。

## 项目概览

该实现通过优化车辆轨迹和设施选择，解决了城市垃圾管理中的物流挑战。项目结构分为三个递进模块：

### 1. 带容量限制的车辆路径问题 (CVRP)
本模块解决了针对 30 个收集点的标准 CVRP。
* **算法序列**：实现 Clarke-Wright 节约算法以获取初始解，随后采用 2-opt 局部搜索及 Google OR-Tools CP-SAT 求解器进行全局优化。
* **性能分析**：该实现将总行驶距离从启发式基准的 **1157.19 km** 缩减至计算得出的最优值 **1108.18 km**。
* **空间建模**：引入曼哈顿距离指标以模拟基于网格的城市街区结构，得出计算路径为 1442 km。

### 2. 异构约束与多车队管理
本模块通过混合求解策略处理四种不同类别垃圾（厨余、可回收、其他及有害垃圾）的配送。
* **高体量垃圾**：在车辆容量和最大作业距离的双重约束下，利用 OR-Tools 路径模型进行求解。
* **有害垃圾**：建模为纯旅行商问题（TSP），并使用包含 Miller-Tucker-Zemlin (MTZ) 子回路消除机制的 CP-SAT 模型进行处理。

### 3. 选址-路径问题 (LRP) 与非对称网络
第三个模块整合了设施选址和非均匀道路条件。
* **设施选址**：利用 CP-SAT 根据固定建设成本和运输距离识别中转站的最佳子集。
* **非对称路由**：考虑城市环境中 $D(i,j) \neq D(j,i)$ 的有向移动，模拟单行道和交通限制。
* **可视化**：实现有向图绘制器，用以展示非对称网络中的车辆流向。

## 技术特性

* **混合优化**：结合元启发式算法（引导局部搜索）与约束规划，以平衡计算效率和精度。
* **非对称矩阵支持**：能够建模有向交通限制。
* **数据缩放**：实现自动的浮点数转整数缩放，以便在 CP-SAT 环境中进行高精度求解。
* **可视化实现**：自定义脚本用于路径动画展示、多车队颜色编码以及有向路径表示。

## 📁 仓库结构

```text
.
├── data/                  # 物流数据，包括收集点坐标
├── docs/                  # 技术报告与数学推导 (PDF)
├── src/                   
│   ├── q1_cvrp_solver.py            # CVRP 与启发式算法实现
│   ├── q2_multitype_solver.py       # 多车队与混合策略逻辑
│   └── q3_lrp_asymmetric_solver.py # LRP 与有向图优化
└── README.md
```
