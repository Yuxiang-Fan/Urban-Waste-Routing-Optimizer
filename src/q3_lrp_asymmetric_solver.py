"""
Urban Waste Collection VRP Optimizer - Question 3 (Ultimate LRP Edition)
Author: Yuxiang Fan
Description: Solves the Capacitated Location-Routing Problem (CLRP) with Asymmetric Networks.
Stage 1: Facility Location (Transfer Station Selection) via CP-SAT.
Stage 2: Asymmetric CVRP (Directed Graph) via OR-Tools Routing Model.
Includes directed arrow visualization for asymmetric paths.
"""

import math
import time
import random
import matplotlib.pyplot as plt
from ortools.sat.python import cp_model
from ortools.constraint_solver import pywrapcp, routing_enums_pb2

plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False

class UrbanWasteLRPOptimizer:
    def __init__(self, seed=42):
        self.scale_waste = 10     
        self.scale_dist = 1000    
        random.seed(seed) # 固定随机种子以保证非对称路网的可复现性
        
        # 1. 垃圾收集点集合 (Collection Points)
        self.collection_points = [
            {"id": 1, "x": 12, "y": 8, "waste": 1.2}, {"id": 2, "x": 5, "y": 15, "waste": 2.3},
            {"id": 3, "x": 20, "y": 30, "waste": 1.8}, {"id": 4, "x": 25, "y": 10, "waste": 3.1},
            {"id": 5, "x": 35, "y": 22, "waste": 2.7}, {"id": 6, "x": 18, "y": 5, "waste": 1.5},
            {"id": 7, "x": 30, "y": 35, "waste": 2.9}, {"id": 8, "x": 10, "y": 25, "waste": 1.1},
            {"id": 9, "x": 22, "y": 18, "waste": 2.4}, {"id": 10, "x": 38, "y": 15, "waste": 3.0},
            {"id": 11, "x": 5, "y": 8, "waste": 1.7}, {"id": 12, "x": 15, "y": 32, "waste": 2.1},
            {"id": 13, "x": 28, "y": 5, "waste": 3.2}, {"id": 14, "x": 30, "y": 12, "waste": 2.6},
            {"id": 15, "x": 10, "y": 10, "waste": 1.9}, {"id": 16, "x": 20, "y": 20, "waste": 2.5},
            {"id": 17, "x": 35, "y": 30, "waste": 3.3}, {"id": 18, "x": 8, "y": 22, "waste": 1.3},
            {"id": 19, "x": 25, "y": 25, "waste": 2.8}, {"id": 20, "x": 32, "y": 8, "waste": 3.4},
            {"id": 21, "x": 15, "y": 5, "waste": 1.6}, {"id": 22, "x": 28, "y": 20, "waste": 2.2},
            {"id": 23, "x": 38, "y": 25, "waste": 3.5}, {"id": 24, "x": 10, "y": 30, "waste": 1.4},
            {"id": 25, "x": 20, "y": 10, "waste": 2.0}, {"id": 26, "x": 30, "y": 18, "waste": 3.6},
            {"id": 27, "x": 5, "y": 25, "waste": 1.0}, {"id": 28, "x": 18, "y": 30, "waste": 2.3},
            {"id": 29, "x": 35, "y": 10, "waste": 3.7}, {"id": 30, "x": 22, "y": 35, "waste": 1.9}
        ]
        
        # 2. 候选处理中转站 (Candidate Transfer Stations)
        self.candidate_stations = {
            101: {"x": 50, "y": 50, "capacity": 25, "cost": 500},
            102: {"x": 60, "y": 60, "capacity": 30, "cost": 600},
            103: {"x": 70, "y": 70, "capacity": 20, "cost": 450},
            104: {"x": 80, "y": 80, "capacity": 35, "cost": 700},
            105: {"x": 90, "y": 90, "capacity": 15, "cost": 400}
        }
        
        # 建立全局索引字典
        self.all_nodes = {p["id"]: {"x": p["x"], "y": p["y"], "waste": p["waste"]} for p in self.collection_points}
        for ts_id, ts_data in self.candidate_stations.items():
            self.all_nodes[ts_id] = {"x": ts_data["x"], "y": ts_data["y"], "waste": 0}
            
        self.asymmetric_dist_matrix = self._build_asymmetric_matrix()

    def _build_asymmetric_matrix(self):
        """
        构建非对称距离矩阵。
        模拟城市道路中的单行道或拥堵：D(i,j) 不一定等于 D(j,i)。
        """
        dist_matrix = {}
        nodes = list(self.all_nodes.keys())
        for i in nodes:
            for j in nodes:
                if i == j:
                    dist_matrix[(i, j)] = 0.0
                else:
                    base_dist = math.hypot(self.all_nodes[i]["x"] - self.all_nodes[j]["x"], 
                                           self.all_nodes[i]["y"] - self.all_nodes[j]["y"])
                    # 模拟 10% 的路段存在单向拥堵或绕行，距离增加 30%~80%
                    asymmetry_factor = 1.0
                    if random.random() < 0.1:
                        asymmetry_factor = random.uniform(1.3, 1.8)
                    dist_matrix[(i, j)] = base_dist * asymmetry_factor
        return dist_matrix

    # ================== 第一阶段：CP-SAT 设施选址 (Facility Location) ==================
    def select_transfer_stations(self):
        """
        根据各中转站容量和建造成本，将垃圾收集点分配给最优的中转站组合。
        返回: 选中的中转站列表，以及分配给各中转站的收集点列表。
        """
        print("\n[*] 正在运行第一阶段：CP-SAT 容量约束设施选址模型...")
        model = cp_model.CpModel()
        
        points = [p["id"] for p in self.collection_points]
        stations = list(self.candidate_stations.keys())
        
        # 决策变量：y[s] 表示是否启用中转站 s
        y = {s: model.NewBoolVar(f"y_{s}") for s in stations}
        # 决策变量：x[p, s] 表示收集点 p 是否分配给中转站 s
        x = {(p, s): model.NewBoolVar(f"x_{p}_{s}") for p in points for s in stations}
        
        # 1. 约束：每个收集点必须且只能分配给一个启用的中转站
        for p in points:
            model.AddExactlyOne([x[(p, s)] for s in stations])
            for s in stations:
                model.AddImplication(x[(p, s)], y[s])
                
        # 2. 约束：中转站容量限制
        for s in stations:
            capacity_int = int(self.candidate_stations[s]["capacity"] * self.scale_waste)
            model.Add(sum(int(self.all_nodes[p]["waste"] * self.scale_waste) * x[(p, s)] for p in points) <= capacity_int * y[s])
            
        # 3. 目标：最小化距离成本 + 建设成本
        dist_cost = sum(int(self.asymmetric_dist_matrix[(p, s)] * self.scale_dist) * x[(p, s)] 
                        for p in points for s in stations)
        build_cost = sum(int(self.candidate_stations[s]["cost"] * self.scale_dist) * y[s] for s in stations)
        model.Minimize(dist_cost + build_cost)
        
        solver = cp_model.CpSolver()
        solver.parameters.max_time_in_seconds = 10
        status = solver.Solve(model)
        
        if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
            active_stations = [s for s in stations if solver.Value(y[s]) == 1]
            assignments = {s: [] for s in active_stations}
            for p in points:
                for s in active_stations:
                    if solver.Value(x[(p, s)]) == 1:
                        assignments[s].append(p)
            print(f"[+] 选址完成！启用的中转站: {active_stations}")
            for s, assigned_pts in assignments.items():
                load = sum(self.all_nodes[p]["waste"] for p in assigned_pts)
                print(f"    -> 中转站 {s} (负载: {load:.2f}/{self.candidate_stations[s]['capacity']}吨) 负责节点: {assigned_pts}")
            return assignments
        else:
            print("[-] 选址失败，检查约束。")
            return None

    # ================== 第二阶段：非对称路网 VRP (Asymmetric CVRP) ==================
    def solve_asymmetric_routing_for_station(self, station_id, assigned_points, vehicle_capacity=8):
        """
        为某个激活的中转站及其负责的节点，在非对称路网上规划车队路线。
        """
        # 构建局部映射（OR-Tools 需要从 0 开始连续的索引）
        local_nodes = [station_id] + assigned_points
        num_local_nodes = len(local_nodes)
        
        # 构建局部非对称距离矩阵
        local_dist_matrix = [[int(self.asymmetric_dist_matrix[(i, j)] * self.scale_dist) for j in local_nodes] for i in local_nodes]
        local_demands = [0] + [int(self.all_nodes[p]["waste"] * self.scale_dist) for p in assigned_points]
        
        num_vehicles = len(assigned_points) # 最大允许车辆数
        manager = pywrapcp.RoutingIndexManager(num_local_nodes, num_vehicles, 0)
        routing = pywrapcp.RoutingModel(manager)

        def distance_callback(from_index, to_index):
            return local_dist_matrix[manager.IndexToNode(from_index)][manager.IndexToNode(to_index)]
        transit_callback_index = routing.RegisterTransitCallback(distance_callback)
        routing.SetArcCostEvaluatorOfAllVehicles(transit_callback_index)

        def demand_callback(from_index):
            return local_demands[manager.IndexToNode(from_index)]
        demand_callback_index = routing.RegisterUnaryTransitCallback(demand_callback)
        routing.AddDimensionWithVehicleCapacity(
            demand_callback_index, 0, [int(vehicle_capacity * self.scale_dist)] * num_vehicles, True, 'Capacity'
        )

        search_parameters = pywrapcp.DefaultRoutingSearchParameters()
        search_parameters.first_solution_strategy = routing_enums_pb2.FirstSolutionStrategy.PATH_CHEAPEST_ARC
        search_parameters.local_search_metaheuristic = routing_enums_pb2.LocalSearchMetaheuristic.GUIDED_LOCAL_SEARCH
        search_parameters.time_limit.FromSeconds(3) # 每个中转站分配 3 秒求解

        solution = routing.SolveWithParameters(search_parameters)

        if solution:
            routes = []
            for vehicle_id in range(num_vehicles):
                index = routing.Start(vehicle_id)
                route = []
                while not routing.IsEnd(index):
                    node_idx = manager.IndexToNode(index)
                    route.append(local_nodes[node_idx])
                    index = solution.Value(routing.NextVar(index))
                route.append(local_nodes[manager.IndexToNode(index)]) # 回到中转站
                if len(route) > 2:
                    routes.append(route)
            return routes
        return []

    # ================== 箭头可视化模块 (Directed Visualization) ==================
    def plot_lrp_directed_network(self, all_routes, activated_stations):
        """绘制包含带方向箭头的非对称路径和中转站状态的高级路网图"""
        plt.figure(figsize=(14, 10))
        
        # 1. 绘制普通收集点
        x_pts = [p["x"] for p in self.collection_points]
        y_pts = [p["y"] for p in self.collection_points]
        plt.scatter(x_pts, y_pts, color='royalblue', s=40, zorder=3, label='收集点')
        for p in self.collection_points:
            plt.text(p["x"]+0.5, p["y"]+0.5, str(p["id"]), fontsize=8)

        # 2. 绘制中转站 (区分未启用和已启用)
        for ts_id, ts_data in self.candidate_stations.items():
            if ts_id in activated_stations:
                plt.scatter(ts_data["x"], ts_data["y"], color='red', marker='^', s=250, zorder=5, edgecolors='black')
                plt.text(ts_data["x"]+1, ts_data["y"]-1, f"HUB-{ts_id}", fontweight='bold', color='red')
            else:
                plt.scatter(ts_data["x"], ts_data["y"], color='lightgray', marker='^', s=150, zorder=2)
                
        # 增加图例代理
        plt.scatter([], [], color='red', marker='^', s=150, edgecolors='black', label='启用中转站')
        plt.scatter([], [], color='lightgray', marker='^', s=150, label='未启用候选站')

        # 3. 绘制带箭头的有向路径 (展示非对称特征)
        colors = ['#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
        color_idx = 0
        
        for ts_id, routes in all_routes.items():
            for route in routes:
                c = colors[color_idx % len(colors)]
                for i in range(len(route) - 1):
                    n1, n2 = route[i], route[i+1]
                    x1, y1 = self.all_nodes[n1]["x"], self.all_nodes[n1]["y"]
                    x2, y2 = self.all_nodes[n2]["x"], self.all_nodes[n2]["y"]
                    # 使用 annotate 画箭头
                    plt.annotate('', xy=(x2, y2), xytext=(x1, y1),
                                 arrowprops=dict(arrowstyle="->", color=c, lw=1.5, alpha=0.8))
                color_idx += 1

        plt.title("问题三：多中转站选址与非对称路网有向路径规划 (LRP)", fontsize=16, fontweight='bold')
        plt.xlabel("X 坐标 (km)")
        plt.ylabel("Y 坐标 (km)")
        plt.grid(True, linestyle='--', alpha=0.4)
        plt.legend(loc='upper left', fontsize=10)
        plt.show()

if __name__ == "__main__":
    start_time = time.time()
    optimizer = UrbanWasteLRPOptimizer(seed=123)
    
    # 阶段 1：设施选址分配
    station_assignments = optimizer.select_transfer_stations()
    
    if station_assignments:
        print("\n[*] 正在运行第二阶段：非对称路网下的多站点 VRP (Asymmetric CVRP)...")
        all_final_routes = {}
        total_lrp_distance = 0
        
        # 阶段 2：并行计算每个激活站点的子路线
        for station_id, assigned_pts in station_assignments.items():
            if not assigned_pts: continue
            
            routes = optimizer.solve_asymmetric_routing_for_station(
                station_id, assigned_pts, vehicle_capacity=8
            )
            all_final_routes[station_id] = routes
            
            # 核算真实非对称距离
            for r in routes:
                r_dist = sum(optimizer.asymmetric_dist_matrix[(r[i], r[i+1])] for i in range(len(r)-1))
                total_lrp_distance += r_dist
                print(f"    [HUB-{station_id}] 路线: {r} (单向实测距离: {r_dist:.2f} km)")
                
        print(f"\n[$$$] LRP 网络总实际行驶距离: {total_lrp_distance:.2f} km")
        
        # 阶段 3：高级有向图可视化
        optimizer.plot_lrp_directed_network(all_final_routes, list(station_assignments.keys()))

    print(f"\n[+] Total Execution Time: {time.time() - start_time:.2f} seconds")