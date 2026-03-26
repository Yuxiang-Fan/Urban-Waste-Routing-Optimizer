"""
城市固体废弃物选址-路径问题优化系统
本模块实现两阶段优化策略：
第一阶段：利用 CP-SAT 求解带容量约束的设施选址模型，确定中转站开启状态及节点分配。
第二阶段：利用 OR-Tools 路径模型在非对称有向路网上，为各中转站规划清运路线。
"""

import math
import time
import random
import matplotlib.pyplot as plt
from ortools.sat.python import cp_model
from ortools.constraint_solver import pywrapcp, routing_enums_pb2

# 配置绘图字体
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False

class WasteLRPOptimizer:
    def __init__(self, seed=42):
        """
        初始化 LRP 优化器
        :param seed: 随机种子，用于复现非对称路网环境
        """
        self.scale_load = 1000     
        self.scale_dist = 1000    
        random.seed(seed)
        
        # 定义垃圾收集点数据
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
        
        # 定义候选垃圾中转站数据
        self.candidate_stations = {
            101: {"x": 50, "y": 50, "capacity": 25, "cost": 500},
            102: {"x": 60, "y": 60, "capacity": 30, "cost": 600},
            103: {"x": 70, "y": 70, "capacity": 20, "cost": 450},
            104: {"x": 80, "y": 80, "capacity": 35, "cost": 700},
            105: {"x": 90, "y": 90, "capacity": 15, "cost": 400}
        }
        
        # 整合所有节点信息
        self.all_nodes = {p["id"]: {"x": p["x"], "y": p["y"], "waste": p["waste"]} for p in self.collection_points}
        for ts_id, ts_data in self.candidate_stations.items():
            self.all_nodes[ts_id] = {"x": ts_data["x"], "y": ts_data["y"], "waste": 0}
            
        self.asymmetric_dist_matrix = self._build_asymmetric_matrix()

    def _build_asymmetric_matrix(self):
        """
        构建非对称距离矩阵。
        模拟城市环境中的单行道或单向交通拥堵，即 D(i,j) 与 D(j,i) 不相等。
        """
        dist_matrix = {}
        nodes = list(self.all_nodes.keys())
        for i in nodes:
            for j in nodes:
                if i == j:
                    dist_matrix[(i, j)] = 0.0
                else:
                    # 计算基础欧氏距离
                    base_dist = math.hypot(self.all_nodes[i]["x"] - self.all_nodes[j]["x"], 
                                           self.all_nodes[i]["y"] - self.all_nodes[j]["y"])
                    # 模拟部分路段存在非对称特征，增加 30% 到 80% 的行驶成本
                    factor = 1.0
                    if random.random() < 0.1:
                        factor = random.uniform(1.3, 1.8)
                    dist_matrix[(i, j)] = base_dist * factor
        return dist_matrix

    def run_facility_location(self):
        """
        第一阶段：基于 CP-SAT 的设施选址优化。
        目标是在满足中转站载荷能力的前提下，最小化建设成本与分配距离成本之和。
        """
        print("执行第一阶段：中转站选址与节点分配...")
        model = cp_model.CpModel()
        
        points = [p["id"] for p in self.collection_points]
        stations = list(self.candidate_stations.keys())
        
        # 决策变量：y 为布尔变量，表示是否启用该中转站
        y = {s: model.NewBoolVar(f"y_{s}") for s in stations}
        # 决策变量：x 为分配变量，表示收集点是否由该中转站负责
        x = {(p, s): model.NewBoolVar(f"x_{p}_{s}") for p in points for s in stations}
        
        # 每个收集点必须分配给且仅分配给一个已启用的中转站
        for p in points:
            model.AddExactlyOne([x[(p, s)] for s in stations])
            for s in stations:
                model.AddImplication(x[(p, s)], y[s])
                
        # 中转站容量约束：分配的垃圾总量不能超过该站额定容量
        for s in stations:
            cap_int = int(self.candidate_stations[s]["capacity"] * self.scale_load)
            point_loads = [int(self.all_nodes[p]["waste"] * self.scale_load) * x[(p, s)] for p in points]
            model.Add(sum(point_loads) <= cap_int * y[s])
            
        # 目标函数：最小化总运输距离 + 设施建设固定成本
        transport_dist = sum(int(self.asymmetric_dist_matrix[(p, s)] * self.scale_dist) * x[(p, s)] 
                             for p in points for s in stations)
        fixed_cost = sum(int(self.candidate_stations[s]["cost"] * self.scale_dist) * y[s] for s in stations)
        model.Minimize(transport_dist + fixed_cost)
        
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
            
            print(f"选址决策完成。启用站点: {active_stations}")
            return assignments
        return None

    def solve_routing_for_station(self, station_id, assigned_points, vehicle_cap=8):
        """
        第二阶段：为各中转站所辖节点进行非对称 CVRP 路径规划。
        """
        # 建立局部索引映射
        nodes_list = [station_id] + assigned_points
        n_count = len(nodes_list)
        
        # 提取非对称距离矩阵并放大
        dist_mtx = [[int(self.asymmetric_dist_matrix[(i, j)] * self.scale_dist) for j in nodes_list] for i in nodes_list]
        demands = [0] + [int(self.all_nodes[p]["waste"] * self.scale_load) for p in assigned_points]
        
        # 初始化管理器与模型，起始点设定为中转站
        manager = pywrapcp.RoutingIndexManager(n_count, len(assigned_points), 0)
        routing = pywrapcp.RoutingModel(manager)

        def distance_callback(from_idx, to_idx):
            return dist_mtx[manager.IndexToNode(from_idx)][manager.IndexToNode(to_idx)]
        
        transit_callback_idx = routing.RegisterTransitCallback(distance_callback)
        routing.SetArcCostEvaluatorOfAllVehicles(transit_callback_idx)

        def demand_callback(from_idx):
            return demands[manager.IndexToNode(from_idx)]
        
        demand_callback_idx = routing.RegisterUnaryTransitCallback(demand_callback)
        routing.AddDimensionWithVehicleCapacity(
            demand_callback_idx, 0, [int(vehicle_cap * self.scale_load)] * len(assigned_points), True, 'Capacity'
        )

        search_params = pywrapcp.DefaultRoutingSearchParameters()
        search_params.first_solution_strategy = routing_enums_pb2.FirstSolutionStrategy.PATH_CHEAPEST_ARC
        search_params.local_search_metaheuristic = routing_enums_pb2.LocalSearchMetaheuristic.GUIDED_LOCAL_SEARCH
        search_params.time_limit.FromSeconds(3)

        solution = routing.SolveWithParameters(search_params)

        if solution:
            routes = []
            for v_id in range(len(assigned_points)):
                idx = routing.Start(v_id)
                path = []
                while not routing.IsEnd(idx):
                    path.append(nodes_list[manager.IndexToNode(idx)])
                    idx = solution.Value(routing.NextVar(idx))
                path.append(nodes_list[manager.IndexToNode(idx)])
                if len(path) > 2:
                    routes.append(path)
            return routes
        return []

    def plot_asymmetric_lrp_network(self, routes_dict, active_hubs):
        """
        绘制 LRP 优化结果。
        通过带箭头的有向边展示非对称路网的流向特征。
        """
        plt.figure(figsize=(12, 10))
        
        # 绘制收集点
        px = [p["x"] for p in self.collection_points]
        py = [p["y"] for p in self.collection_points]
        plt.scatter(px, py, color='cornflowerblue', s=30, zorder=3, label='收集点')

        # 绘制中转站
        for s_id, s_data in self.candidate_stations.items():
            if s_id in active_hubs:
                plt.scatter(s_data["x"], s_data["y"], color='crimson', marker='^', s=200, zorder=5, label='已启用中转站' if s_id == active_hubs[0] else "")
                plt.text(s_data["x"]+1, s_data["y"]+1, f"中转站-{s_id}", color='crimson', fontweight='bold')
            else:
                plt.scatter(s_data["x"], s_data["y"], color='lightgrey', marker='^', s=100, zorder=2)

        # 绘制有向清运路径
        color_cycle = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
        c_idx = 0
        
        for hub_id, routes in routes_dict.items():
            for r in routes:
                color = color_cycle[c_idx % len(color_cycle)]
                for i in range(len(r) - 1):
                    start_node, end_node = r[i], r[i+1]
                    x1, y1 = self.all_nodes[start_node]["x"], self.all_nodes[start_node]["y"]
                    x2, y2 = self.all_nodes[end_node]["x"], self.all_nodes[end_node]["y"]
                    # 利用 annotate 函数绘制带箭头的实线，体现有向性
                    plt.annotate('', xy=(x2, y2), xytext=(x1, y1),
                                 arrowprops=dict(arrowstyle="->", color=color, lw=1.2, alpha=0.6))
                c_idx += 1

        plt.title("选址-路径问题 (LRP) 优化方案可视化 - 非对称有向路网", fontsize=15)
        plt.grid(True, linestyle=':', alpha=0.5)
        plt.legend(loc='upper left')
        plt.show()

if __name__ == "__main__":
    start_timer = time.time()
    lrp_solver = WasteLRPOptimizer(seed=123)
    
    # 第一阶段：设施选址
    station_map = lrp_solver.run_facility_location()
    
    if station_map:
        print("执行第二阶段：基于局部子集的非对称路径规划...")
        final_results = {}
        total_dist_sum = 0
        
        for hub, pts in station_map.items():
            if not pts: continue
            
            paths = lrp_solver.solve_routing_for_station(hub, pts)
            final_results[hub] = paths
            
            # 统计实际行驶距离（基于非对称矩阵）
            for p in paths:
                d = sum(lrp_solver.asymmetric_dist_matrix[(p[i], p[i+1])] for i in range(len(p)-1))
                total_dist_sum += d
                print(f"中转站 {hub} 规划路径: {p}，实测距离: {d:.2f} km")
                
        print(f"\n系统总清运里程: {total_dist_sum:.2f} km")
        lrp_solver.plot_asymmetric_lrp_network(final_results, list(station_map.keys()))

    print(f"任务总执行耗时: {time.time() - start_timer:.2f} 秒")
