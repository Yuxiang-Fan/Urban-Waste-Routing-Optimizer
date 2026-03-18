"""
Urban Waste Collection VRP Optimizer - Question 2 (Ultimate Multi-Fleet Edition)
Author: Y. Fan
Description: Solves the multi-commodity CVRP for different waste types.
Features a hybrid solving strategy:
1. OR-Tools Routing Model for high-volume waste (Food, Recyclable, Other) with capacity/time constraints.
2. Exact CP-SAT Model with MTZ formulation for low-volume waste (Hazardous) pure TSP.
"""

import math
import time
import matplotlib.pyplot as plt
from ortools.constraint_solver import pywrapcp, routing_enums_pb2
from ortools.sat.python import cp_model

plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False

class MultiFleetWasteOptimizer:
    def __init__(self, metric='euclidean'):
        self.metric = metric
        self.depot = 0
        self.scale_waste = 10     
        self.scale_dist = 1000    
        
        self.points = [
            {"id": 0, "x": 0, "y": 0}, {"id": 1, "x": 12, "y": 8}, {"id": 2, "x": 5, "y": 15},
            {"id": 3, "x": 20, "y": 30}, {"id": 4, "x": 25, "y": 10}, {"id": 5, "x": 35, "y": 22},
            {"id": 6, "x": 18, "y": 5}, {"id": 7, "x": 30, "y": 35}, {"id": 8, "x": 10, "y": 25},
            {"id": 9, "x": 22, "y": 18}, {"id": 10, "x": 38, "y": 15}, {"id": 11, "x": 5, "y": 8},
            {"id": 12, "x": 15, "y": 32}, {"id": 13, "x": 28, "y": 5}, {"id": 14, "x": 30, "y": 12},
            {"id": 15, "x": 10, "y": 10}, {"id": 16, "x": 20, "y": 20}, {"id": 17, "x": 35, "y": 30},
            {"id": 18, "x": 8, "y": 22}, {"id": 19, "x": 25, "y": 25}, {"id": 20, "x": 32, "y": 8},
            {"id": 21, "x": 15, "y": 5}, {"id": 22, "x": 28, "y": 20}, {"id": 23, "x": 38, "y": 25},
            {"id": 24, "x": 10, "y": 30}, {"id": 25, "x": 20, "y": 10}, {"id": 26, "x": 30, "y": 18},
            {"id": 27, "x": 5, "y": 25}, {"id": 28, "x": 18, "y": 30}, {"id": 29, "x": 35, "y": 10},
            {"id": 30, "x": 22, "y": 35}
        ]
        
        self.nodes = [p["id"] for p in self.points]
        self.coords = {p["id"]: (p["x"], p["y"]) for p in self.points}
        self.dist_matrix = self._build_distance_matrix()

    def _calculate_distance(self, node1, node2):
        x1, y1 = self.coords[node1]
        x2, y2 = self.coords[node2]
        if self.metric == 'euclidean':
            return math.hypot(x1 - x2, y1 - y2)
        elif self.metric == 'manhattan':
            return abs(x2 - x1) + abs(y2 - y1)

    def _build_distance_matrix(self):
        return {(i, j): self._calculate_distance(i, j) for i in self.nodes for j in self.nodes}

    # ------------------ 大运量垃圾 VRP 求解器 ------------------
    def solve_specific_waste_type(self, waste_name, demands, capacity, max_distance, max_vehicles=20, time_limit=5):
        print(f"\n[*] 正在调度 [{waste_name}] 收运车队 (Routing Model)...")
        data = {
            'distance_matrix': [[int(self.dist_matrix[(i, j)] * self.scale_dist) for j in self.nodes] for i in self.nodes],
            'demands': [int(demands.get(i, 0) * self.scale_dist) for i in self.nodes],
            'vehicle_capacities': [int(capacity * self.scale_dist)] * max_vehicles,
            'num_vehicles': max_vehicles,
            'depot': self.depot
        }

        manager = pywrapcp.RoutingIndexManager(len(data['distance_matrix']), data['num_vehicles'], data['depot'])
        routing = pywrapcp.RoutingModel(manager)

        def distance_callback(from_index, to_index):
            return data['distance_matrix'][manager.IndexToNode(from_index)][manager.IndexToNode(to_index)]
        transit_callback_index = routing.RegisterTransitCallback(distance_callback)
        routing.SetArcCostEvaluatorOfAllVehicles(transit_callback_index)

        def demand_callback(from_index):
            return data['demands'][manager.IndexToNode(from_index)]
        demand_callback_index = routing.RegisterUnaryTransitCallback(demand_callback)
        routing.AddDimensionWithVehicleCapacity(demand_callback_index, 0, data['vehicle_capacities'], True, 'Capacity')

        routing.AddDimension(transit_callback_index, 0, int(max_distance * self.scale_dist), True, 'Distance')

        search_parameters = pywrapcp.DefaultRoutingSearchParameters()
        search_parameters.first_solution_strategy = routing_enums_pb2.FirstSolutionStrategy.PATH_CHEAPEST_ARC
        search_parameters.local_search_metaheuristic = routing_enums_pb2.LocalSearchMetaheuristic.GUIDED_LOCAL_SEARCH
        search_parameters.time_limit.FromSeconds(time_limit)
        search_parameters.log_search = False

        solution = routing.SolveWithParameters(search_parameters)

        if solution:
            routes = []
            total_dist = 0
            for vehicle_id in range(data['num_vehicles']):
                index = routing.Start(vehicle_id)
                route = []
                route_dist = 0
                while not routing.IsEnd(index):
                    node = manager.IndexToNode(index)
                    route.append(node)
                    previous_index = index
                    index = solution.Value(routing.NextVar(index))
                    route_dist += routing.GetArcCostForVehicle(previous_index, index, vehicle_id)
                route.append(self.depot)
                
                if len(route) > 2:
                    routes.append(route)
                    total_dist += route_dist
                    print(f"    -> 车辆 {len(routes)} 路线: {route} (距离: {route_dist / self.scale_dist:.2f} km)")
            return routes, total_dist / self.scale_dist
        return None, 0

    # ------------------ 小运量垃圾 TSP 求解器  ------------------
    def solve_hazardous_tsp_cpsat(self, waste_name="有害垃圾", time_limit=10):
        print(f"\n[*] 正在调度 [{waste_name}] 专车 (CP-SAT 精确 TSP 模型)...")
        n = len(self.nodes)
        model = cp_model.CpModel()
        
        # 决策变量
        x = {(i, j): model.NewBoolVar(f"x_{i}_{j}") for i in range(n) for j in range(n) if i != j}
        # MTZ 消除子回路变量
        u = {i: model.NewIntVar(1, n-1, f"u_{i}") for i in range(1, n)}

        # 节点出入度约束 (每个点进一次，出一次)
        for i in range(n):
            model.Add(sum(x[(i, j)] for j in range(n) if i != j) == 1)
            model.Add(sum(x[(j, i)] for j in range(n) if i != j) == 1)
            
        # MTZ 约束 (防止未包含起点的局部小圈)
        for i in range(1, n):
            for j in range(1, n):
                if i != j:
                    model.Add(u[i] - u[j] + n * x[(i, j)] <= n - 1)

        # 目标函数：最小化总距离
        dist_int = {(i, j): int(self.dist_matrix[(i, j)] * self.scale_dist) for i in range(n) for j in range(n)}
        model.Minimize(sum(dist_int[(i, j)] * x[(i, j)] for i in range(n) for j in range(n) if i != j))

        solver = cp_model.CpSolver()
        solver.parameters.max_time_in_seconds = time_limit
        status = solver.Solve(model)

        if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
            total_dist = solver.ObjectiveValue() / self.scale_dist
            # 提取路径
            route = [self.depot]
            current = self.depot
            visited = {self.depot}
            while len(route) < n:
                for j in range(n):
                    if current != j and solver.Value(x[(current, j)]) == 1:
                        route.append(j)
                        visited.add(j)
                        current = j
                        break
            route.append(self.depot)
            print(f"    -> 专车路线: {route} (距离: {total_dist:.2f} km)")
            return [route], total_dist
        else:
            print(f"[-] [{waste_name}] CP-SAT 求解失败。")
            return None, 0

    # ------------------ 可视化模块 ------------------
    def plot_multi_fleet_routes(self, fleet_routes_dict):
        plt.figure(figsize=(12, 10))
        for i in self.nodes:
            plt.scatter(self.coords[i][0], self.coords[i][1], color='gray', zorder=2)
            plt.text(self.coords[i][0] + 0.5, self.coords[i][1] + 0.5, str(i), fontsize=8)
        plt.scatter(self.coords[self.depot][0], self.coords[self.depot][1], color='red', marker='*', s=200, zorder=5, label='处理中心')

        color_map = {'厨余垃圾': 'green', '可回收垃圾': 'blue', '其他垃圾': 'orange', '有害垃圾': 'red'}
        
        for waste_name, routes in fleet_routes_dict.items():
            base_color = color_map.get(waste_name, 'black')
            for idx, route in enumerate(routes):
                x_coords = [self.coords[node][0] for node in route]
                y_coords = [self.coords[node][1] for node in route]
                # 有害垃圾使用虚线加粗突出显示
                line_style = '--' if waste_name == '有害垃圾' else '-'
                line_width = 3 if waste_name == '有害垃圾' else 2
                label = f"{waste_name}" if idx == 0 else None
                plt.plot(x_coords, y_coords, color=base_color, linestyle=line_style, linewidth=line_width, alpha=0.7, label=label)

        plt.title("多车型协同垃圾清运 (含 MTZ-TSP 精确解)", fontsize=16, fontweight='bold')
        plt.grid(True, linestyle='--', alpha=0.6)
        plt.legend(loc='upper right', fontsize=10)
        plt.show()

if __name__ == "__main__":
    start_time = time.time()
    optimizer = MultiFleetWasteOptimizer(metric='euclidean')
    
    # [模拟数据] 
    sample_demands = {i: 1.0 for i in range(1, 31)}
    fleet_routes = {}
    
    # 1. 厨余垃圾 (VRP)
    routes_food, _ = optimizer.solve_specific_waste_type("厨余垃圾", sample_demands, capacity=5, max_distance=80)
    if routes_food: fleet_routes["厨余垃圾"] = routes_food
        
    # 2. 其他垃圾 (VRP)
    routes_other, _ = optimizer.solve_specific_waste_type("其他垃圾", sample_demands, capacity=6, max_distance=90)
    if routes_other: fleet_routes["其他垃圾"] = routes_other

    # 3. 有害垃圾 (CP-SAT MTZ-TSP 精确求解)
    routes_hazard, _ = optimizer.solve_hazardous_tsp_cpsat("有害垃圾", time_limit=10)
    if routes_hazard: fleet_routes["有害垃圾"] = routes_hazard

    # 绘制结果
    optimizer.plot_multi_fleet_routes(fleet_routes)
    print(f"\n[+] Total Execution Time: {time.time() - start_time:.2f} seconds")
