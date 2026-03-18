"""
Urban Waste Collection VRP Optimizer - Question 1 (Complete Edition)
Author: Yuxiang Fan
Description: Solves the Capacitated Vehicle Routing Problem (CVRP) using:
1. Clarke-Wright Savings Algorithm
2. 2-opt Heuristic Local Search
3. OR-Tools Routing Model with Dynamic Vehicle Sweeping
4. Single-Vehicle TSP Baseline (Theoretical Lower Bound)
Includes built-in routing evaluation and matplotlib EDA tools with Chinese font support.
"""

import math
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from itertools import combinations, permutations
from ortools.sat.python import cp_model
from ortools.constraint_solver import pywrapcp, routing_enums_pb2

# 配置中文字体，支持图表中文字符显示 (衍生自 原始作图.py)
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False

class UrbanWasteOptimizer:
    def __init__(self, capacity=5, metric='euclidean'):
        self.capacity = capacity
        self.metric = metric
        self.depot = 0
        self.scale_waste = 10     
        self.scale_dist = 1000    
        
        # 统一的数据定义
        self.points = [
            {"id": 0, "x": 0, "y": 0, "waste": 0},
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
        
        self.nodes = [p["id"] for p in self.points]
        self.coords = {p["id"]: (p["x"], p["y"]) for p in self.points}
        self.waste = {p["id"]: p["waste"] for p in self.points}
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

    # -------------------- 评估与基础工具 --------------------
    def evaluate_route(self, route):
        total_dist = 0.0
        total_load = 0.0
        for i in range(len(route) - 1):
            total_dist += self.dist_matrix[(route[i], route[i+1])]
            if route[i+1] != self.depot:
                total_load += self.waste[route[i+1]]
        return total_dist, total_load

    # -------------------- 启发式算法 --------------------
    def solve_clarke_wright(self):
        savings = []
        for i, j in combinations(set(self.nodes) - {self.depot}, 2):
            s = self.dist_matrix[(self.depot, i)] + self.dist_matrix[(self.depot, j)] - self.dist_matrix[(i, j)]
            savings.append(((i, j), s))
        savings.sort(key=lambda x: x[1], reverse=True)

        tours = {i: [self.depot, i, self.depot] for i in self.nodes if i != self.depot}
        loads = {i: self.waste[i] for i in self.nodes if i != self.depot}
        route_map = {i: i for i in self.nodes if i != self.depot}

        for (i, j), _ in savings:
            r1, r2 = route_map[i], route_map[j]
            if r1 != r2 and loads[r1] + loads[r2] <= self.capacity:
                tours[r1].pop()
                tours[r1] += tours[r2][1:]
                loads[r1] += loads[r2]
                for node in tours[r2][1:-1]:
                    route_map[node] = r1
                del tours[r2]
                del loads[r2]
        return list(tours.values())

    def solve_2_opt(self, routes):
        optimized_routes = []
        for route in routes:
            best = route
            improved = True
            while improved:
                improved = False
                for i in range(1, len(best) - 2):
                    for j in range(i + 1, len(best) - 1):
                        new_route = best[:i] + best[i:j + 1][::-1] + best[j + 1:]
                        if self.evaluate_route(new_route)[0] < self.evaluate_route(best)[0]:
                            best = new_route
                            improved = True
            optimized_routes.append(best)
        return optimized_routes

    # -------------------- OR-Tools 高级求解器 --------------------
    def solve_tsp_baseline(self):
        """
        求解单车遍历所有节点的旅行商问题 (TSP) 
        作为项目的极限基线 (衍生自 一条.py)
        """
        print("\n[*] Calculating Single-Vehicle TSP Baseline (Infinite Capacity)...")
        data = {
            'distance_matrix': [[int(self.dist_matrix[(i, j)] * self.scale_dist) for j in self.nodes] for i in self.nodes],
            'num_vehicles': 1,
            'depot': self.depot
        }

        manager = pywrapcp.RoutingIndexManager(len(data['distance_matrix']), data['num_vehicles'], data['depot'])
        routing = pywrapcp.RoutingModel(manager)

        def distance_callback(from_index, to_index):
            return data['distance_matrix'][manager.IndexToNode(from_index)][manager.IndexToNode(to_index)]
        transit_callback_index = routing.RegisterTransitCallback(distance_callback)
        routing.SetArcCostEvaluatorOfAllVehicles(transit_callback_index)

        search_parameters = pywrapcp.DefaultRoutingSearchParameters()
        search_parameters.first_solution_strategy = routing_enums_pb2.FirstSolutionStrategy.PATH_CHEAPEST_ARC
        search_parameters.local_search_metaheuristic = routing_enums_pb2.LocalSearchMetaheuristic.GUIDED_LOCAL_SEARCH
        search_parameters.time_limit.FromSeconds(5)

        solution = routing.SolveWithParameters(search_parameters)
        if solution:
            index = routing.Start(0)
            route = []
            while not routing.IsEnd(index):
                route.append(manager.IndexToNode(index))
                index = solution.Value(routing.NextVar(index))
            route.append(self.depot)
            dist = solution.ObjectiveValue() / self.scale_dist
            print(f"[+] TSP Baseline Distance: {dist:.2f} km")
            return route, dist
        return None, None

    def sweep_ortools_routing(self, min_vehicles=2, max_vehicles=20, time_limit_per_sweep=5):
        best_overall_dist = float('inf')
        best_overall_routes = None
        best_v_count = 0

        print(f"\n[*] Starting Dynamic Vehicle Sweep ({min_vehicles} to {max_vehicles} vehicles)...")

        for v_count in range(min_vehicles, max_vehicles + 1):
            data = {
                'distance_matrix': [[int(self.dist_matrix[(i, j)] * self.scale_dist) for j in self.nodes] for i in self.nodes],
                'demands': [int(self.waste[i] * self.scale_dist) for i in self.nodes],
                'vehicle_capacities': [int(self.capacity * self.scale_dist)] * v_count,
                'num_vehicles': v_count,
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

            search_parameters = pywrapcp.DefaultRoutingSearchParameters()
            search_parameters.first_solution_strategy = routing_enums_pb2.FirstSolutionStrategy.PATH_CHEAPEST_ARC
            search_parameters.time_limit.FromSeconds(time_limit_per_sweep)
            search_parameters.log_search = False 

            solution = routing.SolveWithParameters(search_parameters)

            if solution:
                current_dist = solution.ObjectiveValue() / self.scale_dist
                print(f"  -> {v_count} Vehicles: Distance = {current_dist:.2f} km")
                
                if current_dist < best_overall_dist:
                    best_overall_dist = current_dist
                    best_v_count = v_count
                    
                    routes = []
                    for vehicle_id in range(data['num_vehicles']):
                        index = routing.Start(vehicle_id)
                        route = []
                        while not routing.IsEnd(index):
                            node = manager.IndexToNode(index)
                            route.append(node)
                            index = solution.Value(routing.NextVar(index))
                        route.append(self.depot)
                        if len(route) > 2:
                            routes.append(route)
                    best_overall_routes = routes

        print(f"[+] Sweep Complete! Best Configuration: {best_v_count} Vehicles, Min Distance: {best_overall_dist:.2f} km")
        return best_overall_routes

    # -------------------- 可视化模块 --------------------
    def plot_single_route_detail(self, route, distance=None):
        """
        绘制单条路线的详细特写 (衍生自 原始作图.py)
        包含中文字体支持和线段标注
        """
        plt.figure(figsize=(10, 8))
        x_coords = [self.coords[node][0] for node in route]
        y_coords = [self.coords[node][1] for node in route]
        
        # 绘制线段与散点
        plt.plot(x_coords, y_coords, color='royalblue', linewidth=2, marker='o', markersize=8, markerfacecolor='red')
        plt.scatter(self.coords[self.depot][0], self.coords[self.depot][1], color='gold', marker='*', s=300, zorder=5, label='车库/起点')
        
        # 标注文本
        for node in set(route):
            plt.text(self.coords[node][0] + 0.5, self.coords[node][1] + 0.5, str(node), fontsize=10, fontweight='bold')
            
        title = "单条清运路线特写分析"
        if distance:
            title += f" (总距离: {distance:.2f} km)"
            
        plt.title(title, fontsize=15, fontweight='bold')
        plt.xlabel("X 坐标 (km)", fontsize=12)
        plt.ylabel("Y 坐标 (km)", fontsize=12)
        plt.grid(True, linestyle='--', alpha=0.6)
        plt.legend(loc='upper right', fontsize=12)
        plt.show()

if __name__ == "__main__":
    start_time = time.time()
    optimizer = UrbanWasteOptimizer(capacity=5, metric='euclidean')
    
    # 1. 求解单车无容量限制的基线 (对应 一条.py)
    tsp_route, tsp_dist = optimizer.solve_tsp_baseline()
    
    # 2. 动态车辆数全局寻优
    best_routes = optimizer.sweep_ortools_routing(min_vehicles=5, max_vehicles=10, time_limit_per_sweep=3)
    
    # 3. 针对其中最长的一条路线进行详细特写绘图 (对应 原始作图.py)
    if best_routes:
        longest_route = max(best_routes, key=lambda r: optimizer.evaluate_route(r)[0])
        longest_dist, _ = optimizer.evaluate_route(longest_route)
        print("\n[*] Plotting detailed single route analysis...")
        optimizer.plot_single_route_detail(longest_route, distance=longest_dist)
        
    print(f"\n[+] Total Execution Time: {time.time() - start_time:.2f} seconds")