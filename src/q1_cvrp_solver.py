"""
基于启发式算法与 OR-Tools 的城市垃圾清运路径优化 (CVRP)
包含 CW 节约算法、2-opt 局部搜索以及 OR-Tools 动态车辆数寻优
"""

import math
import time
import matplotlib.pyplot as plt
from itertools import combinations
from ortools.constraint_solver import pywrapcp, routing_enums_pb2

# 配置 matplotlib 中文字体，避免图表乱码
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False

class UrbanWasteOptimizer:
    def __init__(self, capacity=5, metric='euclidean'):
        """
        初始化 CVRP 优化器
        :param capacity: 单辆垃圾车的最大载重量
        :param metric: 距离度量方式 ('euclidean' 或 'manhattan')
        """
        self.capacity = capacity
        self.metric = metric
        self.depot = 0
        
        # OR-Tools 只能处理整数，这里设置距离和载重的放大倍数以保留精度
        self.scale_dist = 1000    
        self.scale_waste = 1000     
        
        # 节点数据：id为0的是车库，其余为垃圾收集点，waste为当前点的垃圾量
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
        
        # 预处理数据结构，方便后续快速索引
        self.nodes = [p["id"] for p in self.points]
        self.coords = {p["id"]: (p["x"], p["y"]) for p in self.points}
        self.waste = {p["id"]: p["waste"] for p in self.points}
        
        # 初始化距离矩阵
        self.dist_matrix = self._build_distance_matrix()

    def _calculate_distance(self, node1, node2):
        """计算两点间的物理距离"""
        x1, y1 = self.coords[node1]
        x2, y2 = self.coords[node2]
        if self.metric == 'euclidean':
            return math.hypot(x1 - x2, y1 - y2)
        elif self.metric == 'manhattan':
            return abs(x2 - x1) + abs(y2 - y1)
        else:
            raise ValueError("Unsupported distance metric")

    def _build_distance_matrix(self):
        """构建全连接的距离矩阵"""
        return {(i, j): self._calculate_distance(i, j) for i in self.nodes for j in self.nodes}

    def evaluate_route(self, route):
        """
        计算单条路径的行驶距离和总装载量
        :param route: 节点列表，例如 [0, 1, 2, 0]
        """
        total_dist = 0.0
        total_load = 0.0
        for i in range(len(route) - 1):
            total_dist += self.dist_matrix[(route[i], route[i+1])]
            current_node = route[i+1]
            if current_node != self.depot:
                total_load += self.waste[current_node]
        return total_dist, total_load

    def solve_clarke_wright(self):
        """使用 Clarke-Wright 节约算法生成初始解"""
        savings = []
        # 计算所有非起点节点对的节约值：S(i,j) = D(depot,i) + D(depot,j) - D(i,j)
        for i, j in combinations(set(self.nodes) - {self.depot}, 2):
            s = self.dist_matrix[(self.depot, i)] + self.dist_matrix[(self.depot, j)] - self.dist_matrix[(i, j)]
            savings.append(((i, j), s))
        
        # 按节约值降序排列
        savings.sort(key=lambda x: x[1], reverse=True)

        # 初始化：每个节点单独分配一辆车
        tours = {i: [self.depot, i, self.depot] for i in self.nodes if i != self.depot}
        loads = {i: self.waste[i] for i in self.nodes if i != self.depot}
        route_map = {i: i for i in self.nodes if i != self.depot}

        # 贪心合并路径
        for (i, j), _ in savings:
            r1, r2 = route_map[i], route_map[j]
            # 如果不在同一条路径，且合并后不超过车辆容量
            if r1 != r2 and loads[r1] + loads[r2] <= self.capacity:
                tours[r1].pop()  # 移除 r1 结尾的 depot
                tours[r1] += tours[r2][1:]  # 接入 r2 (跳过 r2 开头的 depot)
                loads[r1] += loads[r2]
                
                # 更新路由映射表
                for node in tours[r2][1:-1]:
                    route_map[node] = r1
                    
                del tours[r2]
                del loads[r2]
                
        return list(tours.values())

    def solve_2_opt(self, routes):
        """对给定的路径集合执行 2-opt 局部搜索优化"""
        optimized_routes = []
        for route in routes:
            best = route
            improved = True
            while improved:
                improved = False
                # 遍历所有可能的边翻转组合
                for i in range(1, len(best) - 2):
                    for j in range(i + 1, len(best) - 1):
                        # 翻转 i 到 j 之间的路径
                        new_route = best[:i] + best[i:j + 1][::-1] + best[j + 1:]
                        if self.evaluate_route(new_route)[0] < self.evaluate_route(best)[0]:
                            best = new_route
                            improved = True
            optimized_routes.append(best)
        return optimized_routes

    def solve_tsp_baseline(self):
        """
        计算单车遍历所有节点（不考虑容量限制）的 TSP 模型。
        该结果作为多车 CVRP 的理论距离下界参考。
        """
        print("\n--- 正在计算单车 TSP 极限基线 (无容量限制) ---")
        
        # 数据转换：将浮点数按比例放大为整数
        data = {
            'distance_matrix': [[int(self.dist_matrix[(i, j)] * self.scale_dist) for j in self.nodes] for i in self.nodes],
            'num_vehicles': 1,
            'depot': self.depot
        }

        manager = pywrapcp.RoutingIndexManager(len(data['distance_matrix']), data['num_vehicles'], data['depot'])
        routing = pywrapcp.RoutingModel(manager)

        def distance_callback(from_index, to_index):
            from_node = manager.IndexToNode(from_index)
            to_node = manager.IndexToNode(to_index)
            return data['distance_matrix'][from_node][to_node]
            
        transit_callback_index = routing.RegisterTransitCallback(distance_callback)
        routing.SetArcCostEvaluatorOfAllVehicles(transit_callback_index)

        # 设置启发式搜索策略：GUIDED_LOCAL_SEARCH 往往能得到更好的全局解
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
            
            # 将距离还原回公里单位
            dist = solution.ObjectiveValue() / self.scale_dist
            print(f"TSP 基准距离: {dist:.2f} km")
            return route, dist
        return None, None

    def sweep_ortools_routing(self, min_vehicles=2, max_vehicles=20, time_limit_per_sweep=5):
        """
        动态遍历可用车辆数，利用 OR-Tools 寻找整体距离最短的车队配置。
        """
        best_overall_dist = float('inf')
        best_overall_routes = None
        best_v_count = 0

        print(f"\n--- 启动动态车辆数寻优 ({min_vehicles} 到 {max_vehicles} 辆车) ---")

        for v_count in range(min_vehicles, max_vehicles + 1):
            data = {
                'distance_matrix': [[int(self.dist_matrix[(i, j)] * self.scale_dist) for j in self.nodes] for i in self.nodes],
                'demands': [int(self.waste[i] * self.scale_waste) for i in self.nodes],
                'vehicle_capacities': [int(self.capacity * self.scale_waste)] * v_count,
                'num_vehicles': v_count,
                'depot': self.depot
            }

            manager = pywrapcp.RoutingIndexManager(len(data['distance_matrix']), data['num_vehicles'], data['depot'])
            routing = pywrapcp.RoutingModel(manager)

            # 距离回调
            def distance_callback(from_index, to_index):
                from_node = manager.IndexToNode(from_index)
                to_node = manager.IndexToNode(to_index)
                return data['distance_matrix'][from_node][to_node]
                
            transit_callback_index = routing.RegisterTransitCallback(distance_callback)
            routing.SetArcCostEvaluatorOfAllVehicles(transit_callback_index)

            # 容量约束回调
            def demand_callback(from_index):
                from_node = manager.IndexToNode(from_index)
                return data['demands'][from_node]
                
            demand_callback_index = routing.RegisterUnaryTransitCallback(demand_callback)
            routing.AddDimensionWithVehicleCapacity(
                demand_callback_index, 
                0,  # 无容量松弛
                data['vehicle_capacities'], 
                True,  # 累计需求从 0 开始
                'Capacity'
            )

            search_parameters = pywrapcp.DefaultRoutingSearchParameters()
            search_parameters.first_solution_strategy = routing_enums_pb2.FirstSolutionStrategy.PATH_CHEAPEST_ARC
            search_parameters.time_limit.FromSeconds(time_limit_per_sweep)
            search_parameters.log_search = False 

            solution = routing.SolveWithParameters(search_parameters)

            if solution:
                current_dist = solution.ObjectiveValue() / self.scale_dist
                print(f"尝试分配 {v_count} 辆车: 总行驶距离 = {current_dist:.2f} km")
                
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
                        
                        # 过滤掉未分配任务的空车 (仅包含起点和终点)
                        if len(route) > 2:
                            routes.append(route)
                    best_overall_routes = routes

        print(f"\n寻优完成！最优配置: {best_v_count} 辆车, 最小总距离: {best_overall_dist:.2f} km")
        return best_overall_routes

    def plot_single_route_detail(self, route, distance=None):
        """
        绘制指定路线的二维坐标图
        """
        plt.figure(figsize=(10, 8))
        x_coords = [self.coords[node][0] for node in route]
        y_coords = [self.coords[node][1] for node in route]
        
        # 绘制路线轨迹
        plt.plot(x_coords, y_coords, color='royalblue', linewidth=2, marker='o', markersize=8, markerfacecolor='red')
        
        # 突出显示车库
        plt.scatter(self.coords[self.depot][0], self.coords[self.depot][1], color='gold', marker='*', s=300, zorder=5, label='车库/起点')
        
        # 为途径的节点添加编号标注
        for node in set(route):
            plt.text(self.coords[node][0] + 0.5, self.coords[node][1] + 0.5, str(node), fontsize=10, fontweight='bold')
            
        title = "单条清运路线特写分析"
        if distance:
            title += f" (路线距离: {distance:.2f} km)"
            
        plt.title(title, fontsize=15, fontweight='bold')
        plt.xlabel("X 坐标 (km)", fontsize=12)
        plt.ylabel("Y 坐标 (km)", fontsize=12)
        plt.grid(True, linestyle='--', alpha=0.6)
        plt.legend(loc='upper right', fontsize=12)
        
        # 关闭阻塞，确保在某些IDE中正常运行
        plt.show(block=True)

if __name__ == "__main__":
    start_time = time.time()
    
    # 实例化优化器
    optimizer = UrbanWasteOptimizer(capacity=5, metric='euclidean')
    
    # 1. 计算单车 TSP 作为理想下界
    tsp_route, tsp_dist = optimizer.solve_tsp_baseline()
    
    # 2. 动态车辆数全局寻优 (测试 5-10 辆车的情况)
    best_routes = optimizer.sweep_ortools_routing(min_vehicles=5, max_vehicles=10, time_limit_per_sweep=3)
    
    # 3. 对规划出距离最长的一条路线进行绘图分析
    if best_routes:
        longest_route = max(best_routes, key=lambda r: optimizer.evaluate_route(r)[0])
        longest_dist, _ = optimizer.evaluate_route(longest_route)
        print("\n--- 正在生成最长路线的可视化图表 ---")
        optimizer.plot_single_route_detail(longest_route, distance=longest_dist)
        
    print(f"\n总运行时间: {time.time() - start_time:.2f} 秒")
