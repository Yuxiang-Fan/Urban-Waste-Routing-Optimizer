"""
城市多类别垃圾协同收运优化系统
针对不同类型的垃圾（厨余、可回收、有害等）采用不同的优化策略：
1. 大运量垃圾：利用 OR-Tools 路径模型处理带容量和最大里程约束的 CVRP。
2. 有害垃圾：利用 CP-SAT 构建带有 MTZ 约束的精确 TSP 模型进行专项收运。
"""

import math
import time
import matplotlib.pyplot as plt
from ortools.constraint_solver import pywrapcp, routing_enums_pb2
from ortools.sat.python import cp_model

# 绘图环境配置
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False

class MultiFleetWasteOptimizer:
    def __init__(self, metric='euclidean'):
        """
        初始化多车型优化器
        :param metric: 距离度量标准
        """
        self.metric = metric
        self.depot = 0
        # 放大系数：将浮点数转化为整数以适配求解器精度要求
        self.scale_dist = 1000    
        self.scale_load = 1000    
        
        # 垃圾收集点坐标数据
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
        """计算节点间欧氏距离或曼哈顿距离"""
        x1, y1 = self.coords[node1]
        x2, y2 = self.coords[node2]
        if self.metric == 'euclidean':
            return math.hypot(x1 - x2, y1 - y2)
        elif self.metric == 'manhattan':
            return abs(x2 - x1) + abs(y2 - y1)

    def _build_distance_matrix(self):
        """预计算全节点的距离矩阵"""
        return {(i, j): self._calculate_distance(i, j) for i in self.nodes for j in self.nodes}

    # ------------------ 大运量垃圾：CVRP 模型 ------------------
    def solve_standard_vrp(self, waste_type, demands, capacity, max_dist, vehicle_num=10):
        """
        基于 OR-Tools Routing 框架求解带载重和里程限制的路径问题
        """
        print(f"正在优化 [{waste_type}] 清运路线...")
        
        # 数据整理与放大处理
        data = {
            'distance_matrix': [[int(self.dist_matrix[(i, j)] * self.scale_dist) for j in self.nodes] for i in self.nodes],
            'demands': [int(demands.get(i, 0) * self.scale_load) for i in self.nodes],
            'vehicle_capacities': [int(capacity * self.scale_load)] * vehicle_num,
            'num_vehicles': vehicle_num,
            'depot': self.depot
        }

        manager = pywrapcp.RoutingIndexManager(len(data['distance_matrix']), data['num_vehicles'], data['depot'])
        routing = pywrapcp.RoutingModel(manager)

        # 注册距离回调函数，用于计算路径总长度
        def dist_callback(from_index, to_index):
            return data['distance_matrix'][manager.IndexToNode(from_index)][manager.IndexToNode(to_index)]
        transit_idx = routing.RegisterTransitCallback(dist_callback)
        routing.SetArcCostEvaluatorOfAllVehicles(transit_idx)

        # 注册需求回调函数，用于处理载重约束
        def demand_callback(from_index):
            return data['demands'][manager.IndexToNode(from_index)]
        demand_idx = routing.RegisterUnaryTransitCallback(demand_callback)
        routing.AddDimensionWithVehicleCapacity(demand_idx, 0, data['vehicle_capacities'], True, 'Capacity')

        # 添加车辆最大行驶距离约束
        routing.AddDimension(transit_idx, 0, int(max_dist * self.scale_dist), True, 'Distance')

        # 设置搜索策略：采用启发式初始解 + 引导局部搜索
        search_params = pywrapcp.DefaultRoutingSearchParameters()
        search_params.first_solution_strategy = routing_enums_pb2.FirstSolutionStrategy.PATH_CHEAPEST_ARC
        search_params.local_search_metaheuristic = routing_enums_pb2.LocalSearchMetaheuristic.GUIDED_LOCAL_SEARCH
        search_params.time_limit.FromSeconds(5)

        solution = routing.SolveWithParameters(search_params)

        if solution:
            final_routes = []
            total_km = 0
            for v_id in range(data['num_vehicles']):
                idx = routing.Start(v_id)
                route = []
                while not routing.IsEnd(idx):
                    route.append(manager.IndexToNode(idx))
                    idx = solution.Value(routing.NextVar(idx))
                route.append(self.depot)
                
                # 仅记录实际分配了任务的车辆路径
                if len(route) > 2:
                    final_routes.append(route)
                    d_dim = routing.GetDimensionOrDie('Distance')
                    total_km += solution.Value(d_dim.CumulVar(idx))
            
            print(f"[{waste_type}] 求解完成，共使用车辆: {len(final_routes)}，总里程: {total_km / self.scale_dist:.2f} km")
            return final_routes
        return []

    # ------------------ 有害垃圾：精确 TSP 模型 ------------------
    def solve_hazardous_tsp(self):
        """
        针对有害垃圾，利用 CP-SAT 构建精确数学模型。
        使用 MTZ 约束消除子回路。
        """
        print("正在为 [有害垃圾] 构建精确 TSP 求解模型...")
        n = len(self.nodes)
        model = cp_model.CpModel()
        
        # 决策变量：x[i,j] 表示是否经过边 (i, j)
        x = {}
        for i in range(n):
            for j in range(n):
                if i != j:
                    x[i, j] = model.NewBoolVar(f'x_{i}_{j}')

        # 辅助变量：u[i] 用于 MTZ 约束中的序列排序
        u = [model.NewIntVar(0, n - 1, f'u_{i}') for i in range(n)]

        # 约束条件 1：每个点必须有且仅有一条入边和一条出边
        for i in range(n):
            model.Add(sum(x[i, j] for j in range(n) if i != j) == 1)
            model.Add(sum(x[j, i] for j in range(n) if i != j) == 1)
            
        # 约束条件 2：MTZ (Miller-Tucker-Zemlin) 子回路消除约束
        # 对于所有不包含仓库的边 (i, j)，满足 u[i] - u[j] + n*x[i,j] <= n-1
        for i in range(1, n):
            for j in range(1, n):
                if i != j:
                    model.Add(u[i] - u[j] + n * x[i, j] <= n - 1)

        # 目标函数：最小化总行驶距离
        dist_costs = []
        for (i, j), var in x.items():
            dist_costs.append(int(self.dist_matrix[(i, j)] * self.scale_dist) * var)
        model.Minimize(sum(dist_costs))

        solver = cp_model.CpSolver()
        solver.parameters.max_time_in_seconds = 10
        status = solver.Solve(model)

        if status == cp_model.OPTIMAL or status == cp_model.FEASIBLE:
            # 解析求解结果，重构路径序列
            route = [self.depot]
            curr = self.depot
            while len(route) < n:
                for j in range(n):
                    if curr != j and solver.Value(x[curr, j]) == 1:
                        route.append(j)
                        curr = j
                        break
            route.append(self.depot)
            dist = solver.ObjectiveValue() / self.scale_dist
            print(f"[有害垃圾] 精确求解完成，路径长度: {dist:.2f} km")
            return [route]
        return []

    def plot_all_results(self, all_routes_map):
        """可视化多车队清运路径图"""
        plt.figure(figsize=(10, 8))
        
        # 绘制所有基础收集点
        for i in self.nodes:
            plt.scatter(self.coords[i][0], self.coords[i][1], color='silver', s=20, zorder=1)
        plt.scatter(self.coords[self.depot][0], self.coords[self.depot][1], color='black', marker='s', s=100, label='收运中心')

        colors = {'厨余': 'forestgreen', '可回收': 'royalblue', '其他': 'darkorange', '有害': 'crimson'}
        
        for waste, routes in all_routes_map.items():
            c = colors.get(waste, 'black')
            for i, path in enumerate(routes):
                px = [self.coords[node][0] for node in path]
                py = [self.coords[node][1] for node in path]
                # 有害垃圾加粗显示，区分不同性质的作业
                lw = 3 if waste == '有害' else 1.5
                ls = '-' if waste != '有害' else '--'
                label = waste if i == 0 else None
                plt.plot(px, py, color=c, linewidth=lw, linestyle=ls, alpha=0.8, label=label)

        plt.title("多类别废弃物协同清运路径规划示意图", fontsize=14)
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.show()

if __name__ == "__main__":
    optimizer = MultiFleetWasteOptimizer()
    
    # 构建测试需求数据
    test_demands = {i: 1.2 for i in range(1, 31)}
    results = {}
    
    # 策略 1：厨余垃圾收运 (大运量 VRP)
    results['厨余'] = optimizer.solve_standard_vrp("厨余", test_demands, capacity=5, max_dist=100)
    
    # 策略 2：可回收垃圾收运 (大运量 VRP)
    results['可回收'] = optimizer.solve_standard_vrp("可回收", test_demands, capacity=8, max_dist=120)

    # 策略 3：有害垃圾收运 (TSP 精确解)
    results['有害'] = optimizer.solve_hazardous_tsp()

    # 结果展示
    optimizer.plot_all_results(results)
