"""
Flowsheet class
"""

# 导包
import numpy as np
from polymer_model.unit_model.block import Block     # 单元模块
from polymer_model.unit_model.stream import Stream, MaterialStream   # 流股
from polymer_model.utility import SolutionError      # 异常类
import prettytable as pt  # 表格工具, 在终端以表格形式输出结果
from polymer_model.flowsheet import sequence

# 从数据库中导入纯组分和二元交互物性参数
from polymer_model.database.pure_components import args as pure_args
from polymer_model.database.binary_interaction import args as binary_args

# 引入聚合物分布相关函数
from polymer_model.kinetics.ziegler_nat.cstr import get_points


class Flowsheet:
    """
    Define material flow attributes.
    
    Attribute
    ----------


    """
    
    def __init__(self, components, polymers, property_method, dynamic=False):
        self.status = False    # 流程的状态

        self.set_components(components)    # 流程的组分
        self.set_polymers(polymers)        # 设置聚合物相关
        self.property_method = property_method   # 流程所使用的物性方法
        self.get_property_parameters()  # 获取物性参数
        if self.polymers != None:         
            self.get_distribution_setup()  # 获取聚合物分布的配置信息
        self.dynamic = dynamic    # False——稳态, True——动态
        self.blocks : list[Block] = []       # 该集合储存流程中的所有单元模块
        self.streams : list[MaterialStream] = []      # 该集合储存流程中的所有流股


    def set_components(self, components):
        """ 
        设置流程的组分

        Params
        
        """
        self.components = components

        # 创建组分列表
        self.component_list = []  # 小分子和聚合物
        self.polymer_list = []   # 包含聚合物
        self.segment_list = []   # 包含链段
        for component in components:
            if components[component]["type"] == "conventional":
                self.component_list.append(component)
            elif components[component]["type"] == "polymer":
                self.component_list.append(component)
                self.polymer_list.append(component)
            elif components[component]["type"] == "segment":
                self.segment_list.append(component)
    
    def get_components(self):
        """ 获取流程的组分 """
        return self.components

    def set_polymers(self, polymers):
        """ 配置流程的聚合物 """
        self.polymers = polymers
        # 获取催化剂列表
        if polymers != None:
            self.catalyst_list = list(polymers["catalyst"].keys())
        else:
            self.catalyst_list = []

    def set_property_method(self, property_method):
        """ 设置流程的物性方法 """
        self.property_method = property_method
    
    def get_property_method(self):
        """ 获取流程的物性方法 """
        return self.property_method

    def get_property_parameters(self):
        """ 从数据库中获取物性参数 """

        # 物性参数包括：纯组分参数、二元交互参数
        property_parameters = {"Pure Components": {}, }
        for i in self.components:
            # 添加纯组分参数
            property_parameters["Pure Components"][i] = pure_args[i]
            
            # 添加二元交互参数, 还需完善
            # 如果物性方法为“PC-SAFT”, 则需要从数据库中获取二元交互参数PCSKij
            if self.property_method == "PC-SAFT":
                property_parameters["PCSKij"] = dict()
                for j in self.components:
                    if (i, j) in binary_args["PCSKij"]:
                        property_parameters["PCSKij"][(i, j)] = binary_args["PCSKij"][(i,j)]

        self.property_parameters = property_parameters
        return self.property_parameters


    def print_property_parameters(self):
        """ 以表格形式输出物性参数 """

        # 以表格的形式输出
        tb = pt.PrettyTable()

        # 创建列标题(参考Aspen)
        column_header = ["m", "r", "s", "e", "w", "Tc", "Pc", "mw", "enth_mol_form_ig_ref",
                        "entr_mol_form_ig_ref"]   
        # 设置表格的列标题
        tb.add_column(" ", column_header)
        
        # 创建物性数据表
        pure_params = self.property_parameters["Pure Components"]
        for component in self.components:
            column_data = [pure_params[component].get(i, None) for i in column_header]
            # 添加流股的数据
            tb.add_column(component, column_data)

        # 输出表格
        print("Pure component scalar parameters: ")
        print(tb)


    def get_mw(self):
        """ 获取流股各组分的相对分子质量 """

        pure_params = self.property_parameters["Pure Components"]
        NC = len(self.component_list)  # 组分数
        self.mw = np.zeros(NC)
        for i in range(NC):
            self.mw[i] = pure_params[self.component_list[i]]["mw"]
        return self.mw


    def get_distribution_setup(self):
        """ 
        获取聚合物分布图的配置； 横坐标的离散点, 是否执行GPC
        创建一个长度为Np的数组,数值从1离散到upper 
        """
        Np = self.polymers["Distribution"]["Number of points"]  # 点的个数
        upper = self.polymers["Distribution"]["Upper limit"]    # 横坐标上限
        self.DPN_points = get_points(Np, upper)    # 离散的数均聚合度数组(链长分布图的横轴)
        poly_index = self.component_list.index(self.polymer_list[0])
        self.MWN_points = self.DPN_points*self.get_mw()[poly_index]    # 离散的数均分子量数组(分子量分布的横轴)
        self.GPC = self.polymers["Distribution"]["GPC"]  # 是否执行GPC
        return self.DPN_points, self.MWN_points, self.GPC
    

    # 添加单元模块
    def add_block(self, block: Block):
        self.blocks.append(block)
        # 给模块设置组分
        block.set_components(self.components, self.component_list, self.polymer_list, self.segment_list)  
        block.set_property_method(self.property_method)  # 给模块设置物性方法
        block.set_property_parameters(self.property_parameters)  # 给模块设置物性参数
        block.set_polymers(self.polymers, self.catalyst_list)   # 配置聚合物和催化剂
        if self.polymers != None:
            block.set_distribution_setup(self.DPN_points, self.MWN_points, self.GPC)

        
    # 添加流股
    def add_stream(self, stream: Stream):
        self.streams.append(stream)
        stream.set_components(self.components, self.component_list, self.polymer_list, self.segment_list)
        stream.set_property_method(self.property_method)
        stream.set_property_parameters(self.property_parameters)
        stream.set_polymers(self.polymers, self.catalyst_list)   # 配置聚合物和催化剂
        if self.polymers != None:
            stream.set_distribution_setup(self.DPN_points, self.MWN_points, self.GPC)


    # # IDAES添加单元模块的代码类似: fs.b1 = Heater()
    # # 如果使用这样的代码, 需要调用如下的函数来进行配置
    # def get_blocks(self):
    #     """ 获取流程的所有单元模块 """
    #     attr_names = dir(self)    # 获取所有属性名
    #     for name in attr_names:
    #         attr = getattr(self, name)    # 获取属性名对应的值
    #         if isinstance(attr, Block):   # 如果该属性是Block, 则添加到列表中
    #             self.blocks.append(attr)
    #             # 设置单元模块的组分、物性方法和物性参数
    #             attr.set_components(self.components)
    #             attr.set_property_method(self.property_method)
    #             attr.set_property_parameters(self.property_parameters)
    #     return self.blocks
    
    # # IDAES添加流股的代码类似: fs.s1 = Stream()
    # # 如果使用这样的代码, 需要调用如下的函数来进行配置
    # def get_streams(self):
    #     """ 获取流程的所有单元模块 """
    #     attr_names = dir(self)    # 获取所有属性名
    #     for name in attr_names:
    #         attr = getattr(self, name)    # 获取属性名对应的值
    #         if isinstance(attr, Stream):   # 如果该属性是stream, 则添加到列表中
    #             self.streams.append(attr)
    #             # 设置流股的组分、物性方法和物性参数
    #             attr.set_components(self.components)
    #             attr.set_property_method(self.property_method)
    #             attr.set_property_parameters(self.property_parameters)
    #     return self.streams

    
    def run(self, tear_streams=None, calc_sequence=None, estimate=None, 
            convergence_method="Wegstein", tol=0.0001, max_iter=100, options={"q_min":-5, "q_max":0}):
        """ 运行流程计算程序 """

        # 构建流程的邻接矩阵
        n = len(self.blocks)   # 模块的数量
        self.block_matrix = [[0 for i in range(n)] for i in range(n)]
        self.stream_matrix = [[[] for i in range(n)] for i in range(n)]
        for stream in self.streams:
            if stream.source != None and stream.destination != None:
                start = self.blocks.index(stream.source)
                end = self.blocks.index(stream.destination)
                self.block_matrix[start][end] = 1
                self.stream_matrix[start][end].append(stream)
                
        # 计算顺序(其中模块用Block对象表示)
        calc_sequence = sequence.calculate_sequence(self.blocks, self.block_matrix)[0]

        # 撕裂物流(其中物流用起终模块的索引表示, 即[index1, index2], 我们要将其转化为Stream对象)
        tear_streams = sequence.calculate_sequence_index(self.block_matrix)[1]
        tear = []
        for subsystem in tear_streams:
            temp = []
            for stream in subsystem:
                start = stream[0]
                end = stream[1]
                temp.extend(self.stream_matrix[start][end])
            tear.append(temp)
        tear_streams = tear

        # 开始计算
        iter_record = {}   # 迭代记录器, 记录每条撕裂物流的迭代过程
        i = -1
        for subsystem in calc_sequence:
            # 如果遍历到的是模块，则直接计算
            if isinstance(subsystem, Block):
                subsystem.run()
            # 如果遍历到的是列表(即回路)，则需要给撕裂物流赋值，再计算
            else:
                i = i+1 
                # 确定撕裂物流，并且给撕裂物流赋初值
                streams = tear_streams[i]
                for stream in streams:
                    # 如果用户提供了流股的估计值(即estimate中提供了数据), 则使用用户的值
                    if stream in estimate:
                        T = estimate[stream][0]
                        P = estimate[stream][1]
                        mole_flow = estimate[stream][2:]
                    # 如果用户没有提供估计值, 则使用默认值
                    else:
                        T = 298.15
                        P = 1e5
                        mole_flow = np.zeros(len(self.component_list))
                        for i in range(len(self.component_list)):
                            component = self.component_list[i]
                            if component not in self.polymer_list and component not in self.catalyst_list:
                                mole_flow[i] = 1
                    stream.input(T, P, None, None, mole_flow, "Mole & Mole Flow")
                    iter_record[stream] = []
                    temp = np.append(np.array([T, P]), mole_flow)
                    iter_record[stream].append(temp)

                # 开始迭代
                iter = 0
                while True:
                    iter = iter + 1
                    for block in subsystem:
                        block.run()
                    # 计算撕裂物流的误差
                    max = 0
                    for stream in streams:
                        x0 = iter_record[stream][-1]   # 之前的变量
                        print("x0", x0)
                        f0 = np.append(np.array([stream.temperature, stream.pressure]), stream.component_mole_flow)
                        print("cur_var:", f0)
                        iter_record[stream].append(f0)
                        error = (f0-x0)/f0  # 相对误差
                        error[~ np.isfinite(error)] = 0   # 将nan置为0
                        error = np.max(np.abs(error))     # 取绝对值，再取最大值
                        if error > max:
                            max = error
                        print("max", max)
                        print()
                    if max <= tol:
                        break
                    # 如果采用的是直接迭代法,
                    if convergence_method == "Iter":
                        pass
                    # 如果采用的是Wegstein法, 则需要计算新的撕裂流股变量，对撕裂物流进行赋值
                    elif convergence_method == "Wegstein":
                        # q的上下界
                        q_min = options["q_min"]
                        q_max = options["q_max"]

                        if iter == 1:
                            for stream in streams:
                                x1 = iter_record[stream][-1]    # x1 = f0
                                iter_record[stream].append(x1)
                        else:
                            for stream in streams:
                                f1 = iter_record[stream][-1]
                                x1 = iter_record[stream][-2]
                                f0 = iter_record[stream][-3]
                                x0 = iter_record[stream][-4]

                                s = np.array([0 if x1[i] == x0[i] else (f1[i]-f0[i])/(x1[i]-x0[i]) for i in range(x0.size)])
                                s[s==1] = 0
                                q = s/(s-1)
                                # 根据上下限截断q
                                q = np.clip(q, q_min, q_max)
                                x2 = x1 + (1-q)*(f1-x1)
                                x2[x2<0] = 0
                                iter_record[stream].append(x2)

                                T = x2[0]
                                P = x2[1]
                                mole_flow = x2[2:]
                                stream.input(T, P, None, None, mole_flow, "Mole & Mole Flow")
        self.status = True
        

    def print_streams_results(self):
        if self.status == False:
            raise SolutionError("未运行计算程序！")
        
        # 以表格的形式输出
        tb = pt.PrettyTable()
        # 创建列标题(参考Aspen)
        column_header = ["From", "To", "Phase", "Component Mole Flow"]
        column_header.extend(self.component_list)
        column_header.append("Component Mole Fraction")
        column_header.extend(self.component_list)
        column_header.append("Component Mass Flow")
        column_header.extend(self.component_list)
        column_header.append("Component Mass Fraction")
        column_header.extend(self.component_list)
        column_header.extend(["Mole Flow", "Mass Flow", "Volume Flow", 
            "Temperature", "Pressure", "Vapor Fraction", "Liquid Fraction",
            "Mole Enthalpy", "Mass Enthalpy", "Enthalpy Flow", "Mole Entropy",
            "Mass Entropy", "Molar Density", "Mass Density", "Average Molecular Weight"])

        # 设置表格的列标题
        tb.add_column(" ", column_header)
        
        # 添加流股的数据
        for stream in self.streams:
            column_data = ["", "", stream.phase, ""]
            if stream.source != None:
                column_data[0] = stream.source.name
            if stream.destination != None:
                column_data[1] = stream.destination.name
            column_data.extend(stream.component_mole_flow.tolist())
            column_data.append("")
            column_data.extend(stream.mole_fraction.tolist())
            column_data.append("")
            column_data.extend(stream.component_mass_flow.tolist())
            column_data.append("")
            column_data.extend(stream.mass_fraction.tolist())
            column_data.extend([stream.mole_flow, stream.mass_flow, stream.volume_flow,
                stream.temperature, stream.pressure, stream.vapor_fraction, stream.liquid_fraction,
                stream.mole_enthalpy, stream.mass_enthalpy, stream.enthalpy_flow, 
                stream.mole_entropy, stream.mass_entropy, stream.mole_density, 
                stream.mass_density, stream.avg_mw])
            # 添加流股的数据
            tb.add_column(stream.name, column_data)

        # 输出表格
        print("Stream Results: ")
        print(tb)

        for stream in self.streams:
            print("\n",stream.name, "Component Attributes: ")
            print(stream.component_attributes)


# test
if __name__ == "__main__":

    from polymer_model.unit_model.stream import MaterialStream    

    # 流程的组分：C3H6、C3H8、H2、N2、PP、TiCl4、TEA、C3H6-R, H2O
    components = {
        "Propylene": {"type": "conventional"},
        "Propane": {"type": "conventional"},
        "Hydrogen": {"type": "conventional"},
        "Nitrogen": {"type": "conventional"},
        "Polypropylene": {"type": "polymer"},
        "Titanium Tetrachloride": {"type": "conventional"},
        "Triethyl Aluminum": {"type": "conventional"},
        "Propylene-R": {"type": "segment"},
        "Water": {"type": "conventional"},
    }
    # 物性方法
    property_method = "PC-SAFT"
    fs = Flowsheet(components, property_method)
    print(fs.get_property_method())
    
    

    
    fs.set_components(components)
    print(fs.component_list)
    print(fs.polymer_list)
    fs.get_property_parameters()
    fs.print_property_parameters()
      