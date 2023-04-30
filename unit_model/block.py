"""
Block: namely unit model, includes mixer, splitter, flash, heater, CSTR, FPR, pump, compressor/turbine, valve
"""

# Python导入其他目录下的文件
# 方法1: 将工程所在目录添加到环境变量PYTHONPATH中, 然后即可根据 from 工程名.包 import 文件 导入
# 方法2: 将工程所在目录添加到sys.path中, 然后即可根据 from 工程名.包 import 文件 导入
# import sys
# import os
# print(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))))
# sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))))

import numpy as np
from polymer_model.unit_model.stream import Stream
from polymer_model.utility import SolutionError      # 异常类
import prettytable as pt  # 表格工具, 在终端以表格形式输出结果
from typing import List, Set, Tuple, Dict

# 未来可以将Block类定义为抽象类, 防止其实例化
from abc import ABCMeta, abstractmethod

from polymer_model.unit_model.stream import MaterialStream
from polymer_model.kinetics.ziegler_nat.cstr import *    # 引入动力学和反应器模型

class Block:
    """
    Block: namely unit model
    """

    count = 0
    blocks = []
    
    def __init__(self):
        Block.count += 1
        self.name = "B" + str(Block.count)
        Block.blocks.append(self)
        self.status = False      # 单元模块的状态，False为初始状态，True为计算完毕的状态

        # 储存所有的物料流
        self.material_streams = []    
        
        # 暂时不处理热流和功流
        # self.inlet_heat_streams = []        # 储存所有的入口热流
        # self.outlet_heat_streams = []       # 储存所有的出口热流
        # self.inlet_work_streams = []        # 储存所有的入口功流
        # self.outlet_work_streams = []       # 储存所有的出口功流
    
    def set_name(self, name):
        self.name = name
    
    def get_name(self):
        return self.name

    def set_components(self, components, component_list, polymer_list, segment_list):
        self.components = components
        self.component_list = component_list  # 小分子和聚合物
        self.polymer_list = polymer_list   # 包含聚合物
        self.segment_list = segment_list   # 包含链段
        self.number_of_components = len(component_list) # 组分数
        self.number_of_segments = len(segment_list)     # 链段数
    
    def set_polymers(self, polymers, catalyst_list):
        """ 配置聚合物 """
        self.polymers = polymers
        self.catalyst_list = catalyst_list

    def set_property_method(self, property_method):
        """ 设置流程的物性方法 """
        self.property_method = property_method

    def set_property_parameters(self, property_parameters):
        """ 从数据库中获取属性参数 """
        self.property_parameters = property_parameters

    def set_distribution_setup(self, DPN_points, MWN_points, GPC):
        """ 设置聚合物分布的配置 """
        self.DPN_points = DPN_points
        self.MWN_points = MWN_points
        self.GPC = GPC
    
    def get_mw(self):
        """ 获取流股各组分的相对分子质量 """

        pure_params = self.property_parameters["Pure Components"]
        self.mw = np.zeros(self.number_of_components)
        for i in range(self.number_of_components):
            self.mw[i] = pure_params[self.component_list[i]]["mw"]
        return self.mw

    def get_mw_seg(self):
        """ 获取各个链段的相对分子质量 """

        pure_params = self.property_parameters["Pure Components"]
        self.mw_seg = np.zeros(self.number_of_segments)
        for i in range(self.number_of_segments):
            self.mw_seg[i] = pure_params[self.segment_list[i]]["mw"]
        return self.mw_seg


    # 因为模块可以有多个进料流, 计算时先将它们合并为一条流股，再进行各自的计算
    # 所以创建了“mix”方法, 计算混合后的组成、焓、熵
    def mix(self, streams: List[MaterialStream]):
        """ 混合进料流股、计算摩尔流量、组分属性、焓 """

        NC = self.number_of_components   # 组分数

        Fc = np.zeros(NC)             # 各组分的摩尔流量, mol/s
        P_min = streams[0].pressure   # 最小的入口压力 
        H = 0                         # 总焓值, 混合气入口和出口等焓, 单位: W
        T0 = 0                        # 计算流股的流量加权平均温度,作为出口温度计算的初始值
        T_min = streams[0].temperature   # 出口温度的最小值
        T_max = streams[0].temperature   # 出口温度的最大值

        self.get_mw_seg()
        component_attributes = dict()
        distributions = dict()
        # 获取每条流股的流量、组成、焓值, 并汇总
        for stream in streams:
            
            # 如果计算的状态为False, 说明还未计算, 需要调用计算函数
            if stream.status == False:
                stream.run()

            Fc = Fc + stream.component_mole_flow
            H = H + stream.enthalpy_flow
            T0 = T0 + stream.mole_flow*stream.temperature
            if stream.pressure < P_min:
                P_min = stream.pressure
            if stream.temperature > T_max:
                T_max = stream.temperature
            if stream.temperature < T_min:
                T_min = stream.temperature

            # 汇总流股的催化剂各位点流量，聚合物各阶矩流量
            for component in stream.component_attributes:
                # 将各流股催化剂的各位点流量加和
                if component in self.catalyst_list:
                    if component not in component_attributes:
                        component_attributes[component] = dict()
                        component_attributes[component]["CPSFLOW"] = stream.component_attributes[component]["CPSFLOW"]
                        component_attributes[component]["CDSFLOW"] = stream.component_attributes[component]["CDSFLOW"]
                        component_attributes[component]["CISFLOW"] = stream.component_attributes[component]["CISFLOW"]
                        component_attributes[component]["CVSFLOW"] = stream.component_attributes[component]["CVSFLOW"]
                    else:
                        component_attributes[component]["CPSFLOW"] += stream.component_attributes[component]["CPSFLOW"]
                        component_attributes[component]["CDSFLOW"] += stream.component_attributes[component]["CDSFLOW"]
                        component_attributes[component]["CISFLOW"] += stream.component_attributes[component]["CISFLOW"]
                        component_attributes[component]["CVSFLOW"] += stream.component_attributes[component]["CVSFLOW"]

                # 将各流股聚合物的各位的零阶矩和一阶矩流量加和
                if component in self.polymer_list:
                    if component not in component_attributes:
                        component_attributes[component] = dict()
                        distributions[component] = dict()

                        component_attributes[component]["Mass"] = stream.component_attributes[component]["Mass"]
                        distributions[component]["CLD"] = stream.distributions[component]["CLD"]
                        distributions[component]["MWD"] = stream.distributions[component]["MWD"]

                        component_attributes[component]["SFLOW"] = stream.component_attributes[component]["SFLOW"]
                        component_attributes[component]["SFRAC"] = component_attributes[component]["SFLOW"]/np.sum(component_attributes[component]["SFLOW"])
                        mw_seg_avg = np.sum(self.mw_seg*component_attributes[component]["SFRAC"]) 

                        component_attributes[component]["ZMOM"] = stream.component_attributes[component]["ZMOM"]
                        component_attributes[component]["FMOM"] = stream.component_attributes[component]["FMOM"]
                        component_attributes[component]["SMOM"] = stream.component_attributes[component]["SMOM"]
                        component_attributes[component]["DPN"] = component_attributes[component]["FMOM"]/component_attributes[component]["ZMOM"]
                        component_attributes[component]["DPW"] = component_attributes[component]["SMOM"]/component_attributes[component]["FMOM"]
                        component_attributes[component]["MWN"] = component_attributes[component]["DPN"]*mw_seg_avg
                        component_attributes[component]["MWW"] = component_attributes[component]["DPW"]*mw_seg_avg

                        component_attributes[component]["SZMOM"] = stream.component_attributes[component]["SZMOM"]
                        component_attributes[component]["SFMOM"] = stream.component_attributes[component]["SFMOM"]
                        component_attributes[component]["SSMOM"] = stream.component_attributes[component]["SSMOM"]

                        component_attributes[component]["LSZMOM"] = stream.component_attributes[component]["LSZMOM"]
                        component_attributes[component]["DSZMOM"] = stream.component_attributes[component]["DSZMOM"]
                        component_attributes[component]["LSFMOM"] = stream.component_attributes[component]["LSFMOM"]
                        component_attributes[component]["DSFMOM"] = stream.component_attributes[component]["DSFMOM"]
                        component_attributes[component]["LSSMOM"] = stream.component_attributes[component]["LSSMOM"]
                        component_attributes[component]["DSSMOM"] = stream.component_attributes[component]["DSSMOM"]
                    else:
                        distributions[component]["CLD"] = (component_attributes[component]["Mass"]*distributions[component]["CLD"] +
                                stream.component_attributes[component]["Mass"]*stream.distributions[component]["CLD"])/(component_attributes[component]["Mass"] + stream.component_attributes[component]["Mass"])
                        distributions[component]["MWD"] = (component_attributes[component]["Mass"]*distributions[component]["MWD"] +
                                stream.component_attributes[component]["Mass"]*stream.distributions[component]["MWD"])/(component_attributes[component]["Mass"] + stream.component_attributes[component]["Mass"])
                        component_attributes[component]["Mass"] += stream.component_attributes[component]["Mass"]

                        component_attributes[component]["SFLOW"] += stream.component_attributes[component]["SFLOW"]
                        component_attributes[component]["SFRAC"] = component_attributes[component]["SFLOW"]/np.sum(component_attributes[component]["SFLOW"])
                        mw_seg_avg = np.sum(self.mw_seg*component_attributes[component]["SFRAC"]) 

                        component_attributes[component]["ZMOM"] += stream.component_attributes[component]["ZMOM"]
                        component_attributes[component]["FMOM"] += stream.component_attributes[component]["FMOM"]
                        component_attributes[component]["SMOM"] += stream.component_attributes[component]["SMOM"]
                        component_attributes[component]["DPN"] = component_attributes[component]["FMOM"]/component_attributes[component]["ZMOM"]
                        component_attributes[component]["DPW"] = component_attributes[component]["SMOM"]/component_attributes[component]["FMOM"]
                        component_attributes[component]["MWN"] = component_attributes[component]["DPN"]*mw_seg_avg
                        component_attributes[component]["MWW"] = component_attributes[component]["DPW"]*mw_seg_avg

                        component_attributes[component]["SZMOM"] += stream.component_attributes[component]["SZMOM"]
                        component_attributes[component]["SFMOM"] += stream.component_attributes[component]["SFMOM"]
                        component_attributes[component]["SSMOM"] += stream.component_attributes[component]["SSMOM"]

                        component_attributes[component]["LSZMOM"] += stream.component_attributes[component]["LSZMOM"]
                        component_attributes[component]["DSZMOM"] += stream.component_attributes[component]["DSZMOM"]
                        component_attributes[component]["LSFMOM"] += stream.component_attributes[component]["LSFMOM"]
                        component_attributes[component]["DSFMOM"] += stream.component_attributes[component]["DSFMOM"]
                        component_attributes[component]["LSSMOM"] += stream.component_attributes[component]["LSSMOM"]
                        component_attributes[component]["DSSMOM"] += stream.component_attributes[component]["DSSMOM"]
                                    
        P = P_min   # 最小的入口压力

        F = np.sum(Fc)   # 出口总摩尔流量
        T0 = T0/F   # 计算出口温度的初始值

        # 如果用户提供了温度的估计值, 则使用该值作为初始值
        if hasattr(self, "T_guess") and self.T_guess != None:
            T0 = self.T_guess
        
        # 合并为所有进料流股为 一条虚拟的进料流股, 
        feed = MaterialStream()
        feed.set_components(self.components, self.component_list, self.polymer_list, self.segment_list)
        feed.set_polymers(self.polymers, self.catalyst_list)
        feed.set_property_method(self.property_method)
        feed.set_property_parameters(self.property_parameters)
        feed.distributions = distributions

        # 迭代, 更新出口温度的上下界
        while True:
            feed.input(T_min, P, None, None, Fc, "Mole & Mole Flow", component_attributes,
                        self.valid_phases, self.max_iter, self.tol)
            feed.run()
            H_cal = feed.enthalpy_flow
            if H_cal == H:
                self.status = True
                return 
            elif H_cal < H:
                break
            else:
                T_max = T_min
                T_min = 0.9*T_min
        
        # 二分法迭代, 计算出口温度
        error = 2*self.tol
        T = T0
        while True:
            feed.input(T, P, None, None, Fc, "Mole & Mole Flow", component_attributes,
                        self.valid_phases, self.max_iter, self.tol)
            feed.run()
            H_cal = feed.enthalpy_flow
            error = H_cal - H
            if abs(error/H) <= self.tol:
                break
            if error > 0:
                T_max = T
            else:
                T_min = T
            T = T_min + 0.5*(T_max-T_min)

        return feed
    

    def run(self):
        """ 
        Method of calculating inlet and outlet streams attributes and heat duty.
        Any class that inherits Block must implement this method.
        """
        raise NotImplementedError("Method 'run' is not implemented")


    def print_results():
        """ 
        Print results of the unit model.
        Any class that inherits Block must implement this method. 
        """
        raise NotImplementedError("Method 'run' is not implemented")


    def print_stream_results(self):
        """ 输出模块进出口物流信息 """
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
        for stream in self.material_streams:
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

        for stream in self.material_streams:
            print("\n",stream.name, "Component Attributes: ")
            print(stream.component_attributes)

# (输入或者输出)端口
class Port():
    """
    The inlet or outlet of the unit model, 
    which is used to store the inlet streams or outlet streams
    """

    def __init__(self, block: Block) -> None:
        self.block = block   # Module to which it belongs
        self.streams : List[MaterialStream] = []    # store streams
        # self.streams = []    # store streams

    
    def add_stream(self, stream: Stream) -> None:
        self.streams.append(stream)    # 端口添加物流
        self.block.material_streams.append(stream)  # 同时模块也添加物流



# test
if __name__ == "__main__":
    pass









