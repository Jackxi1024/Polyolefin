"""
混合器: 有1个入口和一个出口, 入口可以有多条进料, 出口只有一条出料
混合器的run方法: 
(1) 通过物料衡算, 计算出料流量和组成
(2) 出口压力为进料压力最小值
(3) 通过能量平衡, 计算出料温度
"""

import numpy as np
from polymer_model.unit_model.block import Block, Port

class Mixer(Block):
    """
    Mixer

    Attributes
    ----------
    Pressure : float, default = 0
        Out pressure if value > 0; Pressure drop if value <= 0.
    valid_phases : {"V", "L", "S" "VL", "VLL"}, default = "vapor-liquid" 
    T_guess : , default = None
    max_iter : , default = 30
    tolerance : , default = 0.0001
    feeds : Set of feed stream
    product : Product stream
    """

    def __init__(self):
        super().__init__()
        self.inlet = Port(self)    # 混合器的入口对象
        self.outlet = Port(self)   # 混合器的出口对象


    def input(self, pressure=0, valid_phases="vapor-liquid", T_guess=None, max_iter=100, tol=1e-6):
        self.pressure = pressure
        self.valid_phases = valid_phases
        self.T_guess = T_guess
        self.max_iter = max_iter
        self.tol = tol
    
    def run(self):
        
        NC = self.number_of_components   # 组分数
        streams = self.inlet.streams     # 所有的进料流股

        Fc = np.zeros(NC)         # 各组分的摩尔流量, mol/s
        P_min = streams[0].pressure   # 出口压力等于最小的入口压力 
        H = 0                     # 总焓值, 混合气入口和出口等焓, 单位: W
        T0 = 0                    # 计算流股的流量加权平均温度,作为出口温度计算的初始值
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
                                    
        # 正数代表指定出口压力, 否则指定压降
        if self.pressure <= 0:
            self.pressure = P_min - self.pressure
        P = self.pressure

        F = np.sum(Fc)   # 出口总摩尔流量
        T0 = T0/F   # 计算出口温度的初始值
        # 如果用户提供了温度的估计值, 则使用该值作为初始值
        if self.T_guess != None:
            T0 = self.T_guess
        
        product = self.outlet.streams[0]  # 出料流股
        product.distributions = distributions

        # 迭代, 更新出口温度的上下界
        while True:
            product.input(T_min, P, None, None, Fc, "Mole & Mole Flow", component_attributes,
                        self.valid_phases, self.max_iter, self.tol)
            product.run()
            H_cal = product.enthalpy_flow
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
            product.input(T, P, None, None, Fc, "Mole & Mole Flow", component_attributes,
                        self.valid_phases, self.max_iter, self.tol)
            product.run()
            H_cal = product.enthalpy_flow
            error = H_cal - H
            if abs(error/H) <= self.tol:
                break
            if error > 0:
                T_max = T
            else:
                T_min = T
            T = T_min + 0.5*(T_max-T_min)

        self.status = True

    
    def print_results(self):
        """ 
        Print results of the mixer.
        """

        product = self.outlet.streams[0]  # 出料流股

        print(self.name, "Results: ")
        print("Outlet temperature: ", product.temperature)
        print("Outlet pressure: ", product.pressure)
        print("Vapor fraction: ", product.vapor_fraction)



if __name__ == "__main__":
    from polymer_model.flowsheet.flowsheet import Flowsheet
    from polymer_model.unit_model.stream import MaterialStream
    import time

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

    # 流程的聚合物
    polymers = {
        # 链段
        "segments": {"Propylene-R":{"type": "repeat"},
                    }, 
        # 催化剂
        "catalyst": {"Titanium Tetrachloride": {"type": "Z-N",       # 催化剂类型
                            "site_types_num":4,  # 位点数
                            "site_conc":0.45,    # 位点浓度, mol/kgcat
                            }
                    },
        # 聚合物分子量分布
        "Distribution": {"Number of points":100, "Upper limit":100000, "GPC":True}
    }

    # 物性方法
    property_method = "PC-SAFT"

    # 创建流程
    fs = Flowsheet(components, polymers, property_method)

    # 添加单元模块
    # 添加混合器
    mix1 = Mixer()
    mix1.set_name("Mix1")
    fs.add_block(mix1)

    # 添加流股
    # 混合器进料流股
    N2_feed = MaterialStream(source=None, destination=mix1.inlet)
    N2_feed.set_name("N2 Feed")
    fs.add_stream(N2_feed)

    H2_feed = MaterialStream(source=None, destination=mix1.inlet)
    H2_feed.set_name("H2 Feed")
    fs.add_stream(H2_feed)

    Cat_feed = MaterialStream(source=None, destination=mix1.inlet)
    Cat_feed.set_name("Cat Feed")
    fs.add_stream(Cat_feed)

    Cocat_feed = MaterialStream(source=None, destination=mix1.inlet)
    Cocat_feed.set_name("Cocat Feed")
    fs.add_stream(Cocat_feed)

    C3_feed = MaterialStream(source=None, destination=mix1.inlet)
    C3_feed.set_name("C3 Feed")
    fs.add_stream(C3_feed)

    cycle_gas = MaterialStream(source=None, destination=mix1.inlet)
    cycle_gas.set_name("Cycle Gas")
    fs.add_stream(cycle_gas)

    rfeed1 = MaterialStream(source=mix1.outlet, destination=None)
    rfeed1.set_name("RFeed1")
    fs.add_stream(rfeed1)

    # 设置混合器进料流股参数
    # C3H6、C3H8、H2、N2、PP、TiCl4、TEA、H2O
    z1 = np.array([0, 0, 0, 1, 0, 0, 0, 0])
    N2_feed.input(303.15, 3000000, None, 0.0665563, z1, "Mass & Mass Frac")
    z2 = np.array([0, 0, 1, 0, 0, 0, 0, 0])
    H2_feed.input(303.15, 3000000, None, 0.0235621, z2, "Mass & Mass Frac")
    z3 = np.array([0.988, 0.002, 0, 0, 0, 0.01, 0, 0])
    component_attribute = {
        "Titanium Tetrachloride": {
            "CPSFRAC": 1,
            "CDSFRAC": 0.,
            "CVSFRAC": np.array([0., 0., 0., 0.,]),
            "CISFRAC": np.array([0., 0., 0., 0.,]),

        },
    }
    Cat_feed.input(303.15, 3000000, None, 0.0833333, z3, "Mass & Mass Frac", component_attribute, valid_phases="liquid")
    z4 = np.array([0, 0, 0, 0, 0, 0, 1, 0])
    Cocat_feed.input(303.15, 3000000, None, 0.00277778, z4, "Mass & Mass Frac", valid_phases="liquid")
    z5 = np.array([0.8265628, 0.1734372, 0, 0, 0, 0, 0, 0])
    C3_feed.input(303.15, 3000000, None, 85.01753, z5, "Mass & Mass Frac")
    z6 = np.array([0.8286228, 0.1584892, 0.0057499, 0.00713806, 0, 0, 0, 0])
    cycle_gas.input(334.8567, 3200000, None, 189.1424, z6, "Mass & Mass Frac")


    # 设置各单元模块操作参数
    mix1.input()

    # 运行
    start = time.time()
    mix1.run()

    end = time.time()


    # 输出结果
    mix1.print_results()
    print("\n")
    rfeed1.print_result()

    print("运行时间："+str(end - start)+"秒")