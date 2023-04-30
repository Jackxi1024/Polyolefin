"""
加热器: 有1个入口和1个出口, 入口可以有多条进料, 出口只有一条出料
加热器的run方法: 
(1) 通过物料衡算, 计算出料流量和组成
(2) 当压力、温度指定时, 计算加热器负荷
(3) 当压力、温度不完全指定时, 通过能量衡算, 迭代计算出温度、压力
"""

import numpy as np
from polymer_model.unit_model.block import Block, Port

class Heater(Block):
    """
    Heater

    Parameters
    ----------
    flash_specs : dict.
        Key: "Temperature"、"Pressure"、"Duty"、"Vapor fraction" and so on 
        You need to provide two key-value pairs, usually "temperature" and "pressure"
    valid_phases : str, "liquid"、"vapor" or "vapor-liquid" can be selected. 
                    default = "vapor-liquid" 
    T_guess : float, default = None
    P_guess : float, default = None
    max_iter : float, default = 100
    tol : float, default = 1e-6


    Attributes
    ----------
    Duty : float, 
            Heat duty of the heater.
    """

    def __init__(self):
        super().__init__()
        self.inlet = Port(self)    # 加热器的入口对象
        self.outlet = Port(self)   # 加热器的出口对象

    def input(self, flash_specs={}, valid_phases="vapor-liquid", T_guess=None, P_guess=None, max_iter=100, tol=1e-6):
        self.flash_specs = flash_specs
        self.valid_phases = valid_phases
        self.T_guess = T_guess
        self.P_guess = P_guess
        self.max_iter = max_iter
        self.tol = tol
    
    def run(self):
        """ 
        Run the heater.
        """

        NC = self.number_of_components   # 组分数
        streams = self.inlet.streams     # 所有的进料流股

        Fc = np.zeros(NC)             # 各组分的摩尔流量, mol/s
        P_min = streams[0].pressure   # 获取最小的入口压力, Pa 
        H_in = 0                         # 入口混合气总焓值, W

        self.get_mw_seg()
        component_attributes = dict()
        distributions = dict()
        # 获取每条流股的流量、组成、焓值, 并汇总
        for stream in streams:
            
            # 如果计算的状态为False, 说明还未计算, 需要调用计算函数
            if stream.status == False:
                stream.run()

            Fc = Fc + stream.component_mole_flow
            H_in = H_in + stream.enthalpy_flow
            # 获取最小的入口压力
            if stream.pressure < P_min:
                P_min = stream.pressure

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
                                  
        product = self.outlet.streams[0]  # 出料流股
        product.distributions = distributions

        # 出口规范为“温度-压力”
        if self.flash_specs.get("Temperature") != None and self.flash_specs.get("Pressure") != None:
            self.temperature = self.flash_specs.get("Temperature")
            self.pressure = self.flash_specs.get("Pressure")

            # 正数代表指定出口压力, 否则指定压降
            if self.pressure <= 0:
                self.pressure = P_min - self.pressure
 
            product.input(self.temperature, self.pressure, None, None, Fc, "Mole & Mole Flow", component_attributes,
                        self.valid_phases, self.max_iter, self.tol)
            product.run()
            H_out = product.enthalpy_flow
            self.duty = H_out - H_in

        # To be completed
        else:
            raise NotImplementedError("The flash specifications has not been implemented!")

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
        print("Heat Duty: ", self.duty)


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
    heater = Heater()
    heater.set_name("Heater")
    fs.add_block(heater)

    # 添加流股
    # 混合器进料流股
    S5 = MaterialStream(source=None, destination=heater.inlet)
    S5.set_name("S5")
    fs.add_stream(S5)

    S6 = MaterialStream(source=None, destination=heater.inlet)
    S6.set_name("S6")
    fs.add_stream(S6)

    S7 = MaterialStream(source=heater.outlet, destination=None)
    S7.set_name("S7")
    fs.add_stream(S7)

    # 设置混合器进料流股参数
    # C3H6、C3H8、H2、N2、PP、TiCl4、TEA、H2O
    z1 = np.array([0.5, 0.5, 0, 0, 0, 0, 0, 0])
    S5.input(300, 3000000, None, 1, z1, "Mole & Mole Frac")
    z2 = np.array([0, 0, 0, 0, 0, 1, 0, 0])
    component_attribute = {
        "Titanium Tetrachloride": {
            "CPSFRAC": 1,
            "CDSFRAC": 0.,
            "CVSFRAC": np.array([0., 0., 0., 0.,]),
            "CISFRAC": np.array([0., 0., 0., 0.,]),

        },
    }
    S6.input(300, 3000000, None, 0.01, z2, "Mole & Mole Frac", component_attribute, valid_phases="liquid")


    # 设置各单元模块操作参数
    flash_specs = {"Temperature": 320, "Pressure": 0}
    heater.input(flash_specs)

    # 运行
    start = time.time()
    heater.run()

    end = time.time()


    # 输出结果
    heater.print_results()
    print("\n")
    heater.print_stream_results()

    print("运行时间："+str(end - start)+"秒")

