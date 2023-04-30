"""
闪蒸罐: 有1个入口、1个气相出口和1个液相出口, 入口可以有多条进料, 每个出口只有一条出料
闪蒸罐的run方法: 
(1) 物料衡算, 计算所有进料混合后的流量和组成
(2) TP闪蒸: 通过TP闪蒸计算气液相流股组成
    其他闪蒸: 待完成

注意：
如果进料中有聚合物,则需要提供聚合物完整的属性,否则程序会提示一些错误.
如果仅是闪蒸, 则只需提供DPN或者MWN即可计算相平衡, 其余属性随意给出
"""

from dis import dis
import numpy as np
from polymer_model.unit_model.block import Block, Port

class Flash(Block):
    """
    Flash

    Parameters
    ----------
    flash_specs : dict.
        Key: "Temperature"、"Pressure"、"Duty" and "Vapor fraction"
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
            Heat duty of the flash.
    """

    def __init__(self):
        super().__init__()
        self.inlet = Port(self)          # 闪蒸罐的入口
        self.vapor_outlet = Port(self)   # 闪蒸罐的气相出口
        self.liquid_outlet = Port(self)  # 闪蒸罐的液相出口


    def input(self, flash_specs={}, valid_phases="vapor-liquid", T_guess=None, P_guess=None, max_iter=100, tol=1e-7):
        self.flash_specs = flash_specs
        self.valid_phases = valid_phases
        self.T_guess = T_guess
        self.P_guess = P_guess
        self.max_iter = max_iter
        self.tol = tol
    
    def run(self):
        """ 
        Run the flash.
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
                                    
        F = np.sum(Fc)   # 总摩尔流量
        vapor = self.vapor_outlet.streams[0]    # 气相出料流股
        liquid = self.liquid_outlet.streams[0]  # 液相出料流股

        # 出口规范为“温度-压力”
        if self.flash_specs.get("Temperature") != None and self.flash_specs.get("Pressure") != None:
            self.temperature = self.flash_specs.get("Temperature")
            self.pressure = self.flash_specs.get("Pressure")

            # 将vapor作为未分相的出料
            vapor.input(self.temperature, self.pressure, None, None, Fc, "Mole & Mole Flow", 
                        component_attributes, self.valid_phases, self.max_iter, self.tol)
            vapor.run()

            self.vapor_fraction = vapor.vapor_fraction
            # print("beta: ", self.vapor_fraction)
            self.duty = vapor.enthalpy_flow - H_in

            Fl = F*vapor.liquid_fraction      # 液相摩尔流量
            zl = vapor.liquid_mole_fraction   # 液相摩尔分数
            Fv = F*vapor.vapor_fraction       # 气相摩尔流量
            zv = vapor.vapor_mole_fraction    # 气相摩尔分数
            # print("zv: ", zv)
            zv[zv<1e-7] = 0    

            # 分相, 分别给气液流股赋值，以及运行计算程序
            liquid.input(self.temperature, self.pressure, None, Fl, zl, "Mole & Mole Frac", 
                        component_attributes, "liquid", self.max_iter, self.tol)
            liquid.distributions = distributions
            liquid.run()

            vapor.input(self.temperature, self.pressure, None, Fv, zv, "Mole & Mole Frac", 
                        {}, "vapor", self.max_iter, self.tol)
            vapor.run()
            
        # To be completed
        else:
            raise NotImplementedError("The flash specifications has not been implemented!")

        self.status = True

    
    def print_results(self):
        """ 
        Print results of the mixer.
        """

        print(self.name, "Results: ")
        print("Outlet temperature: ",  self.temperature)
        print("Outlet pressure: ",  self.pressure)
        print("Vapor fraction: ", self.vapor_fraction)
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
    flash1 = Flash()
    flash1.set_name("flash1")
    fs.add_block(flash1)

    flash2 = Flash()
    flash2.set_name("flash2")
    fs.add_block(flash2)

    # 添加流股
    # 混合器进料流股
    powder1 = MaterialStream(source=None, destination=flash1.inlet)
    powder1.set_name("Powder1")
    fs.add_stream(powder1)

    gas1 = MaterialStream(source=flash1.vapor_outlet, destination=None)
    gas1.set_name("Gas1")
    fs.add_stream(gas1)

    powder2 = MaterialStream(source=flash1.liquid_outlet, destination=flash2.inlet)
    powder2.set_name("Powder2")
    fs.add_stream(powder2)

    N2 = MaterialStream(source=None, destination=flash2.inlet)
    N2.set_name("N2")
    fs.add_stream(N2)

    gas2 = MaterialStream(source=flash2.vapor_outlet, destination=None)
    gas2.set_name("Gas2")
    fs.add_stream(gas2)

    product = MaterialStream(source=flash2.liquid_outlet, destination=None)
    product.set_name("Product")
    fs.add_stream(product)


    # 设置混合器进料流股参数
    # C3H6、C3H8、H2、N2、PP、TiCl4、TEA、H2O
    z1 = np.array([69.21661, 14.74536, 0.0235413, 0.0665563, 1.138072, 0.000833333, 0.00277773, 0])
    component_attribute = {
        "Titanium Tetrachloride": {
            "CPSFRAC": 0.4649537,
            "CDSFRAC": 0.1223559,
            "CVSFRAC": np.array([2.83208E-05, 0.000229901, 1.01128E-05, 6.95432E-05]),
            "CISFRAC": np.array([0., 0., 0., 0.,]),
        },
        "Polypropylene": {"SFRAC":np.array([1.]), "DPN": 2024.781,
            # 以下属性随意给出
            "Mass": 0,  "DPW": 0, "MWN":0, "MWW":0, "SMOM":0,
            "SZMOM":np.zeros(4), "SFMOM":np.zeros(4), "SSMOM":np.zeros(4), 
            "LSZMOM":np.zeros(4), "LSFMOM":np.zeros(4), "LSSMOM":np.zeros(4), 
            "DSZMOM":np.zeros(4), "DSFMOM":np.zeros(4), "DSSMOM":np.zeros(4)}
    } 
    # 以下属性随意给出
    distributions = {"Polypropylene":{"CLD":np.zeros(100), "MWD":np.zeros(100)}}

    powder1.input(333.15, 3000000, None, None, z1, "Mass & Mass Flow", component_attribute, valid_phases="liquid")
    powder1.distributions = distributions
    z2 = np.array([0, 0, 0, 0.1388889, 0, 0, 0, 0])
    N2.input(303.15, 300000, None, None, z2, "Mass & Mass Flow", valid_phases="vapor")

    # 设置各单元模块操作参数
    flash_specs = {"Temperature": 338.15, "Pressure": 500000}
    flash1.input(flash_specs, max_iter=300)

    flash_specs = {"Temperature": 333.15, "Pressure": 100000}
    flash2.input(flash_specs, max_iter=500)

    # 运行
    start = time.time()
    flash1.run()
    flash2.run()

    end = time.time()

    # 输出结果
    flash1.print_results()
    print("\n")
    flash1.print_stream_results()

    flash2.print_results()
    print("\n")
    flash2.print_stream_results()

    print("运行时间："+str(end - start)+"秒")



