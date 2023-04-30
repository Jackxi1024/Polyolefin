"""
CSTR反应器: 有1个入口和一个出口, 入口可以有多条进料, 出口只有一条出料
CSTR反应器的run方法: 
(1) 通过物料衡算, 计算出料流量和组成
(2) 出口压力为进料压力最小值
(3) 通过能量平衡, 计算出料温度

"""

import numpy as np
from polymer_model.unit_model.block import Block, Port
from polymer_model.kinetics.ziegler_nat.cstr import *    # 引入动力学和反应器模型
import prettytable as pt  # 表格工具, 在终端以表格形式输出结果


class CSTR(Block):
    """
    CSTR

    Parameters
    ----------
    spec : dict
            Keywords: "Pressure"、"Temperature"、"Duty"、"Vapor fraction"
            Pressure must be provided, and one keyword of the other three
    holdup : dict
            Keywords
            "Valid phases" : 
            "Specification type" : 
            "Reactor volume":
            "Vapor volume":
            "Liquid volume":

    Attributes
    ----------
    duty : float


    """

    def __init__(self):
        super().__init__()
        self.inlet = Port(self)      # 混合器的入口对象
        # 当出口为1个时, 出口为outlet, 出口流股的相态根据valid_phases和计算得出
        # 当出口为2个时, 1个气相出口，1个液相出口
        self.outlet = Port(self)     # 混合器的出口对象
        self.vapor_outlet = Port(self)    # 混合器的气相出口对象
        self.liquid_outlet = Port(self)   # 混合器的液相出口对象


    def input(self, spec={}, holdup={}, reactions={}, flash_params={"Max iter":100, "Tol":1e-7}):
        """ 
        输入
        """
        self.spec = spec
        self.holdup = holdup
        self.reactions = reactions
        self.flash_params = flash_params

    
    def run(self):
        """ 
        运行
        """

        self.get_mw()
        self.get_mw_seg()

        # 如果反应类型为"Ziegler-Nat"
        if self.reactions["Type"] == "Ziegler-Nat":

            NC = self.number_of_components   # 组分数
            NS = self.number_of_segments     # 链段的种类数

            # 确定所用Z-N催化剂活性位点数
            catalyst = self.polymers["catalyst"]
            ns = catalyst[self.catalyst_list[0]]["site_types_num"]  # 位点数目
            
            P_min = 5e10     # 获取最小的入口压力, Pa 
            H_in = 0         # 入口混合气总焓值, W

            # 设置进料
            feed = {}
            feed["Component Mole Flow"] = np.zeros(NC)    # 各组分的摩尔流量, mol/s
            
            # 催化剂各位点流量
            feed["CPSFLOW"] = 0
            feed["CDSFLOW"] = 0
            feed["CISFLOW"] = np.zeros(ns)
            feed["CVSFLOW"] = np.zeros(ns)

            # 获取聚合物各阶矩流量
            feed["SFLOW"] = np.zeros(NS)    # 链段流率
            feed["LSZMOM"] = np.zeros(ns)
            feed["DSZMOM"] = np.zeros(ns)
            feed["LSFMOM"] = np.zeros(ns)
            feed["DSFMOM"] = np.zeros(ns)
            feed["LSSMOM"] = np.zeros(ns)
            feed["DSSMOM"] = np.zeros(ns)

            # 聚合物分布
            feed_polymer_mass = 0    # 进料聚合物总质量kg/s
            feed_mass_cld = np.zeros((ns+1, len(self.DPN_points)))   # 各进料流股(聚合物质量*链长分布)求和
            feed_mass_mwd = np.zeros((ns+1, len(self.MWN_points)))   # 各进料流股(聚合物质量*分子量分布)求和

            # 获取每条流股的流量、焓值、催化剂各位点流量、聚合物各位点各阶矩流量，并汇总
            streams = self.inlet.streams     
            for stream in streams:
                
                # 如果计算的状态为False, 说明还未计算, 需要调用计算函数
                if stream.status == False:
                    stream.run()

                feed["Component Mole Flow"] += stream.component_mole_flow
                H_in = H_in + stream.enthalpy_flow
                # 获取最小的入口压力
                if stream.pressure < P_min:
                    P_min = stream.pressure

                # 汇总流股的催化剂各位点流量，聚合物各阶矩流量
                for component in stream.component_attributes:
                    # 将各流股催化剂的各位点流量加和
                    if component in self.catalyst_list:
                        feed["CPSFLOW"] += stream.component_attributes[component]["CPSFLOW"]
                        feed["CDSFLOW"] += stream.component_attributes[component]["CDSFLOW"]
                        feed["CISFLOW"] += stream.component_attributes[component]["CISFLOW"]
                        feed["CVSFLOW"] += stream.component_attributes[component]["CVSFLOW"]
                    
                    # 将各流股聚合物的各位的零阶矩和一阶矩流量加和
                    if component in self.polymer_list:
                        feed["SFLOW"] += stream.component_attributes[component]["SFLOW"]
                        feed["LSZMOM"] += stream.component_attributes[component]["LSZMOM"]
                        feed["DSZMOM"] += stream.component_attributes[component]["DSZMOM"]
                        feed["LSFMOM"] += stream.component_attributes[component]["LSFMOM"]
                        feed["DSFMOM"] += stream.component_attributes[component]["DSFMOM"]
                        feed["LSSMOM"] += stream.component_attributes[component]["LSSMOM"]
                        feed["DSSMOM"] += stream.component_attributes[component]["DSSMOM"]
                        feed_polymer_mass += stream.component_attributes[component]["Mass"]
                        feed_mass_cld += stream.component_attributes[component]["Mass"]*stream.distributions[component]["CLD"]
                        feed_mass_mwd += stream.component_attributes[component]["Mass"]*stream.distributions[component]["MWD"]

            # 确定压力
            if self.spec["Pressure"] <= 0:
                self.pressure = P_min - self.spec["Pressure"]
            else:
                self.pressure = self.spec["Pressure"]

            if self.spec.get("Temperature") != None:
                self.temperature = self.spec.get("Temperature")
            
                # 反应器参数
                specs = {"Temperature": self.temperature, "Pressure": self.pressure, **self.holdup, 
                        "Reactions": self.reactions, "Flash params":self.flash_params}
                out = cstr_run(self.components, self.property_parameters, feed, specs, self.property_method)
                
                Fc_out = out[0: NC]     # 各组分出口摩尔流量
                Fc_out[Fc_out<0] = 0    # 将负数置为0

                # 出料中各种催化剂位点摩尔流量
                F_ps = out[NC]
                F_ds = out[1+NC]
                F_is = np.array(out[2+NC: 2+NC+ns])
                F_vs = np.array(out[2+NC+ns: 2+NC+2*ns])

                # 出料中聚合物摩尔流量的矩(对应不同的催化位点)   
                LSZMOM = np.array(out[2+NC+2*ns: 2+NC+3*ns])   # 活性链零阶矩
                DSZMOM = np.array(out[2+NC+3*ns: 2+NC+4*ns])   # 死聚物零阶矩
                LSFMOM = np.array(out[2+NC+4*ns: 2+NC+5*ns])   # 活性链一阶矩
                DSFMOM = np.array(out[2+NC+5*ns: 2+NC+6*ns])   # 死聚物一阶矩
                LSSMOM = np.array(out[2+NC+6*ns: 2+NC+7*ns])   # 活性链二阶矩
                DSSMOM = np.array(out[2+NC+7*ns: 2+NC+8*ns])   # 死聚物二阶矩
                SZMOM = LSZMOM + DSZMOM
                SFMOM = LSFMOM + DSFMOM
                SSMOM = LSSMOM + DSSMOM
                SZMOM_gen = SZMOM - (feed["LSZMOM"]+feed["DSZMOM"])   # 生成的零阶矩流量
                SFMOM_gen = SFMOM - (feed["LSFMOM"]+feed["DSFMOM"])   # 生成的一阶矩流量
                SSMOM_gen = SSMOM - (feed["LSSMOM"]+feed["DSSMOM"])   # 生成的二阶矩流量
                ZMOM = np.sum(SZMOM)
                FMOM = np.sum(SFMOM)
                SMOM = np.sum(SSMOM)
                DPN = FMOM/ZMOM
                DPW = SMOM/FMOM

                _, τ, SFLOW = model(out, self.components, self.property_parameters, feed, specs, self.property_method)
                
                SFRAC = SFLOW/np.sum(SFLOW)   # 链段分数
                mw_seg_avg = np.sum(self.mw_seg*SFRAC)  # 链段的平均分子量
                MWN = mw_seg_avg*DPN
                MWW = mw_seg_avg*DPW
                polymer_mass = np.sum(SFLOW*self.mw_seg)/1000     # 出料聚合物质量

                self.site_polymer_mass_gen = SFMOM_gen*mw_seg_avg/1000   # 各位点生成的聚合物质量
                self.polymer_mass_gen = np.sum(self.site_polymer_mass_gen)  # 生成的聚合物质量(各位点综合)
                self.SDPN_gen = SFMOM_gen/SZMOM_gen          # 各位点生成的聚合物数均聚合度
                self.SDPW_gen = SSMOM_gen/SFMOM_gen          # 各位点生成的聚合物重均聚合度
                # 解决0/0问题
                self.SDPN_gen[~ np.isfinite(self.SDPN_gen)] = 0
                self.SDPW_gen[~ np.isfinite(self.SDPW_gen)] = 0
                self.SMWN_gen = self.SDPN_gen*mw_seg_avg     # 各位点生成的聚合物数均分子量
                self.SMWW_gen  = self.SDPW_gen*mw_seg_avg    # 各位点生成的聚合物重均分子量
                self.SPDI_gen = self.SDPW_gen/self.SDPN_gen  # 各位点生成的聚合物分散指数
                # 解决0/0问题
                self.SPDI_gen[~ np.isfinite(self.SPDI_gen)] = 0

                self.DPN_gen = np.sum(SFMOM_gen)/np.sum(SZMOM_gen)  # 生成的聚合物数均聚合度(各位点综合)
                self.DPW_gen = np.sum(SSMOM_gen)/np.sum(SFMOM_gen)  # 生成的聚合物重均聚合度(各位点综合)
                self.MWN_gen = self.DPN_gen*mw_seg_avg     # 生成的聚合物数均分子量(各位点综合)
                self.MWW_gen  = self.DPW_gen*mw_seg_avg    # 生成的聚合物重均分子量(各位点综合)
                self.PDI_gen = self.DPW_gen/self.DPN_gen   # 生成的聚合物分散指数(各位点综合)

                local_cld = cld(τ, SFMOM_gen, self.DPN_points, self.GPC)  # 局部链长分布(局部链长分布的纵轴)
                local_mwd = mwd(τ, SFMOM_gen, mw_seg_avg, self.MWN_points)  # 局部的分子量分布
                local_polymer_mass = polymer_mass - feed_polymer_mass  # 反应器生成的聚合物质量
                cumulative_cld = (feed_mass_cld + local_polymer_mass*local_cld)/polymer_mass  # 累计链长分布
                cumulative_mwd = (feed_mass_mwd + local_polymer_mass*local_mwd)/polymer_mass  # 累计分子量分布
                
                self.local_cld = local_cld
                self.local_mwd = local_mwd
                self.cumulative_cld = cumulative_cld
                self.cumulative_mwd = cumulative_mwd

                catalyst_attributes = {"CPSFLOW": F_ps, "CDSFLOW":F_ds, "CVSFLOW":F_vs, "CISFLOW":F_is}
                
                polymer_attributes = {"Mass": polymer_mass, "SFLOW":SFLOW, "SFRAC":SFRAC, 
                                    "DPN": DPN, "DPW": DPW, "MWN":MWN, "MWW":MWW, 
                                    "ZMOM":ZMOM, "FMOM":FMOM, "SMOM":SMOM,
                                    "SZMOM":SZMOM, "SFMOM":SFMOM, "SSMOM":SSMOM, 
                                    "LSZMOM":LSZMOM, "LSFMOM":LSFMOM, "LSSMOM":LSSMOM, 
                                    "DSZMOM":DSZMOM, "DSFMOM":DSFMOM, "DSSMOM":DSSMOM}
                
                polymer_distributions = {"CLD":cumulative_cld, "MWD":cumulative_mwd}

                attributes = {self.catalyst_list[0]: catalyst_attributes,
                              self.polymer_list[0]: polymer_attributes}
                
                distributions = {self.polymer_list[0]: polymer_distributions}
                
                max_iter = self.flash_params["Max iter"] 
                tol = self.flash_params["Tol"]

                # 出口一条流股
                if len(self.outlet.streams) == 1:
                    product = self.outlet.streams[0]   # 出料流股
                    product.input(self.temperature, self.pressure, None, None, Fc_out, "Mole & Mole Flow",
                                component_attributes=attributes, valid_phases=self.holdup["Valid phases"],
                                max_iter=max_iter, tol=tol)
                    product.distributions = distributions
                    product.run()
                    self.vapor_fraction = product.vapor_fraction
                    self.duty = product.enthalpy_flow - H_in
                    
                    # 计算体积流率
                    self.volume_flow = product.volume_flow
                    self.vapor_volume_flow = product.vapor_volume_flow
                    self.liquid_volume_flow = product.liquid_volume_flow

                # 出口一条气相流股、一条液相流股
                elif len(self.vapor_outlet.streams) == 1 and len(self.liquid_outlet.streams) == 1:
                    vapor_product = self.vapor_outlet.streams[0]
                    liquid_product = self.liquid_outlet.streams[0]

                    # 将liquid作为未分相的出料
                    liquid_product.input(self.temperature, self.pressure, None, None, Fc_out, "Mole & Mole Flow", 
                                component_attributes=attributes, valid_phases=self.holdup["Valid phases"],
                                max_iter=max_iter, tol=tol)
                    liquid_product.run()

                    self.vapor_fraction = liquid_product.vapor_fraction
                    self.duty = liquid_product.enthalpy_flow - H_in

                    # 计算体积流率
                    self.volume_flow = liquid_product.volume_flow
                    self.vapor_volume_flow = liquid_product.vapor_volume_flow
                    self.liquid_volume_flow = liquid_product.liquid_volume_flow

                    Fl = liquid_product.mole_flow*liquid_product.liquid_fraction   # 液相摩尔流量
                    zl = liquid_product.liquid_mole_fraction   # 液相摩尔分数
                    Fv = liquid_product.mole_flow*liquid_product.vapor_fraction    # 气相摩尔流量
                    zv = liquid_product.vapor_mole_fraction    # 气相摩尔分数
                    zv[zv<1e-10] = 0    

                    # 分相, 分别给气液流股赋值，以及运行计算程序
                    liquid_product.input(self.temperature, self.pressure, None, Fl, zl, "Mole & Mole Frac", 
                                attributes, "liquid", max_iter, tol)
                    liquid_product.distributions = distributions
                    liquid_product.run()

                    vapor_product.input(self.temperature, self.pressure, None, Fv, zv, "Mole & Mole Frac", 
                                {}, "vapor", max_iter, tol)
                    vapor_product.run()
                
                # 确定反应器体积和各相体积, 以及停留时间
                # 1.反应器为气液两相, 且指定了各相体积
                if self.holdup["Valid phases"] == "vapor-liquid":
                    if self.holdup["Specification type"] == "Reactor volume & Phase volume":
                        self.reactor_volume = self.holdup["Reactor volume"]   # 反应器体积
                        if "Liquid volume" in self.holdup :
                            self.liquid_volume = self.holdup["Liquid volume"]   # 液相体积
                            self.vapor_volume = self.reactor_volume - self.liquid_volume  # 气相体积
                        elif "Vapor volume" in self.holdup :
                            self.vapor_volume = self.holdup["Vapor volume"]    # 气相体积  
                            self.liquid_volume = self.reactor_volume - self.vapor_volume  # 液相体积
                        # 计算停留时间
                        self.reactor_residence_time = self.reactor_volume/self.volume_flow
                        self.vapor_residence_time = self.vapor_volume/self.vapor_volume_flow
                        self.liquid_residence_time = self.liquid_volume/self.liquid_volume_flow
                # 2、如果CSTR内只有液相
                elif self.holdup["Valid phases"] == "liquid":
                    if self.holdup["Specification type"] == "Reactor volume":
                        self.reactor_volume = self.holdup["Reactor volume"]   # 反应器体积
                        self.liquid_volume = self.reactor_volume   # 液相体积
                        self.vapor_volume = None
                        # 计算停留时间
                        self.reactor_residence_time = self.reactor_volume/self.volume_flow
                        self.liquid_residence_time = self.liquid_volume/self.liquid_volume_flow
                        self.vapor_residence_time = None
                elif self.holdup["Valid phases"] == "vapor":
                    if self.holdup["Specification type"] == "Reactor volume":
                        self.reactor_volume = self.holdup["Reactor volume"]   # 反应器体积
                        self.vapor_volume = self.reactor_volume   # 液相体积
                        self.liquid_volume = None
                        # 计算停留时间
                        self.reactor_residence_time = self.reactor_volume/self.volume_flow
                        self.vapor_residence_time = self.vapor_volume/self.vapor_volume_flow
                        self.liquid_residence_time = None
            
        else:
            raise NotImplementedError("The specfication has not been implemented!")
        
        self.status = True

     
    def print_results(self):
        """ 
        Print results of the cstr.
        """
        
        print(self.name, "Results: ")
        print("Outlet temperature:\t", self.temperature)
        print("Outlet pressure:\t", self.pressure)
        print("Vapor fraction:\t", self.vapor_fraction)
        print("Heat duty:\t", self.duty)
        print("Volume")
        print("Reactor:\t", self.reactor_volume)
        print("Vapor phase:\t", self.vapor_volume)
        print("Liquid phase:\t", self.liquid_volume)
        print("Residence time")
        print("Reactor:\t", self.reactor_residence_time)
        print("Vapor phase:\t", self.vapor_residence_time)
        print("Liquid phase:\t", self.liquid_residence_time)
    

    def print_local_polymer_attributes(self):
        """ 
        以表格形式输出反应器生成的聚合物属性:
        各位点生成的聚合物 Mass、DPN、DPW、MWN、MWW、PDI、MWS
        """

        tb = pt.PrettyTable()

        # 设置表格的第一行
        tb.field_names = ["", "Poly Mass Generated", "DPN", "DPW", "MWN", "MWW", "PDI"]
        ns = len(self.SDPN_gen)
        for i in range(ns):
            data = ["Site"+str(i+1), self.site_polymer_mass_gen[i], self.SDPN_gen[i], self.SDPW_gen[i],
                    self.SMWN_gen[i], self.SMWW_gen[i], self.SPDI_gen[i]]
            tb.add_row(data)
        tb.add_row(["Composite", self.polymer_mass_gen, self.DPN_gen, self.DPW_gen,
                    self.MWN_gen, self.MWW_gen, self.PDI_gen])

        # 输出表格
        print(self.name," Generated Polymer Attributes")
        print(tb)


    def print_local_cld(self):
        """ 以表格形式输出局部链长分布"""
        
        tb = pt.PrettyTable()

        # 设置表格的第一列
        tb.add_column("DPN", self.DPN_points)
        ns = len(self.local_cld)-1    # 位点数
        for i in range(1, ns+1):
            tb.add_column("Site_"+str(i), self.local_cld[i])
        tb.add_column("Composite", self.local_cld[0])

        # 输出表格
        print(self.name," Site-based Chain Length Distribution")
        print(tb)


    def print_local_mwd(self):
        """ 以表格形式输出局部分子量分布 """

        tb = pt.PrettyTable()

        # 设置表格的第一列
        tb.add_column("DPN", self.MWN_points)
        ns = len(self.local_mwd)-1    # 位点数
        for i in range(1, ns+1):
            tb.add_column("Site_"+str(i), self.local_mwd[i])
        tb.add_column("Composite", self.local_mwd[0])

        print(self.name," Site-based Molecular Weight Distribution")
        print(tb)


    def plot_local_cld(self):
        " Plot local cld"

        window_title = self.name
        figure_title = "Local Chain Length Distribution"
        xlabel = "log(n)"
        if self.GPC:
            ylabel = "nw(n)"
        else:
            ylabel = "w(n)"
        distribution_plot(self.DPN_points, self.local_cld, window_title, figure_title, xlabel, ylabel)


    def plot_cumu_cld(self):
        " Plot cumulative cld"

        window_title = self.name
        figure_title = "Cumulative Chain Length Distribution"
        xlabel = "log(n)"
        if self.GPC:
            ylabel = "nw(n)"
        else:
            ylabel = "w(n)"
        distribution_plot(self.DPN_points, self.cumulative_cld, window_title, figure_title, xlabel, ylabel)


    def plot_local_mwd(self):
        " Plot local mwd"

        window_title = self.name
        figure_title = "Local Molecular Weight Distribution"
        xlabel = "log(Mn)"
        ylabel = "W(log(Mn))"
        distribution_plot(self.MWN_points, self.local_mwd, window_title, figure_title, xlabel, ylabel)


    def plot_cumu_mwd(self):
        " Plot cumulative mwd"

        window_title = self.name
        figure_title = "Cumulative Molecular Weight Distribution"
        xlabel = "log(Mn)"
        ylabel = "W(log(Mn))"
        distribution_plot(self.MWN_points, self.cumulative_mwd, window_title, figure_title, xlabel, ylabel)



if __name__ == "__main__":
    from polymer_model.flowsheet.flowsheet import Flowsheet
    from polymer_model.unit_model.stream import MaterialStream
    import time
    import numpy as np

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
    # 添加反应器
    cstr1 = CSTR()
    cstr1.set_name("cstr1")
    fs.add_block(cstr1)

    # 添加流股
    feed1 = MaterialStream(source=None, destination=cstr1.inlet)
    feed1.set_name("feed1")
    fs.add_stream(feed1)

    vap = MaterialStream(source=cstr1.vapor_outlet, destination=None)
    vap.set_name("vap")
    fs.add_stream(vap)

    powder = MaterialStream(source=cstr1.liquid_outlet, destination=None)
    powder.set_name("powder")
    fs.add_stream(powder)

    # 设置反应器进料流股参数
    # C3H6、C3H8、H2、N2、PP、TiCl4、TEA、H2O
    z1 = np.array([0.7695529, 0.1446299, 0.0786014, 0.00721171, 0, 6.26484E-07, 3.46973E-06, 0])
    component_attribute = {
        "Titanium Tetrachloride": {
            "CPSFRAC": 1.,
            "CDSFRAC": 0.,
            "CVSFRAC": np.array([0., 0., 0., 0.,]),
            "CISFRAC": np.array([0., 0., 0., 0.,]),
        },
    }
    feed1.input(330.546, 3000000, None, 7012.333, z1, "Mole & Mole Frac", component_attribute)

    # 设置反应
    # 设置Ziegler-Natta反应物种
    species = { "polymer": "Polypropylene",          # 聚合物
                "tdb segment": None,      # 终端双键(用于计算支链)
                "monomers": ["Propylene"],       # 单体
                "segments": {"Propylene":"Propylene-R"},     # 链段
                "precatalyst": None,      # 预催化剂
                "catalyst": ["Titanium Tetrachloride"],      # 催化剂
                "cocatalysts": ["Triethyl Aluminum"],     # 助催化剂
                "solvents": None,         # 溶剂
                "transfer agent": None,   # 链转移剂
                "hydrogens": ["Hydrogen"],  # 氢
                "poisons": None,          # 毒物
                "elec don": None,         # 电子供体
                "byproduct": None,        # 副产物
            }
    
    # 设置Z-N反应：反应类型, 催化位点, 组分1, 组分2, 前指因子, 活化能, 反应级数, 终端双键分数, 参考温度
    r1 = [["Act-Spon", 1, "Titanium Tetrachloride", None, 0.0013, 3.19872e4, 1, None, 343.15],
        ["Act-Spon", 2, "Titanium Tetrachloride", None, 0.0013, 3.19872e4, 1, None, 343.15],
        ["Act-Spon", 3, "Titanium Tetrachloride", None, 0.0013, 3.19872e4, 1, None, 343.15],
        ["Act-Spon", 4, "Titanium Tetrachloride", None, 0.0013, 3.19872e4, 1, None, 343.15],
        ["Chain-Ini", 1, "Propylene", None, 108.85/1000, 3.0145e4, 1, None, 343.15],
        ["Chain-Ini", 2, "Propylene", None, 24.5/1000,   3.0145e4, 1, None, 343.15],
        ["Chain-Ini", 3, "Propylene", None, 170.8/1000,  3.0145e4, 1, None, 343.15],
        ["Chain-Ini", 4, "Propylene", None, 60.55/1000,  3.0145e4, 1, None, 343.15],
        ["Propagation", 1, "Propylene-R", "Propylene", 108.85/1000,  3.0145e4, 1, None, 343.15],
        ["Propagation", 2, "Propylene-R", "Propylene", 24.5/1000,  3.0145e4, 1, None, 343.15],
        ["Propagation", 3, "Propylene-R", "Propylene", 170.8/1000,  3.0145e4, 1, None, 343.15],
        ["Propagation", 4, "Propylene-R", "Propylene", 60.55/1000,  3.0145e4, 1, None, 343.15],
        ["Chat-Mon", 1, "Propylene-R", "Propylene", 0.012/1000,  5.2e4, 1, None, 343.15],
        ["Chat-Mon", 2, "Propylene-R", "Propylene", 0.012/1000,  5.2e4, 1, None, 343.15],
        ["Chat-Mon", 3, "Propylene-R", "Propylene", 0.012/1000,  5.2e4, 1, None, 343.15],
        ["Chat-Mon", 4, "Propylene-R", "Propylene", 0.012/1000,  5.2e4, 1, None, 343.15],
        ["Chat-Cocat", 1, "Propylene-R", "Triethyl Aluminum", 0.12/1000,  5.02416e4, 1, None, 343.15],
        ["Chat-Cocat", 2, "Propylene-R", "Triethyl Aluminum", 0.12/1000,  5.02416e4, 1, None, 343.15],
        ["Chat-Cocat", 3, "Propylene-R", "Triethyl Aluminum", 0.12/1000,  5.02416e4, 1, None, 343.15],
        ["Chat-Cocat", 4, "Propylene-R", "Triethyl Aluminum", 0.12/1000,  5.02416e4, 1, None, 343.15],
        ["Chat-H2", 1, "Propylene-R", "Hydrogen", 4.8/1000,  4.47988e4, 1, None, 343.15],
        ["Chat-H2", 2, "Propylene-R", "Hydrogen", 8.88/1000,  4.47988e4, 1, None, 343.15],
        ["Chat-H2", 3, "Propylene-R", "Hydrogen", 2.64/1000,  4.47988e4, 1, None, 343.15],
        ["Chat-H2", 4, "Propylene-R", "Hydrogen", 6.6/1000,  4.47988e4, 1, None, 343.15],
        ["Deact-Spon", 1, None, None, 0.001,  4.1868e3, 1, None, 343.15],
        ["Deact-Spon", 2, None, None, 0.001,  4.1868e3, 1, None, 343.15],
        ["Deact-Spon", 3, None, None, 0.001,  4.1868e3, 1, None, 343.15],
        ["Deact-Spon", 4, None, None, 0.001,  4.1868e3, 1, None, 343.15],
        ["Atact-Prop", 1, "Propylene-R", "Propylene", 8/1000,  3.0145e4, 1, None, 343.15],
        ["Atact-Prop", 2, "Propylene-R", "Propylene", 1/1000,  3.0145e4, 1, None, 343.15],
        ["Atact-Prop", 3, "Propylene-R", "Propylene", 0.1/1000,  3.0145e4, 1, None, 343.15],
        ["Atact-Prop", 4, "Propylene-R", "Propylene", 0.1/1000,  3.0145e4, 1, None, 343.15]]

    # 设置反应
    reaction={"Type":"Ziegler-Nat", 
              "Species": species, 
              "Reactions": r1, 
              "Reacting phase": "liquid"}

    # 设置反应器体积和各相体积(或者反应器停留时间、各相停留时间)
    holdup = {"Valid phases": "vapor-liquid",
            "Specification type": "Reactor volume & Phase volume", 
            "Reactor volume": 90, 
            "Liquid volume": 60}
    

    # 设置反应器
    cstr1.input(spec={"Pressure":3e6 ,"Temperature":333.15}, holdup=holdup, reactions=reaction,
                flash_params={"Max iter":100, "Tol":1e-7})

    # 运行
    start = time.time()
    cstr1.run()
    end = time.time()


    # 输出结果
    cstr1.print_results()
    print("\n")
    cstr1.print_stream_results()

    # 反应器生成的聚合物属性
    cstr1.print_local_polymer_attributes()

    # 反应器生成的聚合物链长分布和聚合物分布(局部)
    cstr1.print_local_cld()
    cstr1.print_local_mwd()

    # 反应器出料的聚合物链长分布和聚合物分布(累计)
    powder.print_cld()
    powder.print_mwd()

    # 绘图
    cstr1.plot_local_cld()
    cstr1.plot_local_mwd()
    powder.plot_cld()
    powder.plot_mwd()

    print("运行时间："+str(end - start)+"秒")
