"""
Define stream and various attributes
"""

import numpy as np
from polymer_model.utility import SolutionError, InputError      # 异常类
import prettytable as pt    # 表格工具, 在终端以表格形式输出结果
from copy import deepcopy   # 深度复制
from polymer_model.thermodynamics.thermo_property import enth_ig_comp, entr_ig_comp, vapor_pressure
from polymer_model.property_method.property_method_wrapper import PropertyMethod
from polymer_model.thermodynamics.flash import TP_flash 
from polymer_model.utility import R   # 气体常数
# 引入聚合物分布相关函数
from polymer_model.kinetics.ziegler_nat.cstr import get_points, distribution_plot


class Stream:
    """
    Abstract flow, including material flow, work flow and heat flow

    Attribute
    ----------
    name : str
    source : Block
    destination : Block
    """

    count = 0
    streams = []
    
    def __init__(self, source=None, destination=None):
        Stream.count += 1
        self.name = "S" + str(Stream.count)  # 流股的默认名称
        Stream.streams.append(self)
        self.set_source(source)              # 设置流股的来源
        self.set_destination(destination)    # 设置流股的去处
        self.status = False      # 流股的状态，False为初始状态，True为计算完毕的状态


    def set_name(self, name):
        self.name = name
    
    def get_name(self):
        return self.name

    def set_destination(self, block_inlet):
        if block_inlet == None:
            self.destination = None
        else:
            self.destination = block_inlet.block
            block_inlet.add_stream(self) 
        
    def get_destination(self):
        return self.destination

    def set_source(self, block_outlet):
        if block_outlet == None:
            self.source = None
        else:
            self.source = block_outlet.block
            block_outlet.add_stream(self)


class HeatStream(Stream):
    """
    Heat Stream
    """

    def set_duty(self, duty):
        self.duty = duty
    
    def get_duty(self):
        return self.duty

    def set_temperature(self, start_temperature, end_temperature):
        self.start_temperature = start_temperature
        self.end_temperature = end_temperature 

    def get_temperature(self):
        return self.start_temperature, self.end_temperature


class WorkStream(Stream):
    """
    Work Stream
    """

    def set_power(self, power):
        self.power = power
    
    def get_power(self):
        return self.power

    def set_speed(self, speed):
        self.speed = speed

    def get_speed(self):
        return self.speed


class MaterialStream(Stream):
    """
    Define material flow attributes.
    
    Attribute
    ----------
    temperature : float
        stream temperature (K)
    pressure : float
        stream pressure (Pa)
    vapor_fraction : float
        stream mole vapor fraction. Enter None if temperature and pressure are provided.
    flow_basis : str, "mole"/"mass". 
    flow_rate : float, 
        stream flow rate (mol/s or kg/s)
    components : ndarray, shape (n,)
        Names of all components in the system.
    composition : ndarray, shape (n,)
        Mole/Mass fractions of each component in the system. The sum is 1.

    """

    def __init__(self, source=None, destination=None):
        super().__init__(source, destination)
        self.components = None 
        self.property_method = None
        self.property_parameters = None
        self.property_args = None
        self.component_attributes = {}  # 组分属性
        self.distributions = {}     # 聚合物分布
        
    
    def set_components(self, components, component_list, polymer_list, segment_list): 
        """ 配置组分 """
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
    
    def get_property_method(self):
        """ 获取流程的物性方法 """
        return self.property_method

    def set_property_parameters(self, property_parameters):
        """ 设置物性参数 """
        self.property_parameters = property_parameters 

    def get_property_parameters(self):
        """ 获取物性参数 """
        return self.property_parameters 


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
    
    
    def input(self, temperature, pressure, vapor_fraction, flow_rate, composition,
            flow_basis="Mole & Mole Frac", component_attributes=dict(),
            valid_phases="vapor-liquid", max_iter=100, tol=1e-7):
        """
        输入流股的信息
        温度、压力、汽化分数需要给定其中2个, 目前仅支持给定温度和压力
        flow_basis

        Params
        ----------
        temperature : 
        pressure :
        vapor_fraction :
        flow_rate :
        composition :
        flow_basis : 包括"Mole & Mole Frac"、"Mole & Mole Flow"、"Mass & Mass Frac"、"Mass & Mass Flow",
                分别是给定"总摩尔流量和各组分摩尔分数"、"各组分摩尔流量"、"总质量流量和各组分质量分数"、
                "各组分质量流量"
        component_attributes : dict, option
            如果组分中有催化剂, 则需要输入催化剂各位点流量/分数
            如果组分中有聚合物, 则需要输入聚合物链段分数、DPN
            Examples:
            component_attributes = {
                "Titanium Tetrachloride":{
                    "CPSFLOW": 0.000375, 
                    "CDSFLOW": 0., 
                    "CISFLOW": np.array([0., 0., 0., 0.]), 
                    "CVSFLOW": np.array([0., 0., 0., 0.]), 
                }
                "Polypropylene":{
                    "SFRAC": np.array([1.0, 0., 0., 0.])
                    "DPN": 2000
                }

            }
        valid_phases : 
        """

        # 通常流股在计算完毕后, 会将state置为True, 避免重复计算
        # 但在循环流程中, 一条流股可能进行多次赋值和计算
        # 因此在赋值时，需要将state变为False
        self.status = False    

        self.temperature = temperature
        self.pressure = pressure
        self.vapor_fraction = vapor_fraction

        mw = self.get_mw()   # 各组分相对分子质量
        mw_seg = self.get_mw_seg()   # 各链段的相对分子质量
        # 给定总摩尔流量和各组分摩尔分数
        if flow_basis == "Mole & Mole Frac":
            self.mole_flow = flow_rate
            self.mole_fraction = composition
            self.component_mole_flow = flow_rate*composition

            self.component_mass_flow = self.component_mole_flow*mw/1000    # Kg/s
            self.mass_flow = np.sum(self.component_mass_flow)
            self.mass_fraction = self.component_mass_flow/self.mass_flow
        # 给定各组分摩尔流量
        elif flow_basis == "Mole & Mole Flow":
            self.component_mole_flow = composition
            self.mole_flow = np.sum(self.component_mole_flow)
            self.mole_fraction = self.component_mole_flow/self.mole_flow  

            self.component_mass_flow = self.component_mole_flow*mw/1000    # Kg/s
            self.mass_flow = np.sum(self.component_mass_flow)
            self.mass_fraction = self.component_mass_flow/self.mass_flow
        # 给定总质量流率和各组分质量分数
        elif flow_basis == "Mass & Mass Frac":
            self.mass_flow = flow_rate
            self.mass_fraction = composition
            self.component_mass_flow = flow_rate*composition
            
            self.component_mole_flow = self.component_mass_flow/mw*1000   # mol/s
            self.mole_flow = np.sum(self.component_mole_flow)
            self.mole_fraction = self.component_mole_flow/self.mole_flow 
        elif flow_basis == "Mass & Mass Flow":
            self.component_mass_flow = composition
            self.mass_flow = np.sum(self.component_mass_flow)
            self.mass_fraction = self.component_mass_flow/self.mass_flow

            self.component_mole_flow = self.component_mass_flow/mw*1000   # mol/s
            self.mole_flow = np.sum(self.component_mole_flow)
            self.mole_fraction = self.component_mole_flow/self.mole_flow 
        self.avg_mw = np.sum(self.mole_fraction*mw) 

        # 获取组分属性: 催化剂各位点数据、聚合物的数均分子量
        self.component_attributes = component_attributes
        for index in range(self.number_of_components):
            component = self.component_list[index]
            # 如果该组分为催化剂, 并且含量非0, 则需要获取各个位点的流量和分数
            if component in self.catalyst_list and self.component_mole_flow[index] != 0:
                site_conc = self.polymers["catalyst"][component]["site_conc"]  # 位点浓度, mol/kgcat
                site_mole_flow = self.component_mass_flow[index]*site_conc     # 位点流量, mol/s
                # 获取潜在位点的流量和摩尔分数
                if "CPSFRAC" in component_attributes[component] and "CPSFLOW" not in component_attributes[component]:
                    component_attributes[component]["CPSFLOW"] = component_attributes[component]["CPSFRAC"]*site_mole_flow
                elif "CPSFLOW" in component_attributes[component] and "CPSFRAC" not in component_attributes[component]:
                    component_attributes[component]["CPSFRAC"] = component_attributes[component]["CPSFLOW"]/site_mole_flow   
                # 获取死亡位点的流量和摩尔分数
                if "CDSFRAC" in component_attributes[component] and "CDSFLOW" not in component_attributes[component]:
                    component_attributes[component]["CDSFLOW"] = component_attributes[component]["CDSFRAC"]*site_mole_flow
                elif "CDSFLOW" in component_attributes[component] and "CDSFRAC" not in component_attributes[component]:
                    component_attributes[component]["CDSFRAC"] = component_attributes[component]["CDSFLOW"]/site_mole_flow 
                # 获取空位点的流量和摩尔分数
                if "CVSFRAC" in component_attributes[component] and "CVSFLOW" not in component_attributes[component]:
                    component_attributes[component]["CVSFLOW"] = component_attributes[component]["CVSFRAC"]*site_mole_flow
                elif "CVSFLOW" in component_attributes[component] and "CVSFRAC" not in component_attributes[component]:
                    component_attributes[component]["CVSFRAC"] = component_attributes[component]["CVSFLOW"]/site_mole_flow
                # 获取抑制位点的流量和摩尔分数
                if "CISFRAC" in component_attributes[component] and "CISFLOW" not in component_attributes[component]:
                    component_attributes[component]["CISFLOW"] = component_attributes[component]["CISFRAC"]*site_mole_flow
                elif "CISFLOW" in component_attributes[component] and "CISFRAC" not in component_attributes[component]:
                    component_attributes[component]["CISFRAC"] = component_attributes[component]["CISFLOW"]/site_mole_flow 
            
            # 如果该组分为聚合物, 并且含量非0, 则需要获取链段分数和数均聚合度
            elif component in self.polymer_list and self.component_mole_flow[index] != 0:
                # 获取总的链段流率(即FMOM)
                if "FMOM" not in component_attributes[component]:
                    component_attributes[component]["FMOM"] = self.component_mole_flow[index]
                # 获取链段分数和链段流率
                if "SFRAC" in component_attributes[component] and "SFLOW" not in component_attributes[component]:
                    component_attributes[component]["SFLOW"] = component_attributes[component]["SFRAC"]*component_attributes[component]["FMOM"]
                elif "SFLOW" in component_attributes[component] and "SFRAC" not in component_attributes[component]:
                    component_attributes[component]["SFRAC"] = component_attributes[component]["SFLOW"]/component_attributes[component]["FMOM"]   
                # 计算聚合物链段平均分子质量
                mw_seg_avg = np.sum(component_attributes[component]["SFRAC"]*mw_seg) 
                # 获取数均聚合度和数均分子量
                if "DPN" in component_attributes[component] and "MWN" not in component_attributes[component]:
                    component_attributes[component]["MWN"] = component_attributes[component]["DPN"]*mw_seg_avg
                elif "MWN" in component_attributes[component] and "DPN" not in component_attributes[component]:
                    component_attributes[component]["DPN"] = component_attributes[component]["MWN"]/mw_seg_avg  

                # 计算ZMOM
                if "ZMOM" not in component_attributes[component]:
                    component_attributes[component]["ZMOM"] = component_attributes[component]["FMOM"]/component_attributes[component]["DPN"]

        # 用于闪蒸计算的参数
        self.valid_phases = valid_phases
        self.max_iter = max_iter
        self.tol = tol

        
    def convert_params_to_args(self):
        """ 转化参数 """

        # 纯组分参数
        pure_params = self.property_parameters["Pure Components"]
        
        NC = self.number_of_components   # 组分数
        NS = self.number_of_segments     # 链段数

        # 如果物性方法为"PC-SAFT", 则需要从物性参数中获取链段数m, 链段直径s, 能量参数e等参数
        if self.property_method == "PC-SAFT":

            m = np.zeros(NC)      # 各组分链段数, Aspen用PCSFTM表示
            s = np.zeros(NC)      # 各组分链段直径, Aspen用PCSFTV表示
            e = np.zeros(NC)      # 各组分能量参数, Aspen用PCSFTU表示
            e_assoc = np.zeros(NC)  # 关联成分的关联能量,单位K, Aspen用PCSFAU表示
            vol_a = np.zeros(NC)    # 关联成分的有效关联量, Aspen用PCSFAV表示
            w = np.zeros(NC)      # 各组分偏心因子
            Tc = np.zeros(NC)     # 各组分临界温度
            Pc = np.zeros(NC)     # 各组分临界压力

            for i in range(NC):
                component = self.component_list[i]
                
                w[i] = pure_params[component]["w"]
                Tc[i] = pure_params[component]["Tc"]
                Pc[i] = pure_params[component]["Pc"]
                if self.components[component]["type"] == "conventional":
                    m[i] = pure_params[component]["m"]
                    s[i] = pure_params[component]["s"]
                    e[i] = pure_params[component]["e"]
                    if pure_params[component].get("e_assoc") != None:
                        e_assoc[i] = pure_params[component].get("e_assoc")
                    else:
                        e_assoc[i] = 0
                    if pure_params[component].get("vol_a") != None:
                        vol_a[i] = pure_params[component].get("vol_a")
                    else:
                        vol_a[i] = 0
                elif self.components[component]["type"] == "polymer":
                    if self.component_attributes == None or component not in self.component_attributes:
                        m[i] = 1
                        s[i] = 3
                        e[i] = 200
                    else:
                        DPN = self.component_attributes[component]["DPN"]

                        # 聚合物的性质由其组成的链段性质计算
                        sfrac = self.component_attributes[component]["SFRAC"]   # 各链段分数
                        mw_seg = np.zeros(NS)  # 各链段分子质量
                        r_seg = np.zeros(NS)   # 各链段的链段数/数均分子量
                        s_seg = np.zeros(NS)   # 各链段的链段直径
                        e_seg = np.zeros(NS)   # 各链段的能量参数
                        e_assoc_seg = np.zeros(NS)
                        vol_a_seg = np.zeros(NS)

                        for j in range(NS):
                            segment = self.segment_list[j]
                            mw_seg[j] = pure_params[segment]["mw"]
                            r_seg[j] = pure_params[segment]["r"]
                            s_seg[j] = pure_params[segment]["s"]
                            e_seg[j] = pure_params[segment]["e"]
                            if pure_params[segment].get("e_assoc") != None:
                                e_assoc_seg[j] = pure_params[segment].get("e_assoc")
                            else:
                                e_assoc_seg[j] = 0
                            if pure_params[segment].get("vol_a") != None:
                                vol_a_seg[j] = pure_params[segment].get("vol_a")
                            else:
                                vol_a_seg[j] = 0

                        MWN = sfrac*mw_seg*DPN
                        m[i] = np.sum(r_seg*MWN) 
                        s[i] = np.sum(sfrac*s_seg)
                        e[i] = np.sum(sfrac*e_seg)
                        e_assoc[i] = np.sum(sfrac*e_assoc_seg)
                        vol_a[i] = np.sum(sfrac*vol_a_seg)
                
            property_args = {"Tc": Tc, "Pc": Pc, "w": w, "m": m, "s": s, "e": e, 
                            "e_assoc": e_assoc, "vol_a": vol_a}

        self.property_args = property_args

    
    def run(self):
        """ 运行步骤: 提取参数, 运行计算"""

        T = self.temperature   # 温度
        P = self.pressure      # 压力

        # 纯组分参数
        pure_params = self.property_parameters["Pure Components"]
        
        z = deepcopy(self.mole_fraction)   # 摩尔分数, 因为后续要转换聚合物的含量,所以用复制，避免对原变量产生影响
        NC = self.number_of_components     # 组分数量
        NS = self.number_of_segments       # 链段数量
        polymer_index_list = []  # 聚合物索引列表
        polymer_dpn_list = []    # 聚合物DPN列表
        mw_seg_avg_list = []     # 聚合物的链段平均分子质量

        # 计算混合物的理想焓、熵
        H_ig_comp = np.zeros(NC)  # 各组分的理想焓
        S_ig_comp = np.zeros(NC)  # 各组分的理想熵

        for i in range(NC):
            component = self.component_list[i]
            if self.components[component]["type"] == "conventional":
                eqno = pure_params[component]["cp_mol_ig"]["eqno"]
                params = pure_params[component]["cp_mol_ig"]["params"]
                enth_mol_form_ig_ref = pure_params[component]["enth_mol_form_ig_ref"]
                entr_mol_form_ig_ref = pure_params[component]["entr_mol_form_ig_ref"]
                H_ig_comp[i] = enth_ig_comp(T, enth_mol_form_ig_ref, eqno, params, Tref=298.15)
                S_ig_comp[i] = entr_ig_comp(T, P*z[i], entr_mol_form_ig_ref, eqno, params, Tref=298.15, Pref=101325)
            elif self.components[component]["type"] == "polymer":
                if z[i] == 0:
                    H_ig_comp[i] = 0
                    S_ig_comp[i] = 0  
                else:
                    polymer_index_list.append(i)
                    DPN = self.component_attributes[component]["DPN"]
                    polymer_dpn_list.append(DPN)

                    # 聚合物的性质由其组成的链段性质计算
                    sfrac = self.component_attributes[component]["SFRAC"]   # 各链段分数
                    mw_seg = np.zeros(NS)  # 各链段分子质量

                    H_ig_seg = np.zeros(NS)  # 各链段的理想焓
                    S_ig_seg = np.zeros(NS)  # 各链段的理想熵

                    for j in range(NS):
                        segment = self.segment_list[j]
                        mw_seg[j] = pure_params[segment]["mw"]
                        eqno = pure_params[segment]["cp_mol_ig"]["eqno"]
                        params = pure_params[segment]["cp_mol_ig"]["params"]
                        enth_mol_form_ig_ref = pure_params[segment]["enth_mol_form_ig_ref"]
                        entr_mol_form_ig_ref = pure_params[segment]["entr_mol_form_ig_ref"]
                        H_ig_seg[j] = enth_ig_comp(T, enth_mol_form_ig_ref, eqno, params, Tref=298.15)
                        S_ig_seg[j] = entr_ig_comp(T, P*z[i]*sfrac[j], entr_mol_form_ig_ref, eqno, params, Tref=298.15, Pref=101325)

                    mw_seg_avg = sum(sfrac*mw_seg) # 链段的平均相对分子质量
                    mw_seg_avg_list.append(mw_seg_avg)

                    H_ig_comp[i] = np.sum(sfrac*H_ig_seg)*self.mw[i]/mw_seg_avg  # 聚合物的理想焓(以链段表示)
                    S_ig_comp[i] = np.sum(sfrac*S_ig_seg)*self.mw[i]/mw_seg_avg  # 聚合物的理想熵(以链段表示)

        H_ig_mix = np.sum(H_ig_comp*z)   # 混合物的理想焓
        S_ig_mix = np.sum(S_ig_comp*z)   # 混合物的理想熵

        # 计算压力对混合物理想熵的影响
        for i in range(NC):
            # 排除聚合物
            if z[i] != 0 and self.component_list[i] not in self.polymer_list:
                S_ig_mix = S_ig_mix - R*z[i]*np.log(P*z[i]/101325)
        
        # 换算聚合物的含量
        for i in range(len(polymer_index_list)):
            index = polymer_index_list[i]
            DPN = polymer_dpn_list[i]
            mw_seg_avg = mw_seg_avg_list[i]
            z[index] = z[index]*self.mw[index]/mw_seg_avg/DPN
        z = z/sum(z)

        
        # 创建物性方法，并填入相应的物性参数
        self.convert_params_to_args()
        eos = PropertyMethod(self.property_method, self.property_args)

        # 判断相态
        # 如有效相为液相, 汽化分数为0
        if self.valid_phases == "liquid":
            x = z
            y = z
            beta = 0
        # 如有效相为气相, 汽化分数为1
        elif self.valid_phases == "vapor":
            x = z
            y = z
            beta = 1
        # 如有效相为气液, 则要判断
        elif self.valid_phases == "vapor-liquid":
            # 判断是否是纯物质
            index = np.where(z==1)[0]   # 查找组成为1的索引
            # 纯物质, 根据临界温度和饱和蒸气压判断相态
            if len(index) > 0:
                index = index[0]   # 组分的索引
                comp = self.component_list[index]
                Tc = self.property_args["Tc"]
                if self.temperature > Tc[index]:
                    x = z
                    y = z
                    beta = 1
                else:
                    # 计算饱和蒸气压
                    eqno = pure_params[comp]["vapor pressure"]["eqno"]
                    params = pure_params[comp]["vapor pressure"]["params"]
                    vp = vapor_pressure(self.temperature, eqno, params)
                    if self.pressure < vp:
                        x = z
                        y = z
                        beta = 1
                    elif self.pressure > vp:
                        x = z
                        y = z
                        beta = 0
                    else:
                        raise InputError("Please input correct valid phases!")
            # 混合物，通过闪蒸程序判断
            else:            
                beta, x, y, iter = TP_flash(T, P, z, eos, self.tol, self.max_iter)

        # 计算气相的摩尔密度、摩尔剩余焓、摩尔剩余熵
        rho_vap = eos.molar_density(T, P, y, "vap")
        H_res_vap = eos.molar_residual_enthalpy(T, P, y, "vap")
        S_res_vap = eos.molar_residual_entropy(T, P, y, "vap")

        # 计算液相的摩尔密度、摩尔剩余焓、摩尔剩余熵
        rho_liq = eos.molar_density(T, P, x, "liq")
        H_res_liq = eos.molar_residual_enthalpy(T, P, x, "liq")
        S_res_liq = eos.molar_residual_entropy(T, P, x, "liq")

        # 如果存在聚合物, 需要进行换算
        for i in range(len(polymer_index_list)):
            index = polymer_index_list[i]
            DPN = polymer_dpn_list[i]
            mw_seg_avg = mw_seg_avg_list[i]
            x[index] = x[index]*DPN*mw_seg_avg/self.mw[index]

        rho_liq = rho_liq*sum(x)
        H_res_liq = H_res_liq/sum(x)
        S_res_liq = S_res_liq/sum(x)
        beta = beta/(beta+(1-beta)*sum(x)) 

        # 气液混合物的总摩尔焓和熵
        H_mix = H_ig_mix + beta*H_res_vap + (1-beta)*H_res_liq
        S_mix = S_ig_mix + beta*S_res_vap + (1-beta)*S_res_liq  

        self.vapor_mole_fraction = y          # 气相摩尔组成
        self.liquid_mole_fraction = x/sum(x)  # 液相摩尔组成
      
        self.vapor_fraction = beta
        self.liquid_fraction = 1 - self.vapor_fraction
        if beta == 1:
            self.phase = "vapor"
        elif beta == 0:
            self.phase = "liquid"
        else:
            self.phase = "mixed"
        
        self.vapor_volume_flow = self.mole_flow*self.vapor_fraction/rho_vap
        self.liquid_volume_flow = self.mole_flow*self.liquid_fraction/rho_liq
        self.volume_flow = self.mole_flow*(self.vapor_fraction/rho_vap+self.liquid_fraction/rho_liq)
        self.mole_enthalpy = H_mix
        self.enthalpy_flow = self.mole_enthalpy*self.mole_flow
        self.mass_enthalpy = self.enthalpy_flow/self.mass_flow
        self.mole_entropy = S_mix
        self.entropy_flow = self.mole_entropy *self.mole_flow
        self.mass_entropy = self.entropy_flow/self.mass_flow

        self.mole_density = self.mole_flow/self.volume_flow
        self.mass_density = self.mass_flow/self.volume_flow

        self.status = True


    def print_result(self):
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
        column_data = ["", "", self.phase, ""]
        if self.source != None:
            column_data[0] = self.source.name
        if self.destination != None:
            column_data[1] = self.destination.name
        column_data.extend(self.component_mole_flow.tolist())
        column_data.append("")
        column_data.extend(self.mole_fraction.tolist())
        column_data.append("")
        column_data.extend(self.component_mass_flow.tolist())
        column_data.append("")
        column_data.extend(self.mass_fraction.tolist())
        column_data.extend([self.mole_flow, self.mass_flow, self.volume_flow,
            self.temperature, self.pressure, self.vapor_fraction, self.liquid_fraction,
            self.mole_enthalpy, self.mass_enthalpy, self.enthalpy_flow, 
            self.mole_entropy, self.mass_entropy, self.mole_density, 
            self.mass_density, self.avg_mw])
            
        # 添加流股的数据
        tb.add_column(self.name, column_data)

        # 输出表格
        print(tb)
        print("\n",self.name, "Component Attributes: ")
        print(self.component_attributes)


    def print_cld(self):
        """ 以表格形式输出流股中聚合物的链长分布"""

        if self.distributions == {}:
            print(self.name, "has no polymer distribution")
            return

        # 链长分布
        cld = self.distributions[self.polymer_list[0]]["CLD"]
        
        # 创建表格
        tb = pt.PrettyTable()

        # 设置表格的第一列
        tb.add_column("DPN", self.DPN_points)
        ns = len(cld)-1    # 位点数
        # 设置各个位点的数据
        for i in range(1, ns+1):
            tb.add_column("Site_"+str(i), cld[i])
        tb.add_column("Composite", cld[0])

        # 输出表格
        print(self.name," Site-based Chain Length Distribution")
        print(tb)


    def print_mwd(self):
        """ 以表格形式输出流股中聚合物的分子量分布"""

        if self.distributions == {}:
            print(self.name, "has no polymer distribution")
            return

        # 分子量分布
        mwd = self.distributions[self.polymer_list[0]]["MWD"]

        # 创建表格
        tb = pt.PrettyTable()

        # 设置表格的第一列
        tb.add_column("DPN", self.MWN_points)
        ns = len(mwd)-1    # 位点数
        # 设置各个位点的数据
        for i in range(1, ns+1):
            tb.add_column("Site_"+str(i), mwd[i])
        tb.add_column("Composite", mwd[0])

        print(self.name," Site-based Molecular Weight Distribution")
        print(tb)


    def plot_cld(self):
        " Plot Chain Length Distribution"

        if self.distributions == {}:
            print(self.name, "has no polymer distribution")
            return

        window_title = self.name
        figure_title = "Chain Length Distribution"
        xlabel = "log(n)"
        if self.GPC:
            ylabel = "nw(n)"
        else:
            ylabel = "w(n)"
        cld = self.distributions[self.polymer_list[0]]["CLD"]
        distribution_plot(self.DPN_points, cld, window_title, figure_title, xlabel, ylabel)


    def plot_mwd(self):
        " Plot Molecular Weight Distribution"

        if self.distributions == {}:
            print(self.name, "has no polymer distribution")
            return

        window_title = self.name
        figure_title = "Molecular Weight Distribution"
        xlabel = "log(Mn)"
        ylabel = "W(log(Mn))"
        mwd = self.distributions[self.polymer_list[0]]["MWD"]
        distribution_plot(self.MWN_points, mwd, window_title, figure_title, xlabel, ylabel)



# test
if __name__ == "__main__":
    from polymer_model.flowsheet.flowsheet import Flowsheet
    
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
    fs.get_property_parameters()

    s1 = MaterialStream()
    fs.add_stream(s1)

    s2 = MaterialStream()
    fs.add_stream(s2)

    powder = MaterialStream()
    fs.add_stream(powder)
    powder.set_name("Powder")

    powder1 = MaterialStream()
    fs.add_stream(powder1)
    powder1.set_name("Powder1")

    
    z = np.array([0.7695529, 0.1446299, 0.0786014, 0.00721171, 0, 6.26484E-07, 3.46973E-06, 0])
    component_attribute = {
        "Titanium Tetrachloride": {
            "CPSFRAC": 1,
            "CDSFRAC": 0.,
            "CVSFRAC": np.array([0., 0., 0., 0.,]),
            "CISFRAC": np.array([0., 0., 0., 0.,]),

        },
    }
    s1.input(330.546, 3000000, None, 7012.333, z, "Mole & Mole Frac", component_attribute)
    s1.run()
    # s1.print_result()

    z2 = np.array([0.8331688, 0.1668312, 0, 0, 0, 0, 0, 0])
    s2.input(303.15, 3000000, None, 2004.329, z2, "Mole & Mole Frac")
    s2.run()
    # s2.print_result()

    # z3 = np.array([0.0214275, 0.004416, 2.06682E-06, 1.2702E-06, 0.9731196, 0.000158071, 0.000875446, 0])
    z3 = np.array([0.5955161, 0.1227299,5.74411E-05, 3.53016E-05, 27.04503 ,0.00439311, 0.0243304, 0, ])
    component_attribute = {
        "Titanium Tetrachloride": {
            "CPSFRAC": 1,
            "CDSFRAC": 0.,
            "CVSFRAC": np.array([0., 0., 0., 0.,]),
            "CISFRAC": np.array([0., 0., 0., 0.,]),

        },
        "Polypropylene": {
            "SFRAC": np.array([1., ]),
            "DPN": 2024.781, 
        }
    }
    powder.input(338.15, 500000, None, 27.79208, z3, "Mole & Mole Flow", component_attribute, valid_phases="liquid")
    powder.run()


    z4 = np.array([0.2214275, 0.004416, 2.06682E-06, 1.2702E-06, 0.7731196, 0.000158071, 0.000875446, 0])
    component_attribute = {
        "Titanium Tetrachloride": {
            "CPSFRAC": 1,
            "CDSFRAC": 0.,
            "CVSFRAC": np.array([0., 0., 0., 0.,]),
            "CISFRAC": np.array([0., 0., 0., 0.,]),

        },
        "Polypropylene": {
            "SFRAC": np.array([1., ]),
            "DPN": 2024.781, 
        }
    }
    powder1.input(338.15, 3000000, None, 27.79208, z4, "Mole & Mole Frac", component_attribute, max_iter=300)
    powder1.run()

    fs.run()
    fs.print_streams_results()








