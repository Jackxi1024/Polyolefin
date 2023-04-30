"""
Ziegler-Natta催化剂聚合反应动力学 + CSTR反应器模型求解
1. 体系可以是液相, 也可是气液两相(聚合物看做液相)
2. 只考虑均聚 
3. 动力学中不考虑长链支化、副产物
4. 反应器操作条件为指定温度、压力

反应器出口组成通过序贯迭代计算, 以进口组成为出口计算的初始值
并且初始值中必须提供一个位点的各阶矩的初始值, 使得计算收敛

在Aspen中, 链增长为Propagation, 无规链增长为Atact-Prop
链增长包含无规链增长, 无规链增长是链增长的一种, 用于计算无规度(Atact-Prop/Propagation)
"""

# 导包 
from math import exp, log, log10
import math
import numpy as np
from scipy import optimize     # 用于求解CSTR模型(非线性方程组)
from polymer_model.thermodynamics.flash import VLE_TPflash  # 用于计算聚合物体系闪蒸
from polymer_model.utility import SolutionError
import matplotlib.pyplot as plt        # 用于绘制链长分布图和分子量分布图
from polymer_model.property_method.property_method_wrapper import PropertyMethod
from polymer_model.database.pure_components import args as pure_args


R = 8.3145   # 气体常数

def convert_params_to_args(components, params, property_method, DPN, sfrac):
    """
    根据体系包含的组分components, 从物性包params中提取参数, 包装成状态方程所需的参数形式
    """

    # 创建组分列表
    component_list = []  # 小分子和聚合物
    polymer_list = []   # 包含聚合物
    segment_list = []   # 包含链段
    for component in components:
        if components[component]["type"] == "conventional":
            component_list.append(component)
        elif components[component]["type"] == "polymer":
            component_list.append(component)
            polymer_list.append(component)
        elif components[component]["type"] == "segment":
            segment_list.append(component)
    NC = len(component_list)  # 组分数
    NS = len(segment_list)    # 链段数

    # 纯组分参数
    pure_params = params["Pure Components"]
    
    # 如果物性方法为"PC-SAFT", 则需要从物性参数中获取链段数m, 链段直径s, 能量参数e等参数
    if property_method == "PC-SAFT":

        m = np.zeros(NC)      # 各组分链段数, Aspen用PCSFTM表示
        s = np.zeros(NC)      # 各组分链段直径, Aspen用PCSFTV表示
        e = np.zeros(NC)      # 各组分能量参数, Aspen用PCSFTU表示
        e_assoc = np.zeros(NC)  # 关联成分的关联能量,单位K, Aspen用PCSFAU表示
        vol_a = np.zeros(NC)    # 关联成分的有效关联量, Aspen用PCSFAV表示
        w = np.zeros(NC)      # 各组分偏心因子
        Tc = np.zeros(NC)     # 各组分临界温度
        Pc = np.zeros(NC)     # 各组分临界压力
        mw = np.zeros(NC)     # 相对分子质量

        for i in range(NC):
            component = component_list[i]
            
            w[i] = pure_params[component]["w"]
            Tc[i] = pure_params[component]["Tc"]
            Pc[i] = pure_params[component]["Pc"]
            mw[i] = pure_params[component]["mw"]
            if components[component]["type"] == "conventional":
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
            elif components[component]["type"] == "polymer":
                # 聚合物的性质由其组成的链段性质计算
                mw_seg = np.zeros(NS)  # 各链段分子质量
                r_seg = np.zeros(NS)   # 各链段的链段数/数均分子量
                s_seg = np.zeros(NS)   # 各链段的链段直径
                e_seg = np.zeros(NS)   # 各链段的能量参数
                e_assoc_seg = np.zeros(NS)
                vol_a_seg = np.zeros(NS)

                for j in range(NS):
                    segment = segment_list[j]
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
            
        args = {"Tc": Tc, "Pc": Pc, "w": w, "m": m, "s": s, "e": e, 
                "e_assoc": e_assoc, "vol_a": vol_a, "mw":mw}
    
    else:
        raise NotImplementedError("It has not been implemented!")

    return args 


# CSTR反应器, 计算各组分、催化剂各位点和聚合物各阶矩的生成速率
# 返回反应器模型: 进料+生成-出料 = 0
def model(out, components, params, feed, specs, property_method):
    """ 
    Polymerization kinetics of Z-N catalyst
    
    Parameters
    ----------
    out : ndarray

    components : dict

    params : dict

    feed : dict

    specs : dict

        
    Returns
    ----------       
    Reactor model

    """

    # 组分列表
    component_list = []  # 小分子和聚合物
    polymer_list = []   # 包含聚合物
    segment_list = []   # 包含链段
    for component in components:
        if components[component]["type"] == "conventional":
            component_list.append(component)
        elif components[component]["type"] == "polymer":
            component_list.append(component)
            polymer_list.append(component)
        elif components[component]["type"] == "segment":
            segment_list.append(component)

    # 确定各反应物
    species = specs["Reactions"]["Species"]
    
    # 组分摩尔流量
    Fc_in = feed["Component Mole Flow"]
    
    # 催化剂各位点流量
    F_ps_in = feed["CPSFLOW"]
    F_ds_in = feed["CDSFLOW"]
    F_is_in = feed["CISFLOW"]
    F_vs_in = feed["CVSFLOW"]

    # 获取聚合物各阶矩流量
    SFLOW_in = feed["SFLOW"]
    LSZMOM_in = feed["LSZMOM"]
    DSZMOM_in = feed["DSZMOM"]
    LSFMOM_in = feed["LSFMOM"]
    DSFMOM_in = feed["DSFMOM"]
    LSSMOM_in = feed["LSSMOM"]
    DSSMOM_in = feed["DSSMOM"]

    NC = len(Fc_in)      # 体系组分数
    NS = len(SFLOW_in)   # 聚合物链段数
    ns = len(F_vs_in)    # 催化剂位点数


    # 反应器参数
    T = specs["Temperature"]   # 温度
    P = specs["Pressure"]   # 压力

    # 确定反应器体积和各相体积
    # 1、如果CSTR内含有气液两相
    if specs["Valid phases"] == "vapor-liquid":
        if specs["Specification type"] == "Reactor volume & Phase volume":
            Vr = specs["Reactor volume"]   # 反应器体积
            if "Liquid volume" in specs :
                Vl = specs["Liquid volume"]   # 液相体积
                Vv = Vr - Vl                  # 气相体积
            elif "Vapor volume" in specs :
                Vv = specs["Vapor volume"]    # 气相体积  
                Vl = Vr - Vv                  # 液相体积
            # 确定物料的相态和反应相的体积
            reaction_phase = specs["Reactions"]["Reacting phase"]  # 反应的相
            if reaction_phase.lower() in ("l", "liq", "liquid"):
                V = Vl
            else :
                V = Vr
    # 2、如果CSTR内只有液相
    elif specs["Valid phases"] in ["vapor", "liquid"]:
        V = specs["Reactor volume"]   # 反应器体积,也是反应体积

    # 获取出口信息
    Fc_out = out[0: NC]       # 各组分出口摩尔流量
    F_out = np.sum(Fc_out)    # 出口总摩尔流量(混合物)
    z = Fc_out/F_out          # 各组分出口摩尔分数

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

    SFLOW = np.zeros(NS)   # 链段流率
    for mono in species["monomers"]:
        seg = species["segments"][mono]
        mono_index = component_list.index(mono)
        seg_index = segment_list.index(seg)
        # 生成的链段流率 == 对应的单体反应流率
        SFLOW[seg_index] = SFLOW_in[seg_index] + (Fc_in[mono_index]-Fc_out[mono_index])

    if (np.sum(SFLOW) == 0):
        sfrac = np.zeros(NS)
        sfrac.fill(1/NS)
    else:
        sfrac = SFLOW/np.sum(SFLOW)   # 链段分数
    

    # 计算聚合物性质
    DPN = np.sum(LSFMOM+DSFMOM)/np.sum(LSZMOM+DSZMOM)  # 数均聚合度

    # 创建物性方法，并填入相应的物性参数
    property_args = convert_params_to_args(components, params, property_method, DPN, sfrac)
    eos = PropertyMethod(property_method, property_args)
    
    mw = property_args["mw"]

    # 根据DPN和MWN, 换算聚合物的摩尔分数和链段数
    polymer_index = component_list.index(species["polymer"])
    monomers_index = component_list.index(species["monomers"][0])
    z1 = z.copy()          # 复制组分数组, 不可直接修改z
    z1[polymer_index] = z1[polymer_index]/DPN
    z1 = z1/np.sum(z1)

    # 判断是否处于气液两相, 计算各组分液相浓度
    if specs["Valid phases"] == "vapor-liquid":
        # 执行闪蒸程序
        tol = specs["Flash params"]["Tol"]
        max_iter = specs["Flash params"]["Max iter"]
        beta, x, y, iter_total = VLE_TPflash(T, P, z1, eos, tol=tol, max_iter=max_iter) 
        # 计算各组分液相浓度
        rho_l = eos.molar_density(T, P, x, "liq")   # 液相摩尔密度(mol/m3)
        c = rho_l*x     # 各组分液相浓度
        
        x[polymer_index] = x[polymer_index]*DPN
        beta = beta/(beta+(1-beta)*sum(x))
        x = x/sum(x)
        F0 = F_out*(1-beta)*x[monomers_index]  # 单体液相流量
        c0 = c[monomers_index]   # 液相单体浓度
    elif specs["Valid phases"] in ["liquid", "vapor"]:
        rho_l = eos.molar_density(T, P, z1, "liq")   # 液相摩尔密度(mol/m3)
        c = rho_l*z1     # 各组分浓度
        F0 = F_out*z[monomers_index]
        c0 = c[monomers_index]

    # 计算各催化位点的浓度
    C_ps = F_ps/F0*c0
    C_ds = F_ds/F0*c0
    C_is = F_is/F0*c0 
    C_vs = F_vs/F0*c0

    # 计算聚合物浓度的各阶矩
    Y0 = LSZMOM/F0*c0   # 活性链零阶矩
    X0 = DSZMOM/F0*c0   # 死聚物零阶矩
    Y1 = LSFMOM/F0*c0   # 活性链一阶矩
    X1 = DSFMOM/F0*c0   # 死聚物一阶矩
    Y2 = LSSMOM/F0*c0   # 活性链二阶矩
    X2 = DSSMOM/F0*c0   # 死聚物二阶矩

    # 各种净生成速率(mol/(m3*s))
    r = np.zeros(NC)      # 各个组分的净生成速率
    r_ps = 0              # 潜在活性位的净生成速率
    r_ds = 0              # 失活位点的净生成速率
    r_is = np.zeros(ns)   # 抑制位点的净生成速率(对应不同的催化剂位点)
    r_vs = np.zeros(ns)   # 空活性位点的净生成速率(对应不同的催化剂位点)
    r_Y0 = np.zeros(ns)   # 活性链零阶矩的净生成速率(对应不同的催化剂位点)
    r_X0 = np.zeros(ns)   # 死聚物零阶矩的净生成速率(对应不同的催化剂位点)
    r_Y1 = np.zeros(ns)   # 活性链一阶矩的净生成速率(对应不同的催化剂位点)
    r_X1 = np.zeros(ns)   # 死聚物一阶矩的净生成速率(对应不同的催化剂位点)
    r_Y2 = np.zeros(ns)   # 活性链二阶矩的净生成速率(对应不同的催化剂位点)
    r_X2 = np.zeros(ns)   # 死聚物二阶矩的净生成速率(对应不同的催化剂位点)

    K_td = np.zeros(ns)   # 链转移速率+链失活速率
    K_p = np.zeros(ns)    # 链增长速率 
          
    # 各种反应的反应速率
    for reaction in specs["Reactions"]["Reactions"]:
        type = reaction[0]                               # ZN聚合反应类型
        site = reaction[1] - 1                           # 反应位点(python从0开始)
        comp1 = reaction[2]                              # 组分1
        if comp1 != None and comp1 in component_list:    
            comp1_index = component_list.index(comp1)    # 组分1的索引
        comp2 = reaction[3]                              # 组分2
        if comp2 != None and comp2 in component_list:
            comp2_index = component_list.index(comp2)    # 组分2的索引
        k0 = reaction[4]          # 指前因子
        Ea = reaction[5]          # 活化能
        order = reaction[6]       # 反应物的反应级数
        tdb_frac = reaction[7]    # 终端双键分数
        T_ref = reaction[8]       # 参考温度

        # 反应速率常数
        k = k0*exp(-Ea/R*(1/T-1/T_ref))
        
        # Site Activation
        if type in ["Act-Cocat", "Act-Edonor", "Act-H2"] :
            r_a = k*C_ps*c[comp2_index]**order
            r[comp2_index] = r[comp2_index] - r_a  # 助催化剂/电子供体/氢的净生成速率
            r_ps = r_ps - r_a                      # 潜在位点的净生成速率
            r_vs[site] = r_vs[site] + r_a          # 空位点的净生成速率 
            # 聚合物的净生成速率
            r[polymer_index] = r[polymer_index] + r_a*mw[comp2_index]/mw[polymer_index]
        elif type == "Act-Mon":                    # 单体活化反应中单体净生成为0
            r_a = k*C_ps*c[comp2_index]**order
            r_ps = r_ps - r_a
            r_vs[site] = r_vs[site] + r_a
        elif type == "Act-Spon":
            r_a = k*C_ps**order 
            r_ps = r_ps - r_a
            r_vs[site] = r_vs[site] + r_a

        # Chain Initiation
        elif type == "Chain-Ini":
            r_i = k*C_vs[site]*c[comp1_index]**order
            r[comp1_index] = r[comp1_index] - r_i   # 单体净生成速率
            r_vs[site] = r_vs[site] - r_i           # 空活性位点净生成速率
            # 链引发生成链长为1的活性链，对活性链的矩有影响
            r_Y0[site] = r_Y0[site] + r_i   
            r_Y1[site] = r_Y1[site] + r_i
            r_Y2[site] = r_Y2[site] + r_i
            r[polymer_index] = r[polymer_index] + r_i  # 链段流量的净生成速率

        # Chain Propagation
        elif type == "Propagation":
            k_p = k*c[comp2_index]**order            # 单位活性链浓度的链增长速率
            r_p = k_p*Y0[site]                       # 链增长速率 
            r[comp2_index] = r[comp2_index] - r_p    # 单体净生成速率
            r_Y0[site] = r_Y0[site] + 0              # 链增长对零阶矩无影响
            r_Y1[site] = r_Y1[site] + r_p            # 链增长对一阶矩的影响
            r_Y2[site] = r_Y2[site] + k_p*(2*Y1[site]+Y0[site])  # 链增长对二阶矩的影响
            r[polymer_index] = r[polymer_index] + r_p  # 链段流量的净生成速率
            K_p[site] = K_p[site]+k_p
        # 无规链增长, 之后要用于计算无规度(暂不处理)
        elif type == "Atact-Prop":
            pass
            # k_p = k*c[comp2_index]**order            # 单位活性链浓度的链增长速率
            # r_p = k_p*Y0[site]                       # 链增长速率 
            # r[comp2_index] = r[comp2_index] - r_p    # 单体净生成速率
            # r_Y0[site] = r_Y0[site] + 0              # 链增长对零阶矩无影响
            # r_Y1[site] = r_Y1[site] + r_p            # 链增长对一阶矩的影响
            # r_Y2[site] = r_Y2[site] + k_p*(2*Y1[site]+Y0[site])  # 链增长对二阶矩的影响
            # r[polymer_index] = r[polymer_index] + r_p  # 链段流量的净生成速率
            # K_p[site] = K_p[site]+k_p

        # Chain Transfer
        elif type in ["Chat-H2", "Chat-Cocat", "Chat-Sol", "Chat-Agent", "Chat-Edonor"]:
            k_t = k*c[comp2_index]**order            # 单位活性链浓度的链转移速率
            r_t = k_t*Y0[site]                       # 链转移速率
            r[comp2_index] = r[comp2_index] - r_t    # 氢/助催化剂/溶剂/转移剂/电子供体净生成速率
            r_vs[site] = r_vs[site] + r_t            # 空位点的净生成速率
            r_Y0[site] = r_Y0[site] - k_t*Y0[site]   # 链转移对活性链零阶矩的影响
            r_Y1[site] = r_Y1[site] - k_t*Y1[site]   # 链转移对活性链一阶矩的影响
            r_Y2[site] = r_Y2[site] - k_t*Y2[site]   # 链转移对活性链二阶矩的影响
            r_X0[site] = r_X0[site] + k_t*Y0[site]   # 链转移对死聚物零阶矩的影响
            r_X1[site] = r_X1[site] + k_t*Y1[site]   # 链转移对死聚物一阶矩的影响
            r_X2[site] = r_X2[site] + k_t*Y2[site]   # 链转移对死聚物二阶矩的影响
            # 聚合物的净生成速率
            r[polymer_index] = r[polymer_index] + r_t*mw[comp2_index]/mw[polymer_index]
            K_td[site] = K_td[site] + k_t
        elif type == "Chat-Spon":
            k_t = k
            r_t = k_t*Y0[site]
            r_vs[site] = r_vs[site] + r_t            # 空位点的净生成速率
            r_Y0[site] = r_Y0[site] - k_t*Y0[site]   # 链转移对活性链零阶矩的影响
            r_Y1[site] = r_Y1[site] - k_t*Y1[site]   # 链转移对活性链一阶矩的影响
            r_Y2[site] = r_Y2[site] - k_t*Y2[site]   # 链转移对活性链二阶矩的影响
            r_X0[site] = r_X0[site] + k_t*Y0[site]   # 链转移对死聚物零阶矩的影响
            r_X1[site] = r_X1[site] + k_t*Y1[site]   # 链转移对死聚物一阶矩的影响
            r_X2[site] = r_X2[site] + k_t*Y2[site]   # 链转移对死聚物二阶矩的影响 
            K_td[site] = K_td[site] + k_t
        elif type == "Chat-Mon":
            k_t = k*c[comp2_index]**order            # 单位活性链浓度的链转移速率
            r_t = k_t*Y0[site]                       # 链转移速率
            r[comp2_index] = r[comp2_index] - r_t    # 单体净生成速率
            r_Y0[site] = r_Y0[site] - 0              # 链向单体转移对活性链零阶矩无影响
            r_Y1[site] = r_Y1[site] + r_t - k_t*Y1[site]   # 链转移对活性链一阶矩的影响
            r_Y2[site] = r_Y2[site] + r_t - k_t*Y2[site]   # 链转移对活性链二阶矩的影响
            r_X0[site] = r_X0[site] + k_t*Y0[site]         # 链转移对死聚物零阶矩的影响
            r_X1[site] = r_X1[site] + k_t*Y1[site]         # 链转移对死聚物一阶矩的影响
            r_X2[site] = r_X2[site] + k_t*Y2[site]         # 链转移对死聚物二阶矩的影响
            r[polymer_index] = r[polymer_index] + r_t      # 链段流量的净生成速率
            K_td[site] = K_td[site] + k_t

        # Site Inhibition
        elif type == "Fsinh-H2" or type == "Fsinh-Poison":
            r_if = k*C_vs[site]*c[comp1_index]**order
            r[comp1_index] = r[comp1_index] - r_if    # 氢净生成速率
            r_vs[site] = r_vs[site] - r_if            # 空位点的净生成速率
            r_is[site] = r_is[site] + r_if            # 抑制位点的净生成速率
            # 聚合物的净生成速率
            r[polymer_index] = r[polymer_index] + r_if*mw[comp1_index]/mw[polymer_index]
        elif type == "Rsinh-H2" or type == "Rsinh-Poison":
            r_ir = k*C_is[site]**order
            r[comp1_index] = r[comp1_index] + r_ir    # 氢净生成速率
            r_vs[site] = r_vs[site] + r_ir            # 空位点的净生成速率
            r_is[site] = r_is[site] - r_ir            # 抑制位点的净生成速率
            # 聚合物的净生成速率
            r[polymer_index] = r[polymer_index] - r_ir*mw[comp1_index]/mw[polymer_index]

        # 长链支化暂不考虑    
        elif type == "Tdb-Poly":
            r_LCB = k*Y0[site]*c[comp1_index]
            
        # Cocatalyst Poison:暂未考虑
        elif type == "Cocat-Poison":
            r_AX = k*c[comp1_index]*c[comp2_index]

        # Site Deactivation
        elif type in ["Deact-Cocat", "Deact-Edonor", "Deact-H2", "Deact-Poison", "Deact-Mon"]:
            k_d = k*c[comp1_index]**order              # 单位活性位点的失活速率
            r[comp1_index] = r[comp1_index] - k_d*(C_vs[site] + Y0[site])   # 助催化剂/电子供体/氢/毒物/单体净生成速率
            r_ds = r_ds + k_d*(C_vs[site] + Y0[site])  # 死位点的净生成速率
            r_vs[site] = r_vs[site] - k_d*C_vs[site]   # 空位点的净生成速率
            r_Y0[site] = r_Y0[site] - k_d*Y0[site]     # 位点失活对活性链零阶矩的影响
            r_Y1[site] = r_Y1[site] - k_d*Y1[site]     # 位点失活对活性链一阶矩的影响
            r_Y2[site] = r_Y2[site] - k_d*Y2[site]     # 位点失活对活性链二阶矩的影响
            r_X0[site] = r_X0[site] + k_d*Y0[site]     # 位点失活对死聚物零阶矩的影响
            r_X1[site] = r_X1[site] + k_d*Y1[site]     # 位点失活对死聚物一阶矩的影响
            r_X2[site] = r_X2[site] + k_d*Y2[site]     # 位点失活对死聚物二阶矩的影响 
            # 聚合物的净生成速率
            r[polymer_index] = r[polymer_index] + k_d*(C_vs[site] + Y0[site])*mw[comp1_index]/mw[polymer_index]
            K_td[site] = K_td[site] + k_d
        elif type == "Deact-Spon":
            k_d = k
            r_ds = r_ds + k_d*(C_vs[site] + Y0[site])  # 死位点的净生成速率
            r_vs[site] = r_vs[site] - k_d*C_vs[site]   # 空位点的净生成速率
            r_Y0[site] = r_Y0[site] - k_d*Y0[site]     # 位点失活对活性链零阶矩的影响
            r_Y1[site] = r_Y1[site] - k_d*Y1[site]     # 位点失活对活性链一阶矩的影响
            r_Y2[site] = r_Y2[site] - k_d*Y2[site]     # 位点失活对活性链二阶矩的影响
            r_X0[site] = r_X0[site] + k_d*Y0[site]     # 位点失活对死聚物零阶矩的影响
            r_X1[site] = r_X1[site] + k_d*Y1[site]     # 位点失活对死聚物一阶矩的影响
            r_X2[site] = r_X2[site] + k_d*Y2[site]     # 位点失活对死聚物二阶矩的影响 
            K_td[site] = K_td[site] + k_d

    # (链转移速率+链失活速率)/链增长速率, 用于计算链长分布
    # K_td和K_p可能为0，因此将τ中的Nan替换为0
    τ = np.zeros(len(K_td))
    for i in range(len(K_td)):
        if K_p[i] != 0:              
            τ[i] = K_td[i]/K_p[i]

    # 物料衡算方程(以摩尔流量为基准): 输入 + 净生成 = 输出
    equations = np.zeros(len(out))
    for i in range(NC):
        equations[i] = Fc_in[i] + r[i]*V - Fc_out[i]      

    # 催化剂各种位点流量衡算
    equations[NC] = F_ps_in + r_ps*V - F_ps             
    equations[1+NC] = F_ds_in + r_ds*V - F_ds
    equations[2+NC: 2+NC+ns] = F_is_in + r_is*V - F_is
    equations[2+NC+ns: 2+NC+2*ns] = F_vs_in + r_vs*V - F_vs
    # 流量矩衡算
    equations[2+NC+2*ns: 2+NC+3*ns] = LSZMOM_in + r_Y0*V - LSZMOM    # 活性链零阶矩
    equations[2+NC+3*ns: 2+NC+4*ns] = DSZMOM_in + r_X0*V - DSZMOM  # 死聚物零阶矩
    equations[2+NC+4*ns: 2+NC+5*ns] = LSFMOM_in + r_Y1*V - LSFMOM  # 活性链一阶矩
    equations[2+NC+5*ns: 2+NC+6*ns] = DSFMOM_in + r_X1*V - DSFMOM  # 死聚物一阶矩
    equations[2+NC+6*ns: 2+NC+7*ns] = LSSMOM_in + r_Y2*V - LSSMOM  # 活性链二阶矩
    equations[2+NC+7*ns: 2+NC+8*ns] = DSSMOM_in + r_X2*V - DSSMOM  # 死聚物二阶矩

    return equations, τ, SFLOW


def equations(out, components, params, feed, specs, method):
    return model(out, components, params, feed, specs, method)[0]


def cstr_run(components, params, feed, specs, method):
    """ 
    Solve the CSTR model and give the results
    
    Parameters
    ----------
    components : dict

    params : dict

    polymers : dict

    streams : dict

    specs : dict

        
    Returns
    ----------       
    Reactor model

    """

    Fc_in = feed["Component Mole Flow"]         # 各组分的摩尔流量, mol/s
    NC = len(Fc_in)

    # 催化剂各位点流量
    F_ps_in = feed["CPSFLOW"]
    F_ds_in = feed["CDSFLOW"]
    F_is_in = feed["CISFLOW"]
    F_vs_in = feed["CVSFLOW"]
    ns = len(F_vs_in)

    # 链段流率
    SFLOW_in = feed["SFLOW"]   

    # 获取聚合物各阶矩流量
    LSZMOM_in = feed["LSZMOM"]
    DSZMOM_in = feed["DSZMOM"]
    LSFMOM_in = feed["LSFMOM"]
    DSFMOM_in = feed["DSFMOM"]
    LSSMOM_in = feed["LSSMOM"]
    DSSMOM_in = feed["DSSMOM"]


    # 设置出口的初始值(以进口为初始值)
    x0 = np.zeros(2+NC+8*ns)
    x0[0:NC] = Fc_in
    x0[NC] = 1e-6
    x0[1+NC] = F_ds_in
    x0[2+NC: 2+NC+ns] = F_is_in
    x0[2+NC+ns: 2+NC+2*ns] = F_vs_in

    # 给定零阶矩的初值
    if np.all(LSZMOM_in == 0):
        x0[2+NC+2*ns] = 1e-5
    else:
        x0[2+NC+2*ns: 2+NC+3*ns] = LSZMOM_in

    if np.all(DSZMOM_in == 0):
        x0[2+NC+3*ns] = 0.01
    else:
        x0[2+NC+3*ns: 2+NC+4*ns] = DSZMOM_in
    
    # 给定一阶矩的初值    
    if np.all(LSFMOM_in == 0):
        x0[2+NC+4*ns] = 0.1
    else:
        x0[2+NC+4*ns: 2+NC+5*ns] = LSFMOM_in
    
    if np.all(DSFMOM_in == 0):
        x0[2+NC+5*ns] = 10
    else:
        x0[2+NC+5*ns: 2+NC+6*ns] = DSFMOM_in 
    
    # 给定二阶矩的初值, 可加快收敛
    if np.all(LSSMOM_in == 0):
        x0[2+NC+6*ns] = 1000
    else:
        x0[2+NC+6*ns: 2+NC+7*ns] = LSSMOM_in

    if np.all(DSSMOM_in == 0):
        x0[2+NC+7*ns] = 100000  
    else:
        x0[2+NC+7*ns: 2+NC+8*ns] = DSSMOM_in

    sol = optimize.root(equations, x0, args=(components,params, feed, specs, method), method='hybr', tol=1e-8)
    if sol.success :
        print("status: ", sol.status)
        print("fun: ", sol.fun)
        print("nfev: ", sol.nfev)
        return sol.x
    else : 
        raise SolutionError(sol.message)
    

def parse_result(NC:float, ns:float, out:np.ndarray):
    Fc_out = out[0:NC]       # 出口各组分摩尔流量
    F_out = np.sum(Fc_out)   # 出口总摩尔流量
    z = Fc_out/F_out         # 出口各组分摩尔分数
    z[z<1e-17] = 0      # 将组分含量极少的置为0

    # 设置产品
    product = {}
    product["Component Mole Flow"] = Fc_out    # 各组分的摩尔流量, mol/s
    product["Mole Fraction"] = z
    
    # 设置催化剂
    product["CPSFLOW"] = out[NC]
    product["CDSFLOW"] = out[1+NC]
    product["CISFLOW"] = out[2+NC: 2+NC+ns]
    product["CVSFLOW"] = out[2+NC+ns: 2+NC+2*ns]

    # 设置聚合物
    product["LSZMOM"] = out[2+NC+2*ns: 2+NC+3*ns]
    product["DSZMOM"] = out[2+NC+3*ns: 2+NC+4*ns]
    product["LSFMOM"] = out[2+NC+4*ns: 2+NC+5*ns]
    product["DSFMOM"] = out[2+NC+5*ns: 2+NC+6*ns]
    product["LSSMOM"] = out[2+NC+6*ns: 2+NC+7*ns]
    product["DSSMOM"] = out[2+NC+7*ns: 2+NC+8*ns]

    product["ZMOM"] = np.sum(product["LSZMOM"]) + np.sum(product["DSZMOM"]) 
    product["FMOM"] = np.sum(product["LSFMOM"]) + np.sum(product["DSFMOM"]) 
    product["SMOM"] = np.sum(product["LSSMOM"]) + np.sum(product["DSSMOM"]) 

    product["DPN"] = product["FMOM"]/product["ZMOM"]
    product["DPW"] = product["SMOM"]/product["FMOM"]
    return product


def get_points(Np=100, upper=100000):
    """
    生成聚合物分布的离散链长数组
    使用以下公式，在1和指定上限之间以对数步长间隔分布链长:
    ni = max(i, pow(10, i*math.log10(upper)/N))

    其中i在1和Np之间变化，该间距在整个分子量谱上提供了良好的分辨率，
    重点放在离散过程中容易丢失的低分子量物种上。
    为确保准确性，upper应至少设置为比观察到的重量平均聚合度高五倍.

    Paramater
    Np: float, number of points(default 100)
    upper: Upper limit of chain length

    Returns 
    n: ndarray, shape(Np,),  [1, 2, 3, 4, , ,89125, 100000]
    """

    n = np.zeros(Np)
    for i in range(1, Np+1):
        ri = math.floor(math.pow(10, i*math.log10(upper)/Np)) 
        if (i > ri):
            n[i-1] = i
        else :
            n[i-1] = ri
    return n

# 用于和aspen所计算的链长分布数据做比较
def cld(τ, SFMOM, n, GPC=True):
    """
    聚合物链长分布: w[n] = mf[i]*n*τ**2*(1/(1+τ))**(n+1)
    对数形式: w[log(n)] = w[n]*n*ln10 = mf[i]*ln10*n**2*τ**2*(1/(1+τ))**(n+1)

    Paramaters
    τ : ndarray, shape (ns,)
        表示各个活性位点上(链转移速率+链失活速率)/链增长速率, ns为位点个数
    n : ndarray, shape (N,)
        表示离散的链长, N是在链长分布图中的最大链长
    GPC : True代表执行GPC计算, False代表不执行GPC计算
          GPC计算的纵坐标为n*w(n)
    
    Returns
    w : ndarray, shape (ns+1, N)
        聚合物链长分布的二维数据
    """

    ns = len(τ)       # 活性位点数
    mf = np.zeros(ns)  
    w = np.zeros((ns+1, len(n)))
    for i in range(ns):
        mf[i] = SFMOM[i]/np.sum(SFMOM)    # 某个活性位点生成的聚合物占总聚合物的质量分数
        if GPC:
            w[i+1] = mf[i]*n**2*τ[i]**2*(1/(1+τ[i]))**(n+1)  # 各个活性位点上，每个链长对应的质量分数
        else :
            w[i+1] = mf[i]*n*τ[i]**2*(1/(1+τ[i]))**(n+1)  # 各个活性位点上，每个链长对应的质量分数
        w[0] = w[0] + w[i+1]   # 每个链长对应的质量分数(所有活性位点相加)
    return w

def cld_cum(τ1, SFMOM1, τ2, SFMOM2, n, GPC=True):
    """
    累计的链长分布: w[n] = (m1*w1[n] + m2*w2[n])/(m1+m2)
    """

    w1 = cld(τ1, SFMOM1, n, GPC)   # 进料各位点上链长为n的聚合物占总聚合物的质量分数 
    w2 = cld(τ2, SFMOM2, n, GPC)   # 反应产生的各位点上链长为n的聚合物占反应产生的总聚合物的质量分数
    m1 = sum(SFMOM1)
    m2 = sum(SFMOM2)
    w = (m1*w1 + m2*w2)/(m1+m2)
    return w

def mwd(τ, SFMOM, mw, Mn):
    """  
    聚合物分子量分布: w[Mn] = mf[i]*Mn*(τ/mw)**2*(1/(1+τ))**(Mn/mw+1)
    取对数形式: w[logMn] = mf[i]*ln10*Mn**2*(τ/mw)**2*(1/(1+τ))**(Mn/mw+1)

    Paramaters
    τ : ndarray, shape (ns,)
        表示各个活性位点上(链转移速率+链失活速率)/链增长速率, ns为位点个数
    SFMOM : ndarray, shape (ns,)
        各位点(活性链+死聚物)流量一阶矩, 用于计算各个位点所产生的聚合物的质量分数
    mw : float, 链段的平均分子量
    Mn : ndarray, shape (N,)
        表示离散的聚合物分子量, N是在分子量分布图中的最大分子量
    
    Returns
    w : ndarray, shape (ns+1, N)
        聚合物分子量分布的二维数据
    """
    ns = len(τ)                # 活性位点数
    mf = np.zeros(ns)  
    w = np.zeros((ns+1, len(Mn)))
    for i in range(ns):
        mf[i] = SFMOM[i]/np.sum(SFMOM)      # 某个活性位点生成的聚合物占总聚合物的质量分数
        w[i+1] = mf[i]*log(10)*Mn**2*(τ[i]/mw)**2*(1/(1+τ[i]))**(Mn/mw+1)  # 各个活性位点上，每个链长对应的质量分数
        w[0] = w[0] + w[i+1]   # 每个链长对应的质量分数(所有活性位点相加)
    return w 

def mwd_cum(τ1, SFMOM1, τ2, SFMOM2, mw, Mn):
    """
    累计的分子量分布: w[n] = (m1*w1[n] + m2*w2[n])/(m1+m2)
    """

    w1 = mwd(τ1, SFMOM1, mw, Mn)   # 进料各位点上分子量为Mn的聚合物占总聚合物的质量分数 
    w2 = mwd(τ2, SFMOM2, mw, Mn)   # 反应产生的各位点上链长为Mn的聚合物占反应产生的总聚合物的质量分数
    m1 = sum(SFMOM1)
    m2 = sum(SFMOM2)
    w = (m1*w1 + m2*w2)/(m1+m2)
    return w


def distribution_plot(x, y, window_title="", figure_title="", xlabel="", ylabel="", fontsize=16):
    """ 绘制聚合物分布图(对数横坐标)"""

    plt.figure(window_title)               # 窗口标题
    plt.title(figure_title, fontsize=fontsize)   # 图标题
    plt.xlabel(xlabel, fontsize=fontsize)   # 横坐标标签
    plt.ylabel(ylabel, fontsize=fontsize)  # 纵坐标标签
    plt.xticks(fontsize=fontsize)                               # 设置横坐标刻度字体
    plt.yticks(fontsize=fontsize)                               # 设置纵坐标刻度字体
    plt.xlim(0, np.log10(np.max(x)))                      # 设置横坐标刻度上下限 
    plt.ylim(0, np.max(y)*1.1)                            # 设置纵坐标刻度下限 
    # python颜色对照表: https://blog.csdn.net/qq_42612717/article/details/106291819
    color = ['b','g','r','c','m','y','k','w','Gray','Navy']   # 颜色列表
    linestyle=['-','--','-.',':']                       # 线型
    marker=['o','s','d','v','^','<','>','p','*','x']    # 点型
    for i in range(len(y)):
        if i == 0:
            # 绘制链长分布的复合曲线
            plt.plot(np.log10(x), y[i], label="Composite", linestyle="-", color=color[i], marker=marker[i], linewidth=2, markersize=6)    
        else :
            # 绘制各位点的链长分布曲线
            plt.plot(np.log10(x), y[i], label="Site"+str(i), linestyle="-", color=color[i], marker=marker[i],  linewidth=2, markersize=6) 
    plt.legend()
    # plt.grid()   # 网格线
    plt.show()


def cld_plot(n, w):
    """ 绘制聚合物链长分布图(对数横坐标) """

    plt.figure("Chain Length Distribution")               # 窗口标题
    plt.title("Chain Length Distribution", fontsize=16)   # 图标题
    plt.xlabel("log(n)", fontsize=16)   # 横坐标标签
    plt.ylabel('nw(n) for GPC \n w(n) for not GPC', fontsize=16)  # 纵坐标标签
    plt.xticks(fontsize=16)                               # 设置横坐标刻度字体
    plt.yticks(fontsize=16)                               # 设置纵坐标刻度字体
    plt.xlim(0, np.log10(np.max(n)))                      # 设置横坐标刻度上下限 
    plt.ylim(0, np.max(w)*1.1)                            # 设置纵坐标刻度下限 
    # python颜色对照表: https://blog.csdn.net/qq_42612717/article/details/106291819
    color = ['b','g','r','c','m','y','k','w','Gray','Navy']   # 颜色列表
    linestyle=['-','--','-.',':']                       # 线型
    marker=['o','s','d','v','^','<','>','p','*','x']    # 点型
    for i in range(len(w)):
        if i == 0:
            # 绘制链长分布的复合曲线
            plt.plot(np.log10(n), w[i], label="Composite", linestyle="-", color=color[i], marker=marker[i], linewidth=2, markersize=6)    
        else :
            # 绘制各位点的链长分布曲线
            plt.plot(np.log10(n), w[i], label="Site"+str(i), linestyle="-", color=color[i], marker=marker[i],  linewidth=2, markersize=6) 
    plt.legend()
    # plt.grid()   # 网格线
    plt.show()

def mwd_plot(MW, w):
    """ 绘制聚合物链长分布图(对数横坐标) """

    plt.figure("Molecular Weight Distribution")               # 窗口标题
    plt.title("Molecular Weight Distribution", fontsize=16)   # 图标题
    plt.xlabel("log(Mn)", fontsize=16)                    # 横坐标标签
    plt.ylabel("W(log(Mn))", fontsize=16)            # 纵坐标标签
    plt.xticks(fontsize=16)                               # 设置横坐标刻度字体
    plt.yticks(fontsize=16)                               # 设置纵坐标刻度字体
    plt.xlim(0, np.log10(np.max(MW)))                      # 设置横坐标刻度上下限 
    plt.ylim(0, np.max(w)*1.1)                            # 设置纵坐标刻度下限 
    # python颜色对照表: https://blog.csdn.net/qq_42612717/article/details/106291819
    color = ['b','g','r','c','m','y','k','w','Gray','Navy']   # 颜色列表
    linestyle=['-','--','-.',':']                       # 线型
    marker=['o','s','d','v','^','<','>','p','*','x']    # 点型
    for i in range(len(w)):
        if i == 0:
            # 绘制链长分布的复合曲线
            plt.plot(np.log10(MW), w[i], label="Composite", linestyle="-", color=color[i], marker=marker[i], linewidth=2, markersize=6)    
        else :
            # 绘制各位点的链长分布曲线
            plt.plot(np.log10(MW), w[i], label="Site"+str(i), linestyle="-", color=color[i], marker=marker[i],  linewidth=2, markersize=6) 
    plt.legend()
    plt.show()



if __name__ == "__main__":

    # 反应器1的计算
    # 设置组分
    components = {
        "Titanium Tetrachloride": {"type": "conventional"},
        "Triethyl Aluminum": {"type": "conventional"},
        "Ethylene": {"type": "conventional"},
        "Hydrogen": {"type": "conventional"},
        "N-hexane": {"type": "conventional"},
        "HDPE": {"type": "polymer"},
        "Ethylene-R": {"type": "segment"},
    }

    # 设置聚合物链段
    segments = {"Ethylene-R":{"type": "repeat"}}

    # 设置催化剂
    catalyst = {"Titanium Tetrachloride": {"type": "Z-N", 
                        "site_types_num":4,   # 位点数
                        "site_conc":0.1,      # 位点浓度, mol/kgcat
                }}

    ploymers = {
        "segments": segments,
        "catalyst": catalyst,
    }

    feed = {
        "Component Mole Flow": np.array([0.00878623, 0.0218978 ,59.40975, 0.0413384, 174.0107, 0]),
        "CPSFLOW": 0.000166667,
        "CDSFLOW": 0,
        "CVSFLOW": np.zeros(4),
        "CISFLOW": np.zeros(4), 
        "SFLOW": np.zeros(1),
        "LSZMOM": np.zeros(4),
        "DSZMOM": np.zeros(4),
        "LSFMOM": np.zeros(4),
        "DSFMOM": np.zeros(4),
        "LSSMOM": np.zeros(4),
        "DSSMOM": np.zeros(4),
    }

    # 设置Z-N反应物种
    species = { "polymer": "HDPE",
                "tdb segment": None,
                "monomers": ["Ethylene"],
                "segments": {"Ethylene": "Ethylene-R"},
                "precatalyst": None,
                "catalyst": ["Titanium Tetrachloride"],
                "cocatalysts": ["Triethyl Aluminum"],
                "solvents": ["N-hexane"],
                "transfer agent": None,
                "hydrogens": ["Hydrogen"],
                "poisons": None,
                "elec don": None,
                "byproduct": None,
            }


    # 设置Ziegler-Natta反应动力学参数
    # 反应类型, 催化位点, 组分1, 组分2, 前指因子, 活化能, 反应级数, 终端双键分数, 参考温度
    r1 = [["Act-Spon", 1, "Titanium Tetrachloride", None, 0.08, 0, 1, None, 1e35],
        ["Act-Spon", 2, "Titanium Tetrachloride", None, 0.08, 0, 1, None, 1e35],
        ["Act-Spon", 3, "Titanium Tetrachloride", None, 0, 0, 1, None, 1e35],
        ["Act-Spon", 4, "Titanium Tetrachloride", None, 0, 0, 1, None, 1e35],
        ["Act-Cocat", 1, "Titanium Tetrachloride", "Triethyl Aluminum", 0.15/1000, 0, 1, None, 1e35],
        ["Act-Cocat", 2, "Titanium Tetrachloride", "Triethyl Aluminum", 0.15/1000, 0, 1, None, 1e35],
        ["Act-Cocat", 3, "Titanium Tetrachloride", "Triethyl Aluminum", 0, 0, 1, None, 1e35],
        ["Act-Cocat", 4, "Titanium Tetrachloride", "Triethyl Aluminum", 0, 0, 1, None, 1e35],
        ["Chain-Ini", 1, "Ethylene", None, 255/1000, 0, 1, None, 1e35],
        ["Chain-Ini", 2, "Ethylene", None, 90/1000,   0, 1, None, 1e35],
        ["Chain-Ini", 3, "Ethylene", None, 0,  0, 1, None, 1e35],
        ["Chain-Ini", 4, "Ethylene", None, 0,  0, 1, None, 1e35],
        ["Propagation", 1, "Ethylene-R", "Ethylene", 255/1000,  0, 1, None, 1e35],
        ["Propagation", 2, "Ethylene-R", "Ethylene", 90/1000,  0, 1, None, 1e35],
        ["Propagation", 3, "Ethylene-R", "Ethylene", 0,  0, 1, None, 1e35],
        ["Propagation", 4, "Ethylene-R", "Ethylene", 0,  0, 1, None, 1e35],
        ["Chat-Mon", 1, "Ethylene-R", "Ethylene", 0.09/1000,  0, 1, None, 1e35],
        ["Chat-Mon", 2, "Ethylene-R", "Ethylene", 0.24/1000,  0, 1, None, 1e35],
        ["Chat-Mon", 3, "Ethylene-R", "Ethylene", 0,  0, 1, None, 1e35],
        ["Chat-Mon", 4, "Ethylene-R", "Ethylene", 0,  0, 1, None, 1e35],
        ["Chat-H2", 1, "Ethylene-R", "Hydrogen", 5.55/1000,  0, 1, None, 1e35],
        ["Chat-H2", 2, "Ethylene-R", "Hydrogen", 18.5/1000,  0, 1, None, 1e35],
        ["Chat-H2", 3, "Ethylene-R", "Hydrogen", 0,  0, 1, None, 1e35],
        ["Chat-H2", 4, "Ethylene-R", "Hydrogen", 0,  0, 1, None, 1e35],
        ["Chat-Spon", 1, "Ethylene-R", None, 0.004,  0, 1, None, 1e35],
        ["Chat-Spon", 2, "Ethylene-R", None, 0.012,  0, 1, None, 1e35],
        ["Chat-Spon", 3, "Ethylene-R", None, 0,  0, 1, None, 1e35],
        ["Chat-Spon", 4, "Ethylene-R", None, 0,  0, 1, None, 1e35],
        ["Deact-Spon", 1, None, None, 0.0001,  0, 1, None, 1e35],
        ["Deact-Spon", 2, None, None, 0.0006,  0, 1, None, 1e35],
        ["Deact-Spon", 3, None, None, 0,  0, 1, None, 1e35],
        ["Deact-Spon", 4, None, None, 0,  0, 1, None, 1e35]]


    # 设置反应器参数
    specs = {
        "Pressure": 20265000,
        "Temperature": 433.15,
        "Valid phases": "liquid",
        "Specification type": "Reactor volume",
        "Reactor volume": 60,
        "Reactions": {"Type":"Z-N", 
                        "Species":species, 
                        "Reactions":r1,
                        "Reacting phase": "liquid"
                    },

        "Flash params": {"Max iter":100, "Tol":1e-7}
    }
    
    params = {}
    params["Pure Components"] = pure_args
    
    # 输出出口流量、组成、各位点各阶矩
    # 运行计算
    print("\nCSTR Product")
    
    res = cstr_run(components, params, feed, specs, "PC-SAFT")
    print("res: ", res )
    print("")
    product = parse_result(6, 4, res)
    print(product)