"""
聚合物体系的焓/熵计算
聚合物的性质由组成聚合物的链段表示, 由于本程序只考虑均聚物, 因此聚合物的性质可以直接由链段性质替代
(如果是共聚物, 则还需要链段分数, 通过加权计算求得共聚的性质)

物性数据来源:
1.小分子的理想气体摩尔生成焓/熵可以直接从json文件中读取,
2.小分子比热容的数据json文件和Aspen所用的公式和参数并不相同。为了保持与Aspen的一致, 
   案例所用的比热容公式(编号107)和数据来自于Aspen
3.链段的理想气体摩尔生成焓/熵通过基团贡献法获取, 基团数据来自于文献[1]
4.链段比热容的数据来自于Aspen

此程序包含:
1.纯物质焓、熵计算
2.单相混合物焓熵计算
3.多相混合物焓熵计算
4.相态判断(vap, liq, mix)
"""


import numpy as np
from polymer_model.thermodynamics.thermo_property import enth_ig_comp, entr_ig_comp
from polymer_model.thermodynamics.flash_algorithm import *
import pcsaft
from polymer_model.utility import InputError, R
import copy
from polymer_model.thermodynamics.args import args


def transform_args(components, args, stream):
    """
    根据体系包含的组分components, 从物性包args中提取参数, 包装成状态方式包PCSAFT所需的参数形式
    """

    NC = len(stream["Mole Fraction"].keys())    # 组分数
    
    # 将各组分参数包装为np.array的形式, 以方便求和运算
    m = np.zeros(NC)      # 各组分链段数, Aspen用PCSFTM表示
    s = np.zeros(NC)      # 各组分链段直径, Aspen用PCSFTV表示
    e = np.zeros(NC)      # 各组分能量参数, Aspen用PCSFTU表示
    e_assoc = np.zeros(NC)  # 关联成分的关联能量,单位K, Aspen用PCSFAU表示
    vol_a = np.zeros(NC)    # 关联成分的有效关联量, Aspen用PCSFAV表示
    w = np.zeros(NC)      # 各组分偏心因子
    Tc = np.zeros(NC)     # 各组分临界温度
    Pc = np.zeros(NC)     # 各组分临界压力
    mw = np.zeros(NC)     # 相对分子质量

    i = -1
    for component in stream["Mole Fraction"].keys():
        i = i + 1
        w[i] = args[component]["w"]
        Tc[i] = args[component]["Tc"]
        Pc[i] = args[component]["Pc"]
        mw[i] = args[component]["mw"]
        # 如果该组分为普通组分
        if components[component]["type"] == "conventional":
            m[i] = args[component]["m"]
            s[i] = args[component]["s"]
            e[i] = args[component]["e"]
            if args[component].get("e_assoc") != None:
                e_assoc[i] = args[component].get("e_assoc")
            else:
                e_assoc[i] = 0
            if args[component].get("vol_a") != None:
                vol_a[i] = args[component].get("vol_a")
            else:
                vol_a[i] = 0
        # 如果该组分为聚合物
        elif components[component]["type"] == "polymer":
            DPN = stream[component]["DPN"]  # 获取数均聚合度

            Ns = 0   # 链段的种类数
            for segment in stream[component]["SFRAC"].keys():
                if stream[component]["SFRAC"][segment] != 0 or stream[component]["SFRAC"][segment] != None:
                    Ns = Ns + 1

            sfrac = np.zeros(Ns)   # 各链段分数
            mw_seg = np.zeros(Ns)  # 各链段分子质量
            r_seg = np.zeros(Ns)   # 各链段的链段数/数均分子量
            s_seg = np.zeros(Ns)   # 各链段的链段直径
            e_seg = np.zeros(Ns)   # 各链段的能量参数
            e_assoc_seg = np.zeros(Ns)
            vol_a_seg = np.zeros(Ns)

            j = -1
            for segment in stream[component]["SFRAC"].keys():
                if stream[component]["SFRAC"][segment] != 0 or stream[component]["SFRAC"][segment] != None:
                    j = j + 1
                    sfrac[j] = stream[component]["SFRAC"][segment]
                    mw_seg[j] = args[segment]["mw"]
                    r_seg[j] = args[segment]["r"]
                    s_seg[j] = args[segment]["s"]
                    e_seg[j] = args[segment]["e"]
                    if args[segment].get("e_assoc") != None:
                        e_assoc_seg[j] = args[segment].get("e_assoc")
                    else:
                        e_assoc_seg[j] = 0
                    if args[segment].get("vol_a") != None:
                        vol_a_seg[j] = args[segment].get("vol_a")
                    else:
                        vol_a_seg[j] = 0
            MWN = sfrac*mw_seg*DPN
            m[i] = np.sum(r_seg*MWN) 
            s[i] = np.sum(sfrac*s_seg)
            e[i] = np.sum(sfrac*e_seg)
            e_assoc[i] = np.sum(sfrac*e_assoc_seg)
            vol_a[i] = np.sum(sfrac*vol_a_seg)

    pyargs = {"Tc": Tc, "Pc": Pc, "w": w, "m": m, "s": s, "e": e, 
            "e_assoc": e_assoc, "vol_a": vol_a}

    return pyargs


def mix(components, catalyst, args, streams, h_tol=0.01):
    """
    混合器的主要方法：输入多个物流, 计算出口物流的流率、摩尔组成、温度、压力

    h_tol: 摩尔焓的误差值
    """

    NC = len(streams[0]["Mole Fraction"].keys())  # 组分数 
    component_list = list(streams[0]["Mole Fraction"].keys()) 
    if catalyst == None:
        catalyst_list = []
    else:
        catalyst_list = list(catalyst.keys())   # 催化剂列表
    polymer_list = []   # 聚合物列表
    i = 0    # 组分数
    for component in components.keys():
        if components[component]["type"] == "conventional":
            i += 1
        elif components[component]["type"] == "polymer":
            polymer_list.append(component)
            i += 1

    Fc = np.zeros(NC)         # 各组分的摩尔流量, mol/s
    P_min = streams[0]["P"]   # 出口压力等于最小的入口压力 
    H = 0                     # 总焓值, 混合气入口和出口等焓
    T0 = 0                    # 计算出口温度的初始值
    T_min = streams[0]["T"]   # 出口温度的最小值
    T_max = streams[0]["T"]   # 出口温度的最大值

    product = copy.deepcopy(streams[0])
    for stream in streams:
        # 获取每条流股的摩尔流量、摩尔组成、摩尔焓
        Fs = stream["Mole Flow"]
        x = np.array(list(stream["Mole Fraction"].values())) 
        Fc = Fc + Fs*x
        h = molar_enthalpy_of_mixture(components, args, stream)[1]
        H = H + Fs*h
        T0 = T0 + Fs*stream["T"]
        if stream["P"] < P_min:
            P_min = stream["P"]
        if stream["T"] > T_max:
            T_max = stream["T"]
        if stream["T"] < T_min:
            T_min = stream["T"]

        for key in stream.keys():
            # 将各流股催化剂的各位点流量加和
            if key in catalyst_list:
                if key not in product:
                    product[key] = dict()
                    product[key]["CPSFLOW"] = stream[key]["CPSFLOW"]
                    product[key]["CDSFLOW"] = stream[key]["CDSFLOW"]
                    product[key]["CISFLOW"] = stream[key]["CISFLOW"]
                    product[key]["CVSFLOW"] = stream[key]["CVSFLOW"]
                else:
                    product[key]["CPSFLOW"] += stream[key]["CPSFLOW"]
                    product[key]["CDSFLOW"] += stream[key]["CDSFLOW"]
                    product[key]["CISFLOW"] += stream[key]["CISFLOW"]
                    product[key]["CVSFLOW"] += stream[key]["CVSFLOW"]

            # 将各流股聚合物的各位的零阶矩和一阶矩流量加和
            if key in polymer_list:
                if key not in product:
                    product[key] = dict()
                    product[key]["ZMOM"] = stream[key]["ZMOM"]
                    product[key]["FMOM"] = stream[key]["FMOM"]
                    product[key]["DPN"] = product[key]["FMOM"]/product[key]["ZMOM"]
                    product[key]["SFLOW"] = dict()
                    product[key]["SFRAC"] = dict()
                    for segment in stream[key]["SFLOW"].keys():
                        product[key]["SFLOW"][segment] = stream[key]["SFLOW"][segment]
                        product[key]["SFRAC"][segment] = product[key]["SFLOW"][segment]/product[key]["FMOM"]

                else:
                    product[key]["ZMOM"] += stream[key]["ZMOM"]
                    product[key]["FMOM"] += stream[key]["FMOM"]
                    product[key]["DPN"] = product[key]["FMOM"]/product[key]["ZMOM"]
                    for segment in stream[key]["SFLOW"].keys():
                        product[key]["SFLOW"][segment] += stream[key]["SFLOW"][segment]
                        product[key]["SFRAC"][segment] = product[key]["SFLOW"][segment]/product[key]["FMOM"]
        
    F = np.sum(Fc)   # 出口总摩尔流量
    x = Fc/F    # 出口摩尔分数
    P = P_min   # 出口压力
    Hm = H/F    # 出口摩尔焓
    T0 = T0/F   # 计算出口温度的初始值

    product["P"] = P
    product["Mole Flow"] = F
    for i in range(NC):
        product["Mole Fraction"][component_list[i]] = x[i]

    # 迭代, 更新出口温度的上下界
    while True:
        product["T"] = T_min
        Hm_cal = molar_enthalpy_of_mixture(components, args, product)[1]
        if Hm_cal == Hm:
            return product
        elif Hm_cal < Hm:
            break
        else:
            T_max = T_min
            T_min = 0.9*T_min
    
    # 二分法迭代, 计算出口温度
    error = 2*h_tol
    while True:
        T = (T_max+T_min)/2
        product["T"] = T
        Hm_cal = molar_enthalpy_of_mixture(components, args, product)[1]
        error = Hm_cal - Hm
        if abs(error) <= h_tol:
            break
        if error > 0:
            T_max = T
        else:
            T_min = T
    return product


def flash(components, catalyst, args, streams, flash_type, valid_phases="vapor-liquid", tol=1e-7, max_iter=100):
    """
    闪蒸单元的主要方法：输入多个物流, 计算出口气液物流的流率、摩尔组成、温度、压力
    """

    NC = len(streams[0]["Mole Fraction"].keys())  # 组分数 
    component_list = list(streams[0]["Mole Fraction"].keys()) 
    if catalyst == None:
        catalyst_list = []
    else:
        catalyst_list = list(catalyst.keys())   # 催化剂列表
    polymer_list = []   # 聚合物列表
    polymer_index = []  # 聚合物索引
    i = 0    # 组分数
    for component in components.keys():
        if components[component]["type"] == "conventional":
            i += 1
        elif components[component]["type"] == "polymer":
            polymer_list.append(component)
            polymer_index.append(i)
            i += 1

    Fc = np.zeros(NC)    # 各组分的摩尔流量, mol/s
    H_in = 0             # 入口总焓值

    # 遍历各条流股
    feed = dict()
    for stream in streams:
        # 获取每条流股的摩尔流量、摩尔组成、摩尔焓
        Fs = stream["Mole Flow"]
        z = np.array(list(stream["Mole Fraction"].values())) 
        Fc = Fc + Fs*z
        h = molar_enthalpy_of_mixture(components, args, stream)[1]
        H_in = H_in + Fs*h

        for key in stream.keys():
            # 将各流股催化剂的各位点流量加和
            if key in catalyst_list:
                if key not in feed:
                    feed[key] = dict()
                    feed[key]["CPSFLOW"] = stream[key]["CPSFLOW"]
                    feed[key]["CDSFLOW"] = stream[key]["CDSFLOW"]
                    feed[key]["CISFLOW"] = stream[key]["CISFLOW"]
                    feed[key]["CVSFLOW"] = stream[key]["CVSFLOW"]
                else:
                    feed[key]["CPSFLOW"] += stream[key]["CPSFLOW"]
                    feed[key]["CDSFLOW"] += stream[key]["CDSFLOW"]
                    feed[key]["CISFLOW"] += stream[key]["CISFLOW"]
                    feed[key]["CVSFLOW"] += stream[key]["CVSFLOW"]

            # 将各流股聚合物的各位的零阶矩和一阶矩流量加和
            if key in polymer_list:
                if key not in feed:
                    feed[key] = dict()
                    feed[key]["ZMOM"] = stream[key]["ZMOM"]
                    feed[key]["FMOM"] = stream[key]["FMOM"]
                    feed[key]["DPN"] = feed[key]["FMOM"]/feed[key]["ZMOM"]
                    feed[key]["SFLOW"] = dict()
                    feed[key]["SFRAC"] = dict()
                    for segment in stream[key]["SFLOW"].keys():
                        feed[key]["SFLOW"][segment] = stream[key]["SFLOW"][segment]
                        feed[key]["SFRAC"][segment] = feed[key]["SFLOW"][segment]/feed[key]["FMOM"]

                else:
                    feed[key]["ZMOM"] += stream[key]["ZMOM"]
                    feed[key]["FMOM"] += stream[key]["FMOM"]
                    feed[key]["DPN"] = feed[key]["FMOM"]/feed[key]["ZMOM"]
                    for segment in stream[key]["SFLOW"].keys():
                        feed[key]["SFLOW"][segment] += stream[key]["SFLOW"][segment]
                        feed[key]["SFRAC"][segment] = feed[key]["SFLOW"][segment]/feed[key]["FMOM"]
        
    F = np.sum(Fc)   # 入口总摩尔流量
    z = Fc/F    # 出口摩尔分数

    feed["Mole Fraction"] = dict()
    for i in range(NC):
        feed["Mole Fraction"][component_list[i]] = z[i]

    # 换算参数
    flash_args = transform_args(components, args, feed)

    # 闪蒸规范为“温度-压力”
    if flash_type.get("T") != None and flash_type.get("P") != None:
        T = flash_type.get("T")
        P = flash_type.get("P")
        for index in polymer_index:
            # 换算聚合物
            if "DPN" in feed[component_list[index]].keys():
                z[index] = z[index]/feed[component_list[index]]["DPN"]
        z = z/np.sum(z)        
        # 调用闪蒸
        beta, x, y, iter = VLE_TPflash(T, P, z, flash_args, tol, max_iter)

        # 重新换算成链段
        for index in polymer_index:
            # 换算聚合物
            if "DPN" in feed[component_list[index]].keys():
                x[index] = x[index]*feed[component_list[index]]["DPN"]
        beta = beta/(beta+(1-beta)*sum(x))        
        x = x/sum(x)    # 新的液相组成

        # 将气液相变量封装为字典形式
        liquid = copy.deepcopy(feed)
        liquid["T"] = T
        liquid["P"] = P
        liquid["Mole Flow"] = F*(1-beta)
        liquid["Valid phases"] = "liquid"

        vapor = dict()
        vapor["T"] = T
        vapor["P"] = P
        vapor["Mole Flow"] = F*beta
        vapor["Valid phases"] = "vapor"
        vapor["Mole Fraction"] = dict()
        for i in range(NC):
            vapor["Mole Fraction"][component_list[i]] = y[i]
            liquid["Mole Fraction"][component_list[i]] = x[i]
        
    H_out = F*beta*molar_enthalpy_of_mixture(components, args, vapor)[1] + \
            F*(1-beta)*molar_enthalpy_of_mixture(components, args, liquid)[1]

    duty = H_out - H_in
    
    return liquid, vapor, beta, duty


def heat(components, catalyst, args, streams, flash_type, valid_phases="vapor-liquid", max_iter=100, tol=0.001):
    """
    加热器/冷却器的主要方法：给定若干进口物流和加热器的出口规范, 计算出口的温度或者热负荷

    Params:
    ------------
    components: dict

    streams: dict
    
    args: dict

    flash_type: dicts, 
        Keywords include temperature, pressure, heat load, vaporization fraction, 
        temperature change、Degree of Superheating、Degree of Supercooling

    max_iter: int
        max iteration
    tol : float, default 0.001
        tolerance

    Returns:
    ------------
    product: 出口物流
    duty: 加热器/冷却器热负荷
    """

    # 获取入口物流混合后的摩尔流量、摩尔组成、总焓值
    NC = len(streams[0]["Mole Fraction"].keys())  # 组分数 
    component_list = list(streams[0]["Mole Fraction"].keys()) 
    if catalyst == None:
        catalyst_list = []
    else:
        catalyst_list = list(catalyst.keys())   # 催化剂列表
    polymer_list = []   # 聚合物列表
    i = 0    # 组分数
    for component in components.keys():
        if components[component]["type"] == "conventional":
            i += 1
        elif components[component]["type"] == "polymer":
            polymer_list.append(component)
            i += 1

    Fc = np.zeros(NC)         # 各组分的摩尔流量, mol/s
    H_in = 0                  # 入口总焓值
    P_min = streams[0]["P"]   # 获取最小的入口压力 
    
    product = copy.deepcopy(streams[0])
    for stream in streams:
        # 获取每条流股的摩尔流量、摩尔组成、摩尔焓
        Fs = stream["Mole Flow"]
        x = np.array(list(stream["Mole Fraction"].values())) 
        Fc = Fc + Fs*x
        h = molar_enthalpy_of_mixture(components, args, stream)[1]
        H_in = H_in + Fs*h
        if stream["P"] < P_min:
            P_min = stream["P"]
        for key in stream.keys():
            # 将各流股催化剂的各位点流量加和
            if key in catalyst_list:
                if key not in product:
                    product[key] = dict()
                    product[key]["CPSFLOW"] = stream[key]["CPSFLOW"]
                    product[key]["CDSFLOW"] = stream[key]["CDSFLOW"]
                    product[key]["CISFLOW"] = stream[key]["CISFLOW"]
                    product[key]["CVSFLOW"] = stream[key]["CVSFLOW"]
                else:
                    product[key]["CPSFLOW"] += stream[key]["CPSFLOW"]
                    product[key]["CDSFLOW"] += stream[key]["CDSFLOW"]
                    product[key]["CISFLOW"] += stream[key]["CISFLOW"]
                    product[key]["CVSFLOW"] += stream[key]["CVSFLOW"]

            # 将各流股聚合物的各位的零阶矩和一阶矩流量加和
            if key in polymer_list:
                if key not in product:
                    product[key] = dict()
                    product[key]["ZMOM"] = stream[key]["ZMOM"]
                    product[key]["FMOM"] = stream[key]["FMOM"]
                    product[key]["DPN"] = product[key]["FMOM"]/product[key]["ZMOM"]
                    product[key]["SFLOW"] = dict()
                    product[key]["SFRAC"] = dict()
                    for segment in stream[key]["SFLOW"].keys():
                        product[key]["SFLOW"][segment] = stream[key]["SFLOW"][segment]
                        product[key]["SFRAC"][segment] = product[key]["SFLOW"][segment]/product[key]["FMOM"]

                else:
                    product[key]["ZMOM"] += stream[key]["ZMOM"]
                    product[key]["FMOM"] += stream[key]["FMOM"]
                    product[key]["DPN"] = product[key]["FMOM"]/product[key]["ZMOM"]
                    for segment in stream[key]["SFLOW"].keys():
                        product[key]["SFLOW"][segment] += stream[key]["SFLOW"][segment]
                        product[key]["SFRAC"][segment] = product[key]["SFLOW"][segment]/product[key]["FMOM"]
        
    F = np.sum(Fc)   # 出口总摩尔流量
    x = Fc/F    # 出口摩尔分数

    # 初始化出口物流
    product["Mole Flow"] = F
    for i in range(NC):
        product["Mole Fraction"][component_list[i]] = x[i]
    product["Valid phases"] = valid_phases

    # 出口规范为“温度-压力”
    if flash_type.get("T") != None and flash_type.get("P") != None:
        product["T"] = flash_type.get("T")
        # 如果P>0,即给定出口压力；如果P<=0,即给定压降
        if flash_type.get("P") > 0:
            product["P"] = flash_type.get("P")
        else:
            product["P"] = P_min - flash_type.get("P")
        
        H_out = F*molar_enthalpy_of_mixture(components, args, product)[1]
        duty = H_out - H_in
    # 出口规范为“负荷-压力”
    elif flash_type.get("Duty") != None and flash_type.get("P") != None:
        pass
    
    return product, duty


def  multi_stream_heat_exchang(components, catalyst, args, hot_streams, cold_streams, specs, h_tol=0.001):
    """
    多物流热交换器的主要方法：给定若干进口物流和加热器的出口规范, 计算出口的温度或者热负荷
    在指定流股的出口规范时: 必须指定所有热流股的出口, 或者所有冷流股的出口,
    至少有一条流股的出口是待计算的

    Params:
    ------------
    components: dict

    streams: dict
    
    args: dict

    spec: dicts, 
        Keywords include temperature, pressure, heat load, vaporization fraction, 
        temperature change、Degree of Superheating、Degree of Supercooling

    max_iter: int
        max iteration
    tol : float, default 0.001
        tolerance

    Returns:
    ------------
    out: 出口物流
    duty: 加热器/冷却器热负荷
    """

    # 每条冷热流股的出口和传热量 
    hot_out = [0 for i in range(len(hot_streams))]
    hot_duty = [0 for i in range(len(hot_streams))]
    cold_out = [0 for i in range(len(cold_streams))]
    cold_duty = [0 for i in range(len(cold_streams))]
    undetermined = []  # 未指定的流股
    duty = 0  # 换热器总负荷

    # 判断是否所有的热流股都指定了出口
    flag = True
    for spec in specs["hot"]:
        if spec["spec"] == None:
            flag = False
    
    # 如果所有热流股都指定了出口
    if flag:
        T_max = 0   # 入口热流股温度的最大值，也是出口温度的上界
        # 计算每条热流股的换热量
        for i in range(len(hot_streams)):
            # 获取每条热流股的出口规范、压力、进料条件、有效相、最大迭代数、容差
            flash_type = copy.deepcopy(specs["hot"][i]["spec"])
            flash_type["P"] = specs["hot"][i]["P"]
            stream = hot_streams[i]
            valid_phases = specs["hot"][i]["Valid phases"]
            max_iter = specs["hot"][i]["Max iter"]
            tol = specs["hot"][i]["Tol"]
            # 调用加热器
            out, stream_duty = heat(components, catalyst, args, [stream], flash_type, valid_phases, max_iter, tol)
            hot_out[i] = out
            hot_duty[i] = stream_duty
            if stream["T"] > T_max:
                T_max = stream["T"]
        
        T_min = 0  # 未指定的冷流股入口温度最大值，也是其出口温度的下界
        # 获取冷流股的出口规范, 并计算冷流股的吸热量
        for i in range(len(cold_streams)):
            if (specs["cold"][i]["spec"] != None):
                # 获取每条热流股的出口规范、压力、进料条件、有效相、最大迭代数、容差
                flash_type = copy.deepcopy(specs["cold"][i]["spec"])
                flash_type["P"] = specs["cold"][i]["P"]
                stream = cold_streams[i]
                valid_phases = specs["cold"][i]["Valid phases"]
                max_iter = specs["cold"][i]["Max iter"]
                tol = specs["cold"][i]["Tol"]
                # 调用加热器
                out, stream_duty = heat(components, catalyst, args, [stream], flash_type, valid_phases, max_iter, tol)
                cold_out[i] = out
                cold_duty[i] = stream_duty
            else:
                if cold_streams[i]["T"] > T_min:
                    T_min = cold_streams[i]["T"]
                # 第i条流股未指定出口
                undetermined.append(i)  
        duty = np.sum(hot_duty)
        remaining_duty = 0 - np.sum(hot_duty) - np.sum(cold_duty)

        # 二分法迭代, 计算未指定流股的出口温度
        while True:
            T = (T_min+T_max)/2
            calc_duty = 0
            for i in undetermined:
                # 获取每条热流股的出口规范、压力、进料条件、有效相、最大迭代数、容差
                flash_type = dict()
                flash_type["T"] = T
                flash_type["P"] = specs["cold"][i]["P"]
                stream = cold_streams[i]
                valid_phases = specs["cold"][i]["Valid phases"]
                max_iter = specs["cold"][i]["Max iter"]
                tol = specs["cold"][i]["Tol"]
                # 调用加热器
                out, stream_duty = heat(components, catalyst, args, [stream], flash_type, valid_phases, max_iter, tol)
                cold_out[i] = out
                cold_duty[i] = stream_duty
                calc_duty += stream_duty
            error = calc_duty - remaining_duty
            if abs(error) <= h_tol:
                break
            if error > 0:
                T_max = T
            else:
                T_min = T
    
    # 如果有热流股未指定出口，那么所有的冷流股应该指定了出口
    else:
        T_min = 2000   # 入口冷流股温度的最小值，也是出口温度的下界
        # 计算每条冷流股的换热量
        for i in range(len(cold_streams)):
            # 获取每条冷流股的出口规范、压力、进料条件、有效相、最大迭代数、容差
            flash_type = copy.deepcopy(specs["cold"][i]["spec"])
            flash_type["P"] = specs["cold"][i]["P"]
            stream = cold_streams[i]
            valid_phases = specs["cold"][i]["Valid phases"]
            max_iter = specs["cold"][i]["Max iter"]
            tol = specs["cold"][i]["Tol"]
            # 调用加热器
            out, stream_duty = heat(components, catalyst, args, [stream], flash_type, valid_phases, max_iter, tol)
            cold_out[i] = out
            cold_duty[i] = stream_duty
            if stream["T"] <  T_min:
                T_min = stream["T"]
        
        T_max = 2000  # 未指定的热流股入口温度最小值，也是其出口温度的上界
        # 获取冷流股的出口规范, 并计算冷流股的吸热量
        for i in range(len(hot_streams)):
            if (specs["hot"][i]["spec"] != None):
                # 获取每条热流股的出口规范、压力、进料条件、有效相、最大迭代数、容差
                flash_type = copy.deepcopy(specs["hot"][i]["spec"])
                flash_type["P"] = specs["hot"][i]["P"]
                stream = hot_streams[i]
                valid_phases = specs["hot"][i]["Valid phases"]
                max_iter = specs["hot"][i]["Max iter"]
                tol = specs["hot"][i]["Tol"]
                # 调用加热器
                out, stream_duty = heat(components, catalyst, args, [stream], flash_type, valid_phases, max_iter, tol)
                hot_out[i] = out
                hot_duty[i] = stream_duty
            else:
                if hot_streams[i]["T"] < T_max:
                    T_max = hot_streams[i]["T"]
                # 第i条流股未指定出口
                undetermined.append(i)  
        duty = np.sum(cold_duty)
        remaining_duty = 0 - np.sum(hot_duty) - np.sum(cold_duty)

        # 二分法迭代, 计算未指定流股的出口温度
        while True:
            T = (T_min+T_max)/2
            calc_duty = 0
            for i in undetermined:
                # 获取每条热流股的出口规范、压力、进料条件、有效相、最大迭代数、容差
                flash_type = dict()
                flash_type["T"] = T
                flash_type["P"] = specs["hot"][i]["P"]
                stream = hot_streams[i]
                valid_phases = specs["hot"][i]["Valid phases"]
                max_iter = specs["hot"][i]["Max iter"]
                tol = specs["hot"][i]["Tol"]
                # 调用加热器
                out, stream_duty = heat(components, catalyst, args, [stream], flash_type, valid_phases, max_iter, tol)
                hot_out[i] = out
                hot_duty[i] = stream_duty
                calc_duty += stream_duty
            error = calc_duty - remaining_duty
            if abs(error) <= h_tol:
                break
            if error > 0:
                T_max = T
            else:
                T_min = T

    return hot_out, hot_duty, cold_out, cold_duty, -duty 


# 等熵压缩
def compress_isentropic(components, catalyst, args, streams, pressure_spec, valid_phases="vapor", 
            s_tol=0.001, h_tol=0.001, isentropic_efficiency=0.72, mechanical_efficiency=1):
    """
    等熵压缩机的主要方法：输入流股、压缩机参数, 计算出口物流的温度
    
    Params
    ---------
    pressure_spec: dict
        用于计算出口压力,出口压力的规范有3种, 选择一个即可
        "Discharge pressure": 出口压力
        "Pressure increase": 压差
        "Pressure ratio": 压力比
    efficiency: 等熵压缩的效率, 默认0.72
    valid_phases: 有效相, 默认是气体
    s_tol: 摩尔熵的误差值

    Returns
    ---------
    out: 出口流股
    power: 压缩机所需功率
    """

    stream = mix(components, catalyst, args, streams, h_tol) # 混合后的流股

    P_in = stream["P"]   # 入口压力
    if pressure_spec.get("Discharge pressure") != None:
        P_out = pressure_spec.get("Discharge pressure")
    elif pressure_spec.get("Pressure increase") != None:
        P_out = P_in + pressure_spec.get("Pressure increase")
    elif pressure_spec.get("Pressure ratio") != None:
        P_out = P_in*pressure_spec.get("Pressure ratio")
    else:
        raise InputError("Wrong pressure specification! Please enter 'Discharge pressure' \
            or 'Pressure increase' or 'Pressure ratio' ")

    # 设置流股的有效相
    stream["Valid phases"] = valid_phases
    # 计算入口物流的摩尔焓
    _, H_in, S_in = molar_enthalpy_of_mixture(components, args, stream)
    T_in = stream["T"]   # 入口温度
    print("S_in: ", S_in)

    # 初始化出口物流
    out = copy.deepcopy(stream)
    out["P"] = P_out

    # 迭代, 更新等熵出口温度的上下界
    T_min = T_in
    T_max = T_min*1.1
    while True:
        out["T"] = T_max
        S_mix = molar_enthalpy_of_mixture(components, args, out)[2]
        if S_mix < S_in:
            T_min = T_max
            T_max = 1.1*T_max
        else:
            break

    # 二分法迭代, 计算等熵出口温度
    error = 2*s_tol
    while True:
        T_isentropic = (T_min + T_max)/2   # 等熵出口温度
        out["T"] = T_isentropic
        _, H_isentropic, S_isentropic = molar_enthalpy_of_mixture(components, args, out)
        error = S_isentropic - S_in
        if abs(error) <= s_tol:
            break
        if error > 0:
            T_max = T_isentropic
        else:
            T_min = T_isentropic
    
    mole_flow = out["Mole Flow"]   # 摩尔流率, mol/s
    # 计算等熵压缩功(可逆绝热压缩功)
    isentropic_compression_work = (H_isentropic - H_in)*mole_flow
    # 计算真实压缩功, aspen用Indicated horsepower表示
    true_compression_work = isentropic_compression_work/isentropic_efficiency
    # 计算压缩机所需做功, aspen用Brake horsepower表示
    net_work = true_compression_work/mechanical_efficiency
    # 计算出口焓值
    H_out = (H_in*mole_flow + true_compression_work)/mole_flow

    # 根据出口焓值, 迭代计算真实的出口温度
    # 迭代, 更新出口温度的上下界
    T_min = T_isentropic
    T_max = T_min+1
    while True:
        out["T"] = T_max
        _, H_mix, _ = molar_enthalpy_of_mixture(components, args, out)
        if H_mix < H_out:
            T_min = T_max
            T_max = T_min+1
        else:
            break

    # 二分法迭代, 计算出口温度
    error = 2*s_tol
    while True:
        T_out = (T_min + T_max)/2   # 等熵出口温度
        out["T"] = T_out
        _, H_mix, _ = molar_enthalpy_of_mixture(components, args, out)
        error = H_mix - H_out
        if abs(error) <= h_tol:
            break
        if error > 0:
            T_max = T_out
        else:
            T_min = T_out

    return out, T_isentropic, isentropic_compression_work, true_compression_work, net_work


def molar_enthalpy_of_mixture(components, args, stream):
    """
    计算混合物的摩尔焓, 需要先判断流股的相态: 如果是单相, 可直接计算
    如果是气液两相, 则需要进行闪蒸计算, 分成气液两相进行计算, 再加和
    """

    T = stream["T"]   # 温度
    P = stream["P"]   # 压力

    NC = 0    # 组分数
    for component in stream["Mole Fraction"].keys():
        if stream["Mole Fraction"][component] != 0 and stream["Mole Fraction"][component] != None:
            NC += 1
    
    # 将各组分参数包装为np.array的形式, 以方便求和运算
    z = np.zeros(NC)      # 各组分摩尔分数
    m = np.zeros(NC)      # 各组分链段数, Aspen用PCSFTM表示
    s = np.zeros(NC)      # 各组分链段直径, Aspen用PCSFTV表示
    e = np.zeros(NC)      # 各组分能量参数, Aspen用PCSFTU表示
    e_assoc = np.zeros(NC)  # 关联成分的关联能量,单位K, Aspen用PCSFAU表示
    vol_a = np.zeros(NC)    # 关联成分的有效关联量, Aspen用PCSFAV表示
    w = np.zeros(NC)      # 各组分偏心因子
    Tc = np.zeros(NC)     # 各组分临界温度
    Pc = np.zeros(NC)     # 各组分临界压力
    mw = np.zeros(NC)     # 相对分子质量
    H_ig_comp = np.zeros(NC)  # 各组分的理想焓
    S_ig_comp = np.zeros(NC)  # 各组分的理想熵

    i = -1
    has_polymer = False
    for component in stream["Mole Fraction"].keys():
        if stream["Mole Fraction"][component] != 0 and stream["Mole Fraction"][component] != None:
            i = i + 1
            z[i] = stream["Mole Fraction"][component]
            w[i] = args[component]["w"]
            Tc[i] = args[component]["Tc"]
            Pc[i] = args[component]["Pc"]
            mw[i] = args[component]["mw"]
            if components[component]["type"] == "conventional":
                m[i] = args[component]["m"]
                s[i] = args[component]["s"]
                e[i] = args[component]["e"]
                if args[component].get("e_assoc") != None:
                    e_assoc[i] = args[component].get("e_assoc")
                else:
                    e_assoc[i] = 0
                if args[component].get("vol_a") != None:
                    vol_a[i] = args[component].get("vol_a")
                else:
                    vol_a[i] = 0
                eqno = args[component]["cp_mol_ig"]["eqno"]
                params = args[component]["cp_mol_ig"]["params"]
                enth_mol_form_ig_ref = args[component]["enth_mol_form_ig_ref"]
                entr_mol_form_ig_ref = args[component]["entr_mol_form_ig_ref"]
                H_ig_comp[i] = enth_ig_comp(T, enth_mol_form_ig_ref, eqno, params, Tref=298.15)
                S_ig_comp[i] = entr_ig_comp(T, P*z[i], entr_mol_form_ig_ref, eqno, params, Tref=298.15, Pref=101325)
            elif components[component]["type"] == "polymer":
                has_polymer = True
                polymer_index = i
                DPN = stream[component]["DPN"]

                Ns = 0   # 链段的种类数
                for segment in stream[component]["SFRAC"].keys():
                    if stream[component]["SFRAC"][segment] != 0 or stream[component]["SFRAC"][segment] != None:
                        Ns = Ns + 1

                sfrac = np.zeros(Ns)   # 各链段分数
                mw_seg = np.zeros(Ns)  # 各链段分子质量
                r_seg = np.zeros(Ns)   # 各链段的链段数/数均分子量
                s_seg = np.zeros(Ns)   # 各链段的链段直径
                e_seg = np.zeros(Ns)   # 各链段的能量参数
                e_assoc_seg = np.zeros(Ns)
                vol_a_seg = np.zeros(Ns)

                H_ig_seg = np.zeros(Ns)  # 各链段的理想焓
                S_ig_seg = np.zeros(Ns)  # 各链段的理想熵

                j = -1
                for segment in stream[component]["SFRAC"].keys():
                    if stream[component]["SFRAC"][segment] != 0 or stream[component]["SFRAC"][segment] != None:
                        j = j + 1
                        sfrac[j] = stream[component]["SFRAC"][segment]
                        mw_seg[j] = args[segment]["mw"]
                        r_seg[j] = args[segment]["r"]
                        s_seg[j] = args[segment]["s"]
                        e_seg[j] = args[segment]["e"]
                        if args[segment].get("e_assoc") != None:
                            e_assoc_seg[j] = args[segment].get("e_assoc")
                        else:
                            e_assoc_seg[j] = 0
                        if args[segment].get("vol_a") != None:
                            vol_a_seg[j] = args[segment].get("vol_a")
                        else:
                            vol_a_seg[j] = 0
                        eqno = args[segment]["cp_mol_ig"]["eqno"]
                        params = args[segment]["cp_mol_ig"]["params"]
                        enth_mol_form_ig_ref = args[segment]["enth_mol_form_ig_ref"]
                        entr_mol_form_ig_ref = args[segment]["entr_mol_form_ig_ref"]
                        H_ig_seg[j] = enth_ig_comp(T, enth_mol_form_ig_ref, eqno, params, Tref=298.15)
                        S_ig_seg[j] = entr_ig_comp(T, P*z[i]*sfrac[j], entr_mol_form_ig_ref, eqno, params, Tref=298.15, Pref=101325)

                mw_seg_avg = sum(sfrac*mw_seg) # 链段的平均相对分子质量
                MWN = sfrac*mw_seg*DPN
                m[i] = np.sum(r_seg*MWN) 
                s[i] = np.sum(sfrac*s_seg)
                e[i] = np.sum(sfrac*e_seg)
                e_assoc[i] = np.sum(sfrac*e_assoc_seg)
                vol_a[i] = np.sum(sfrac*vol_a_seg)

                H_ig_comp[i] = np.sum(sfrac*H_ig_seg)*mw[i]/mw_seg_avg  # 聚合物的理想焓(以链段表示)
                S_ig_comp[i] = np.sum(sfrac*S_ig_seg)*mw[i]/mw_seg_avg  # 聚合物的理想熵(以链段表示)

    H_ig_mix = np.sum(H_ig_comp*z)   # 混合物的理想焓
    S_ig_mix = np.sum(S_ig_comp*z)   # 混合物的理想熵

    pyargs = {"Tc": Tc, "Pc": Pc, "w": w, "m": m, "s": s, "e": e, 
            "e_assoc": e_assoc, "vol_a": vol_a}
    # 如果存在聚合物, 需要进行换算
    if has_polymer:
        z[polymer_index] = z[polymer_index]*mw[polymer_index]/mw_seg_avg/DPN
        z = z/sum(z)

        # 如果存在聚合物,则体系只会是液相或者气液
        
        Pb, _, state1 = bubble_pressure(T, z, pyargs)
        Tb, _, state2 = bubble_temperature(T, z, pyargs)

        # 大于泡点压力，处于液相
        if stream.get("Valid phases") == "liquid" or (state1 and P >= Pb) or (state2 and T <= Tb) :
            x = z
            y = z
            beta = 0
        else:
            beta, x, y, iter_total = VLE_TPflash(T, P, z, pyargs)
    else:
        # 有效相为液相
        if stream.get("Valid phases") == "liquid":
            x = z
            y = z
            beta = 0
        elif stream.get("Valid phases") == "vapor":
            x = z
            y = z
            beta = 1
        else:
            Pb, _, state1 = bubble_pressure(T, z, pyargs)
            Tb, _, state2 = bubble_temperature(P, z, pyargs)
            Pd, _, state3 = dew_pressure(T, z, pyargs)
            Td, _, state4 = dew_temperature(P, z, pyargs)

            # 考虑纯物质的情况(P=Pb=Pd)
            if abs(Pb-Pd)<0.1 and abs(P-Pb)<0.1:
                if (stream.get("Valid phases") == "vapor"):
                    beta = 1
                elif (stream.get("Valid phases") == "liquid"):
                    beta = 0
                else:
                    raise InputError("Please input valid phases of the stream")
                x = z
                y = z
            else:
                # 高于泡点压力，处于液相
                if (state1 and P >= Pb) or (state2 and T <= Tb):
                    x = z
                    y = z
                    beta = 0
                # 低于露点压力，处于气相
                elif (state3 and P <= Pd) or (state4 and T >= Td) :
                    x = z
                    y = z
                    beta = 1
                else:
                    beta, x, y, iter_total = VLE_TPflash(T, P, z, pyargs)
    
    rho_vap = pcsaft.pcsaft_den(T, P, y, pyargs, "vap")
    H_res_vap = pcsaft.pcsaft_hres(T, rho_vap, y, pyargs)
    S_res_vap = pcsaft.pcsaft_sres(T, rho_vap, y, pyargs)
    # 剩余熵还需减去压力导致的熵变,
    # S_res_vap = S_res_vap - R*np.sum(y*np.log(P*y/101325))
    for i in y:
        if i != 0:
            S_res_vap = S_res_vap - R*i*np.log(P*i/101325)

    rho_liq = pcsaft.pcsaft_den(T, P, x, pyargs, "liq")
    H_res_liq = pcsaft.pcsaft_hres(T, rho_liq, x, pyargs)
    S_res_liq = pcsaft.pcsaft_sres(T, rho_liq, x, pyargs)

    # 如果存在聚合物, 需要进行换算
    if has_polymer:
        x[polymer_index] = x[polymer_index]*DPN*mw_seg_avg/mw[polymer_index]
        H_res_liq = H_res_liq/sum(x)
        S_res_liq = S_res_liq/sum(x)
        beta = beta/(beta+(1-beta)*sum(x)) 

    H_mix = H_ig_mix + beta*H_res_vap + (1-beta)*H_res_liq
    S_mix = S_ig_mix + beta*S_res_vap + (1-beta)*S_res_liq   
      
    return beta, H_mix, S_mix



def properties_of_stream(components, args, stream, method="PC-SAFT"):
    """
    计算混合物的性质, 需要先判断流股的相态: 如果是单相, 可直接计算
    如果是气液两相, 则需要进行闪蒸计算, 分成气液两相进行计算, 再加和

    Params
    ----------
    components : 组分配置
    args : 物性参数
    stream : 流股参数
    method : 物性方法
    """

    T = stream["T"]   # 温度
    P = stream["P"]   # 压力

    NC = 0    # 组分数
    for component in stream["Mole Fraction"].keys():
        if stream["Mole Fraction"][component] != 0 and stream["Mole Fraction"][component] != None:
            NC += 1
    
    # 将各组分参数包装为np.array的形式, 以方便求和运算
    z = np.zeros(NC)      # 各组分摩尔分数
    m = np.zeros(NC)      # 各组分链段数, Aspen用PCSFTM表示
    s = np.zeros(NC)      # 各组分链段直径, Aspen用PCSFTV表示
    e = np.zeros(NC)      # 各组分能量参数, Aspen用PCSFTU表示
    e_assoc = np.zeros(NC)  # 关联成分的关联能量,单位K, Aspen用PCSFAU表示
    vol_a = np.zeros(NC)    # 关联成分的有效关联量, Aspen用PCSFAV表示
    w = np.zeros(NC)      # 各组分偏心因子
    Tc = np.zeros(NC)     # 各组分临界温度
    Pc = np.zeros(NC)     # 各组分临界压力
    mw = np.zeros(NC)     # 相对分子质量
    H_ig_comp = np.zeros(NC)  # 各组分的理想焓
    S_ig_comp = np.zeros(NC)  # 各组分的理想熵

    i = -1
    has_polymer = False
    for component in stream["Mole Fraction"].keys():
        if stream["Mole Fraction"][component] != 0 and stream["Mole Fraction"][component] != None:
            i = i + 1
            z[i] = stream["Mole Fraction"][component]
            w[i] = args[component]["w"]
            Tc[i] = args[component]["Tc"]
            Pc[i] = args[component]["Pc"]
            mw[i] = args[component]["mw"]
            if components[component]["type"] == "conventional":
                m[i] = args[component]["m"]
                s[i] = args[component]["s"]
                e[i] = args[component]["e"]
                if args[component].get("e_assoc") != None:
                    e_assoc[i] = args[component].get("e_assoc")
                else:
                    e_assoc[i] = 0
                if args[component].get("vol_a") != None:
                    vol_a[i] = args[component].get("vol_a")
                else:
                    vol_a[i] = 0
                eqno = args[component]["cp_mol_ig"]["eqno"]
                params = args[component]["cp_mol_ig"]["params"]
                enth_mol_form_ig_ref = args[component]["enth_mol_form_ig_ref"]
                entr_mol_form_ig_ref = args[component]["entr_mol_form_ig_ref"]
                H_ig_comp[i] = enth_ig_comp(T, enth_mol_form_ig_ref, eqno, params, Tref=298.15)
                S_ig_comp[i] = entr_ig_comp(T, P*z[i], entr_mol_form_ig_ref, eqno, params, Tref=298.15, Pref=101325)
            elif components[component]["type"] == "polymer":
                has_polymer = True
                polymer_index = i
                DPN = stream[component]["DPN"]

                Ns = 0   # 链段的种类数
                for segment in stream[component]["SFRAC"].keys():
                    if stream[component]["SFRAC"][segment] != 0 or stream[component]["SFRAC"][segment] != None:
                        Ns = Ns + 1

                sfrac = np.zeros(Ns)   # 各链段分数
                mw_seg = np.zeros(Ns)  # 各链段分子质量
                r_seg = np.zeros(Ns)   # 各链段的链段数/数均分子量
                s_seg = np.zeros(Ns)   # 各链段的链段直径
                e_seg = np.zeros(Ns)   # 各链段的能量参数
                e_assoc_seg = np.zeros(Ns)
                vol_a_seg = np.zeros(Ns)

                H_ig_seg = np.zeros(Ns)  # 各链段的理想焓
                S_ig_seg = np.zeros(Ns)  # 各链段的理想熵

                j = -1
                for segment in stream[component]["SFRAC"].keys():
                    if stream[component]["SFRAC"][segment] != 0 or stream[component]["SFRAC"][segment] != None:
                        j = j + 1
                        sfrac[j] = stream[component]["SFRAC"][segment]
                        mw_seg[j] = args[segment]["mw"]
                        r_seg[j] = args[segment]["r"]
                        s_seg[j] = args[segment]["s"]
                        e_seg[j] = args[segment]["e"]
                        if args[segment].get("e_assoc") != None:
                            e_assoc_seg[j] = args[segment].get("e_assoc")
                        else:
                            e_assoc_seg[j] = 0
                        if args[segment].get("vol_a") != None:
                            vol_a_seg[j] = args[segment].get("vol_a")
                        else:
                            vol_a_seg[j] = 0
                        eqno = args[segment]["cp_mol_ig"]["eqno"]
                        params = args[segment]["cp_mol_ig"]["params"]
                        enth_mol_form_ig_ref = args[segment]["enth_mol_form_ig_ref"]
                        entr_mol_form_ig_ref = args[segment]["entr_mol_form_ig_ref"]
                        H_ig_seg[j] = enth_ig_comp(T, enth_mol_form_ig_ref, eqno, params, Tref=298.15)
                        S_ig_seg[j] = entr_ig_comp(T, P*z[i]*sfrac[j], entr_mol_form_ig_ref, eqno, params, Tref=298.15, Pref=101325)

                mw_seg_avg = sum(sfrac*mw_seg) # 链段的平均相对分子质量
                MWN = sfrac*mw_seg*DPN
                m[i] = np.sum(r_seg*MWN) 
                s[i] = np.sum(sfrac*s_seg)
                e[i] = np.sum(sfrac*e_seg)
                e_assoc[i] = np.sum(sfrac*e_assoc_seg)
                vol_a[i] = np.sum(sfrac*vol_a_seg)

                H_ig_comp[i] = np.sum(sfrac*H_ig_seg)*mw[i]/mw_seg_avg  # 聚合物的理想焓(以链段表示)
                S_ig_comp[i] = np.sum(sfrac*S_ig_seg)*mw[i]/mw_seg_avg  # 聚合物的理想熵(以链段表示)

    H_ig_mix = np.sum(H_ig_comp*z)   # 混合物的理想焓
    S_ig_mix = np.sum(S_ig_comp*z)   # 混合物的理想熵

    pyargs = {"Tc": Tc, "Pc": Pc, "w": w, "m": m, "s": s, "e": e, 
            "e_assoc": e_assoc, "vol_a": vol_a}
    # 如果存在聚合物, 需要进行换算
    if has_polymer:
        z[polymer_index] = z[polymer_index]*mw[polymer_index]/mw_seg_avg/DPN
        z = z/sum(z)

        # 如果存在聚合物,则体系只会是液相或者气液
        
        Pb, _, state1 = bubble_pressure(T, z, pyargs)
        Tb, _, state2 = bubble_temperature(T, z, pyargs)

        # 大于泡点压力，处于液相
        if stream.get("Valid phases") == "liquid" or (state1 and P >= Pb) or (state2 and T <= Tb) :
            x = z
            y = z
            beta = 0
        else:
            beta, x, y, iter_total = VLE_TPflash(T, P, z, pyargs)
    else:
        # 有效相为液相
        if stream.get("Valid phases") == "liquid":
            x = z
            y = z
            beta = 0
        elif stream.get("Valid phases") == "vapor":
            x = z
            y = z
            beta = 1
        else:
            Pb, _, state1 = bubble_pressure(T, z, pyargs)
            Tb, _, state2 = bubble_temperature(P, z, pyargs)
            Pd, _, state3 = dew_pressure(T, z, pyargs)
            Td, _, state4 = dew_temperature(P, z, pyargs)

            # 考虑纯物质的情况(P=Pb=Pd)
            if abs(Pb-Pd)<0.1 and abs(P-Pb)<0.1:
                if (stream.get("Valid phases") == "vapor"):
                    beta = 1
                elif (stream.get("Valid phases") == "liquid"):
                    beta = 0
                else:
                    raise InputError("Please input valid phases of the stream")
                x = z
                y = z
            else:
                # 高于泡点压力，处于液相
                if (state1 and P >= Pb) or (state2 and T <= Tb):
                    x = z
                    y = z
                    beta = 0
                # 低于露点压力，处于气相
                elif (state3 and P <= Pd) or (state4 and T >= Td) :
                    x = z
                    y = z
                    beta = 1
                else:
                    beta, x, y, iter_total = VLE_TPflash(T, P, z, pyargs)
    
    rho_vap = pcsaft.pcsaft_den(T, P, y, pyargs, "vap")
    H_res_vap = pcsaft.pcsaft_hres(T, rho_vap, y, pyargs)
    S_res_vap = pcsaft.pcsaft_sres(T, rho_vap, y, pyargs)
    # 剩余熵还需减去压力导致的熵变,
    # S_res_vap = S_res_vap - R*np.sum(y*np.log(P*y/101325))
    for i in y:
        if i != 0:
            S_res_vap = S_res_vap - R*i*np.log(P*i/101325)

    rho_liq = pcsaft.pcsaft_den(T, P, x, pyargs, "liq")
    H_res_liq = pcsaft.pcsaft_hres(T, rho_liq, x, pyargs)
    S_res_liq = pcsaft.pcsaft_sres(T, rho_liq, x, pyargs)

    # 如果存在聚合物, 需要进行换算
    if has_polymer:
        x[polymer_index] = x[polymer_index]*DPN*mw_seg_avg/mw[polymer_index]
        H_res_liq = H_res_liq/sum(x)
        S_res_liq = S_res_liq/sum(x)
        beta = beta/(beta+(1-beta)*sum(x)) 

    H_mix = H_ig_mix + beta*H_res_vap + (1-beta)*H_res_liq
    S_mix = S_ig_mix + beta*S_res_vap + (1-beta)*S_res_liq   
      
    return beta, H_mix, S_mix


if __name__ == "__main__":
    # 设置参数
    # 设置组分参数：C3H6、C3H8、H2、N2、PP、TiCl4、TEA、C3H6-R
    components = {
        "C3H6": {"type": "conventional"},
        "C3H8": {"type": "conventional"},
        "H2": {"type": "conventional"},
        "N2": {"type": "conventional"},
        "PP": {"type": "polymer"},
        "TiCl4": {"type": "conventional"},
        "TEA": {"type": "conventional"},
        "C3H6-R": {"type": "segment"},
        "H2O": {"type": "conventional"},
    }

    # 设置催化剂
    catalyst = {"TiCl4": {"type": "Z-N", 
                        "site_types_num":4,   # 位点数
                        "site_conc":0.1,      # 位点浓度, mol/kgcat
                }}

    # POWDER1、POWDER2、GAS1是闪蒸罐STRIP1的进出口物流
    # 设置流股参数: 温度、压力、摩尔流量、摩尔分数、
    # 如果流股中有聚合物，则需要设置聚合物的数均聚合度DPN和链段分数SFRAC
    POWDER1 = {
        "T": 342.15, "P": 2800000,  "Mole Flow": 136.2576,  
        "Mole Fraction": { 
            "C3H6": 0.2098286, "C3H8": 0.00408862, "H2": 0.0005843, "N2": 0.000137371, 
            "PP": 0.7851512, "TiCl4": 3.22413E-05, "TEA": 0.000177694, "H2O": 0,
        },
        "PP": {
            "DPN": 535.4335, "SFRAC": {"C3H6-R":1}
        }
    }
    print("\nPOWDER1: ")
    print(molar_enthalpy_of_mixture(components, args, POWDER1))


    POWDER2 = {
        "T": 338.15, "P": 500000,  "Mole Flow": 109.8246,  
        "Mole Fraction": { 
            "C3H6": 0.0251184, "C3H8": 0.00049539, "H2": 1.04542E-06, "N2": 7.41444E-07, 
            "PP": 0.9741239, "TiCl4": 4.00012E-05, "TEA": 0.000220462, "H2O": 0,
        },
        "PP": {
            "DPN": 535.4335, "SFRAC": {"C3H6-R":1}
        }
    }
    print("\nPOWDER2: ")
    print(molar_enthalpy_of_mixture(components, args, POWDER2))

    GAS1 = {
        "T": 338.15, "P": 500000, "Mole Flow": 136.2576,  
        "Mole Fraction": {
            "C3H6": 0.9772706, "C3H8": 0.0190179, "H2": 0.00300639,"N2": 0.000705038, 
            "PP": 0, "TiCl4": 1.8894E-15, "TEA": 1.0413E-14, "H2O": 0,
        },
        "PP": {
            "DPN": 535.4335, "SFRAC": {"C3H6-R":1}
        }
    }
    print("\nGAS1: ")
    print(molar_enthalpy_of_mixture(components, args, GAS1))

    
    # 以下物流是混合器FMIX1的进出口物流
    # 设置流股参数: 温度、压力、摩尔流量、摩尔分数、
    # 如果流股中有聚合物，则需要设置聚合物的数均聚合度DPN和链段分数SFRAC
    C3FEED = {"T": 303.15, "P": 3000000, "Mole Flow": 134.1611,  
        "Mole Fraction": { "C3H6": 0.9958757, "C3H8": 0.00412434, "H2": 0,
            "N2": 0, "PP": 0, "TiCl4": 0, "TEA": 0, "H2O":0}
    }

    CAT = {"T": 303.15, "P": 3000000, "Mole Flow": 1.964734,  
        "Mole Fraction": { "C3H6": 0.9958403, "C3H8": 0.00192372, "H2": 0,
            "N2": 0, "PP": 0, "TiCl4": 0.00223599, "TEA": 0, "H2O":0},
        "Valid phases": "liquid"
    }

    COCAT = {"T": 303.15, "P": 3000000, "Mole Flow": 0.0243309,  
        "Mole Fraction": { "C3H6": 0, "C3H8": 0, "H2": 0,
            "N2": 0, "PP": 0, "TiCl4": 0, "TEA": 1, "H2O":0},
        "Valid phases": "liquid"
    }

    CYCGAGB = {"T": 334.9868, "P": 3050000, "Mole Flow": 5963.128,  
        "Mole Fraction": { "C3H6": 0.8822209, "C3H8": 0.0169148, "H2": 0.0923867,
            "N2": 0.0084775, "PP": 0, "TiCl4": 0, "TEA": 0, "H2O":0}
    }

    H2FEED = {"T": 303.15, "P": 3000000, "Mole Flow": 0.2648948,  
        "Mole Fraction": { "C3H6": 0, "C3H8": 0, "H2": 1,
            "N2": 0, "PP": 0, "TiCl4": 0, "TEA": 0, "H2O":0}
    }

    N2FEED = {"T": 303.15, "P": 3000000, "Mole Flow": 0.0187177,  
        "Mole Fraction": { "C3H6": 0, "C3H8": 0, "H2": 0,
            "N2": 1, "PP": 0, "TiCl4": 0, "TEA": 0, "H2O":0}
    }

    streams = [C3FEED, CAT, COCAT, CYCGAGB, H2FEED, N2FEED]

    print("\nFeed1: ")
    Feed1 = mix(components, catalyst, args, streams)
    print(Feed1)
    print(molar_enthalpy_of_mixture(components, args, Feed1)) 

    # 压缩机测试
    VAP_A = {"T": 342.15, "P": 2800000, "Mole Flow": 5963.128,  
        "Mole Fraction": { "C3H6": 0.8822209, "C3H8": 0.0169148, "H2": 0.0923867,
            "N2": 0.0084775, "PP": 0, "TiCl4": 2.1762E-14, "TEA": 1.1994E-13, "H2O":0}
    }
    print("\nVAP_A")
    print(molar_enthalpy_of_mixture(components, args, VAP_A)) 

    pressure_spec = {"Pressure increase": 250000}

    print("\nVAP-B: ")
    VAP_B, T_isentropic, isentropic_compression_work, true_compression_work, net_work = compress_isentropic(components, catalyst, args, [VAP_A], pressure_spec)
    print(VAP_B)
    print("T_isentropic: ", T_isentropic)
    print(isentropic_compression_work, true_compression_work, net_work)

    print(molar_enthalpy_of_mixture(components, args, VAP_B)) 

    # 加热器测试
    S1 = {"T": 473.15, "P": 100000, "Mole Flow": 0.2777778,  
        "Mole Fraction": { "C3H6": 1, "C3H8": 0}
    }

    S2 = {"T": 453.15, "P": 100000, "Mole Flow": 0.2777778,  
        "Mole Fraction": { "C3H6": 0, "C3H8": 1}
    }

    streams = [S1, S2]
    flash_type = {"T":150+273.15, "P":100000}

    print("S3")
    S3, duty = heat(components, catalyst, args, streams, flash_type)
    print(S3)
    print(duty)

    # 多流股换热器的测试
    S4 = {"T": 120+273.15, "P": 100000, "Mole Flow": 0.2777778,  
        "Mole Fraction": { "C3H6": 1, "C3H8": 0, "H2O":0}
    }

    S5 = {"T": 150+273.15, "P": 100000, "Mole Flow": 0.2777778,  
        "Mole Fraction": { "C3H6": 0, "C3H8": 1 , "H2O":0}
    }

    hot_streams = [S4, S5]

    S6 = {"T": 20+273.15, "P": 100000, "Mole Flow": 0.2777778,  
        "Mole Fraction": { "C3H6": 0, "C3H8": 0, "H2O":1}
    }

    S7 = {"T": 25+273.15, "P": 100000, "Mole Flow": 0.2777778,  
        "Mole Fraction": { "C3H6": 0, "C3H8": 0, "H2O":1}
    }
    
    cold_streams = [S6, S7]

    specs = {
        "hot":[{"spec":{"T":100+273.15}, "P":100000, "Valid phases": "vapor-liquid", "Max iter":100, "Tol":0.01},
               {"spec":{"T":120+273.15}, "P":100000, "Valid phases": "vapor-liquid", "Max iter":100, "Tol":0.01}],
        "cold": [{"spec":{"T":320}, "P":0, "Valid phases": "vapor-liquid", "Max iter":100, "Tol":0.01},
                 {"spec":None, "P":0, "Valid phases": "vapor-liquid", "Max iter":100, "Tol":0.01}]
    }
    print("\n多流股换热器1: ")
    hot_out, hot_duty, cold_out, cold_duty, duty = multi_stream_heat_exchang(components,catalyst, args, hot_streams, cold_streams, specs)
    print(hot_out)
    print(hot_duty)
    print(cold_out)
    print(cold_duty)
    print(duty)

    print("\n多流股换热器2: ")
    specs = {
        "hot":[{"spec":None, "P":100000, "Valid phases": "vapor-liquid", "Max iter":100, "Tol":0.01},
               {"spec":None, "P":100000, "Valid phases": "vapor-liquid", "Max iter":100, "Tol":0.01}],
        "cold": [{"spec":{"T":340}, "P":0, "Valid phases": "vapor-liquid", "Max iter":100, "Tol":0.01},
                 {"spec":{"T":330}, "P":0, "Valid phases": "vapor-liquid", "Max iter":100, "Tol":0.01}]
    }
    hot_out, hot_duty, cold_out, cold_duty, duty = multi_stream_heat_exchang(components, catalyst, args, hot_streams, cold_streams, specs)
    print(hot_out)
    print(hot_duty)
    print(cold_out)
    print(cold_duty)
    print(duty)


    C3FEED = {"T": 303.15, "P": 3000000, "Mole Flow": 2004.329,  
        "Mole Fraction": { "C3H6": 0.8331688, "C3H8": 0.1668312, "H2": 0,
            "N2": 0, "PP": 0, "TiCl4": 0, "TEA": 0, "H2O":0}
    }
    print("Mole entropy")
    print(molar_enthalpy_of_mixture(components, args, C3FEED)[2])
    

    

    





    


    