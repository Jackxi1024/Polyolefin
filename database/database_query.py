import json
import math

# 加载JSON文件
database_json = open('E:\AZFT\polymer_model\database\chemsep1.json').read()
# 将json格式转化为dict格式
database = json.loads(database_json)
# # 读取数据
# H_ig = float(database["Ethylene"]["HeatOfFormation"]["@value"])   
# G_ig = float(database["Ethylene"]["GibbsEnergyOfFormation"]["@value"])  
# Cp_eqno = float(database["Ethylene"]["RPPHeatCapacityCp"]["eqno"]["@value"])
# Cp_A = float(database["Ethylene"]["RPPHeatCapacityCp"]["A"]["@value"])
# Cp_B = float(database["Ethylene"]["RPPHeatCapacityCp"]["B"]["@value"])
# Cp_C = float(database["Ethylene"]["RPPHeatCapacityCp"]["C"]["@value"])
# Cp_D = float(database["Ethylene"]["RPPHeatCapacityCp"]["D"]["@value"])
# Cp_E = float(database["Ethylene"]["RPPHeatCapacityCp"]["E"]["@value"])
# Cp_Tmin = float(database["Ethylene"]["RPPHeatCapacityCp"]["Tmin"]["@value"])
# Cp_Tmax = float(database["Ethylene"]["RPPHeatCapacityCp"]["Tmax"]["@value"])

comps = ["Ethylene", "Hydrogen", "Nitrogen", "1-butene", "Isopentane", "N-hexane"]
for i in comps:
    H_ig = float(database[i]["HeatOfFormation"]["@value"])/1000
    G_ig = float(database[i]["GibbsEnergyOfFormation"]["@value"])/1000
    S_ig = (H_ig-G_ig)/298.15
    Tc = float(database[i]["CriticalTemperature"]["@value"])  
    Pc = float(database[i]["CriticalPressure"]["@value"])  
    w = float(database[i]["AcentricityFactor"]["@value"])
    Cp_eqno = float(database[i]["RPPHeatCapacityCp"]["eqno"]["@value"])
    Cp_A = float(database[i]["RPPHeatCapacityCp"]["A"]["@value"])/1000
    Cp_B = float(database[i]["RPPHeatCapacityCp"]["B"]["@value"])/1000
    Cp_C = float(database[i]["RPPHeatCapacityCp"]["C"]["@value"])/1000
    Cp_D = float(database[i]["RPPHeatCapacityCp"]["D"]["@value"])/1000
    Cp_E = float(database[i]["RPPHeatCapacityCp"]["E"]["@value"])/1000
    Cp_Tmin = float(database[i]["RPPHeatCapacityCp"]["Tmin"]["@value"])
    Cp_Tmax = float(database[i]["RPPHeatCapacityCp"]["Tmax"]["@value"])
    print(i)
    print("H_ig: ", H_ig)
    print("S_ig: ", S_ig)
    print("Tc: ", Tc)
    print("Pc: ", Pc)
    print("w: ", w)
    print("Cp_eqno: ", Cp_eqno)
    print("Cp_A: ", Cp_A)
    print("Cp_B: ", Cp_B)
    print("Cp_C: ", Cp_C)
    print("Cp_D: ", Cp_D)
    print("Cp_E: ", Cp_E)
    print()
    


# RPP_type = set()
# CP_type = set()
# no = {"eqno":{"@value":"no"}}
# for i in database:
#     type = database[i].get("RPPHeatCapacityCp",no)["eqno"]["@value"]
#     RPP_type.add(type)
#     if type == "100":
#         print("100",i)
#     if type == "no":
#         print("no",i)    
#     type = database[i]["IdealGasHeatCapacityCp"]["eqno"]["@value"]
#     CP_type = set(type)
# print(RPP_type)
# print(CP_type)

# print(50*"*")
# print("Propylene")
# enth_ig_form = float(database["Propylene"]["HeatOfFormation"]["@value"])
# gibbs_ig_form = float(database["Propylene"]["GibbsEnergyOfFormation"]["@value"])
# entr_ig_form = (enth_ig_form - gibbs_ig_form)/298.15
# print("数据库")
# print(enth_ig_form)
# print(gibbs_ig_form)
# print(entr_ig_form)

# print(50*"*")
# print("N-hexane")
# enth_ig_form = float(database["N-hexane"]["HeatOfFormation"]["@value"])
# gibbs_ig_form = float(database["N-hexane"]["GibbsEnergyOfFormation"]["@value"])
# entr_ig_form = (enth_ig_form - gibbs_ig_form)/298.15
# print("数据库")
# print(enth_ig_form)  
# print(gibbs_ig_form)
# print(entr_ig_form)
# print("计算")
# P = 20265.88
# R = 8.314
# print("enth_ig_form: ", -166940000)
# print("entr_ig_form: ", -546320 + R*math.log(P/101325)*1000)

# print(50*"*")
# print("TiCL4")
# P = 1653.19
# R = 8.314
# print("enth_ig_form: ", -761660000)
# print("entr_ig_form: ", -83037.74 + R*math.log(P/101325)*1000)

# print(50*"*")
# print("TEA")
# P = 7.308071
# R = 8.314
# print("enth_ig_form: ", -163600000)
# print("entr_ig_form: ", -469420 + R*math.log(P/101325)*1000)

# print("\n")
# print(database.keys())

# print(80*"*")
# print(80*"*")
# comps = ["Ethylene", "Propylene", "Propane", "Hydrogen", "Nitrogen", "N-hexane", "Water"]
# for comp in comps:
#     eqno = float(database[comp]["VaporPressure"]["eqno"]["@value"]) 

#     A = float(database[comp]["VaporPressure"]["A"]["@value"])
#     B = float(database[comp]["VaporPressure"]["B"]["@value"])
#     C = float(database[comp]["VaporPressure"]["C"]["@value"])
#     D = float(database[comp]["VaporPressure"]["D"]["@value"])
#     E = float(database[comp]["VaporPressure"]["E"]["@value"])
#     Tmin = float(database[comp]["VaporPressure"]["Tmin"]["@value"])
#     Tmax = float(database[comp]["VaporPressure"]["Tmax"]["@value"])
#     print("\n", comp)
#     print("eqno: ", eqno)
#     print([A, B, C, D, E, Tmin, Tmax])
    

