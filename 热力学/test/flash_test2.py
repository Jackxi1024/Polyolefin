# 聚丙烯反应器出口气液平衡

# 导包
from polymer_model.thermodynamics.thermo import flash
from polymer_model.thermodynamics.test.args import args

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

F1 = {"T": 342.15, "P": 2800000, "Mole Flow": 6.099386*1000,  
    "Mole Fraction": { "C3H6": 0.8671999, "C3H8": 0.0166283, "H2": 0.0903359,
        "N2": 0.00829114, "PP": 0.0175399, "TiCl4": 7.20255E-07, "TEA": 3.9696E-06},
    "PP":{"ZMOM": 0.000199789*1000, "FMOM": 0.1069736*1000, "DPN":535.4335,
          "SFRAC": {"C3H6-R":1}, "SFLOW":{"C3H6-R":0.1069736*1000}}
    }


streams = [F1]
flash_type = {"T": 342.15, "P": 2800000}

vapor, liquid, beta, duty = flash(components, catalyst, args, streams, flash_type)
print("vapor:")
print(vapor)
print("liquid")
print(liquid)
print("beta", beta)
print("duty:", duty)