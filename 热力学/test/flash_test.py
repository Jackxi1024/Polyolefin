# 闪蒸模块测试
# Aspen测试文件: two_pp_flash.apw

# 导包
from asyncio import streams
from polymer_model.thermodynamics.args import args
from polymer_model.thermodynamics.thermo import flash

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

F1 = {"T": 342.15, "P": 2800000, "Mole Flow": 0.1362576*1000,  
    "Mole Fraction": { "C3H6": 0.2098286, "C3H8": 0.00408862, "H2": 0.0005843,
        "N2": 0.000137371, "PP": 0.7851512, "TiCl4": 3.22413E-05, "TEA": 0.000177694},
    "PP":{"ZMOM": 0.000199806*1000, "FMOM": 0.1069828*1000, "DPN":535.4334,
          "SFRAC": {"C3H6-R":1}, "SFLOW":{"C3H6-R":0.1069828*1000}}
    }

F2 = {"T": 342.15, "P": 2800000, "Mole Flow": 0.1362576*1000,  
    "Mole Fraction": { "C3H6": 0.2, "C3H8": 0, "H2": 0,
        "N2": 0, "PP": 0.8, "TiCl4": 0, "TEA": 0},
    "PP":{"ZMOM": 0.000229352*1000, "FMOM": 0.109006*1000, "DPN":475.2779,
        "SFRAC": {"C3H6-R":1}, "SFLOW":{"C3H6-R":0.109006*1000}}
    }

streams = [F1, F2]
flash_type = {"T": 338.15  , "P": 500000}

vapor, liquid, beta, duty = flash(components, catalyst, args, streams, flash_type)
print("vapor:")
print(vapor)
print("liquid")
print(liquid)
print("beta", beta)
print("duty:", duty)