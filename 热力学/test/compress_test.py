# 压缩计算测试文件

from polymer_model.thermodynamics.thermo import compress_isentropic
from polymer_model.thermodynamics.args import args

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

# 压缩机测试
S1 = {"T": 342.15, "P": 2800000, "Mole Flow": 0.9,  
    "Mole Fraction": { "C3H6": 0.9, "C3H8": 0.0, "H2": 0.0, "N2": 0.1}
}
S2 = {"T": 333.15, "P": 2026500, "Mole Flow": 0.1,  
    "Mole Fraction": { "C3H6": 0.0, "C3H8": 0.5, "H2": 0.5, "N2": 0.0}
}

feed = [S1, S2]

pressure_spec = {"Discharge pressure": 3090410}

print("\nproduct: ")
product, T_isentropic, isentropic_compression_work, true_compression_work, net_work = compress_isentropic(components, args, feed, pressure_spec)
print(product)
print("T_isentropic: ", T_isentropic)
print(isentropic_compression_work, true_compression_work, net_work)
