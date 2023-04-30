# 聚合反应设置

# 设置Z-N反应物种
# 设置Ziegler-Natta反应物种
species = { "polymer": "PP",          # 聚合物
            "tdb segment": None,      # 终端双键(用于计算支链)
            "monomers": "C3H6",       # 单体
            "segments": "C3H6-R",     # 链段
            "precatalyst": None,      # 预催化剂
            "catalyst": "TiCl4",      # 催化剂
            "cocatalysts": "TEA",     # 助催化剂
            "solvents": None,         # 溶剂
            "transfer agent": None,   # 链转移剂
            "hydrogens": "H2",        # 氢
            "poisons": None,          # 毒物
            "elec don": None,         # 电子供体
            "byproduct": None,        # 副产物
        }

# 设置Z-N反应：反应类型, 催化位点, 组分1, 组分2, 前指因子, 活化能, 反应级数, 终端双键分数, 参考温度
r1 = [["Act-Spon", 1, "TiCl4", None, 0.0013, 3.19872e4, 1, None, 343.15],
    ["Act-Spon", 2, "TiCl4", None, 0.0013, 3.19872e4, 1, None, 343.15],
    ["Act-Spon", 3, "TiCl4", None, 0.0013, 3.19872e4, 1, None, 343.15],
    ["Act-Spon", 4, "TiCl4", None, 0.0013, 3.19872e4, 1, None, 343.15],
    ["Chain-Ini", 1, "C3H6", None, 108.85/1000, 3.0145e4, 1, None, 343.15],
    ["Chain-Ini", 2, "C3H6", None, 24.5/1000,   3.0145e4, 1, None, 343.15],
    ["Chain-Ini", 3, "C3H6", None, 170.8/1000,  3.0145e4, 1, None, 343.15],
    ["Chain-Ini", 4, "C3H6", None, 60.55/1000,  3.0145e4, 1, None, 343.15],
    ["Propagation", 1, "C3H6-R", "C3H6", 108.85/1000,  3.0145e4, 1, None, 343.15],
    ["Propagation", 2, "C3H6-R", "C3H6", 24.5/1000,  3.0145e4, 1, None, 343.15],
    ["Propagation", 3, "C3H6-R", "C3H6", 170.8/1000,  3.0145e4, 1, None, 343.15],
    ["Propagation", 4, "C3H6-R", "C3H6", 60.55/1000,  3.0145e4, 1, None, 343.15],
    ["Chat-Mon", 1, "C3H6-R", "C3H6", 0.012/1000,  5.2e4, 1, None, 343.15],
    ["Chat-Mon", 2, "C3H6-R", "C3H6", 0.012/1000,  5.2e4, 1, None, 343.15],
    ["Chat-Mon", 3, "C3H6-R", "C3H6", 0.012/1000,  5.2e4, 1, None, 343.15],
    ["Chat-Mon", 4, "C3H6-R", "C3H6", 0.012/1000,  5.2e4, 1, None, 343.15],
    ["Chat-Cocat", 1, "C3H6-R", "TEA", 0.12/1000,  5.02416e4, 1, None, 343.15],
    ["Chat-Cocat", 2, "C3H6-R", "TEA", 0.12/1000,  5.02416e4, 1, None, 343.15],
    ["Chat-Cocat", 3, "C3H6-R", "TEA", 0.12/1000,  5.02416e4, 1, None, 343.15],
    ["Chat-Cocat", 4, "C3H6-R", "TEA", 0.12/1000,  5.02416e4, 1, None, 343.15],
    ["Chat-H2", 1, "C3H6-R", "H2", 4.8/1000,  4.47988e4, 1, None, 343.15],
    ["Chat-H2", 2, "C3H6-R", "H2", 8.88/1000,  4.47988e4, 1, None, 343.15],
    ["Chat-H2", 3, "C3H6-R", "H2", 2.64/1000,  4.47988e4, 1, None, 343.15],
    ["Chat-H2", 4, "C3H6-R", "H2", 6.6/1000,  4.47988e4, 1, None, 343.15],
    ["Deact-Spon", 1, None, None, 0.001,  4.1868e3, 1, None, 343.15],
    ["Deact-Spon", 2, None, None, 0.001,  4.1868e3, 1, None, 343.15],
    ["Deact-Spon", 3, None, None, 0.001,  4.1868e3, 1, None, 343.15],
    ["Deact-Spon", 4, None, None, 0.001,  4.1868e3, 1, None, 343.15],
    ["Atact-Prop", 1, "C3H6-R", "C3H6", 8/1000,  3.0145e4, 1, None, 343.15],
    ["Atact-Prop", 2, "C3H6-R", "C3H6", 1/1000,  3.0145e4, 1, None, 343.15],
    ["Atact-Prop", 3, "C3H6-R", "C3H6", 0.1/1000,  3.0145e4, 1, None, 343.15],
    ["Atact-Prop", 4, "C3H6-R", "C3H6", 0.1/1000,  3.0145e4, 1, None, 343.15]]

