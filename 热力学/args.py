"""
作为聚丙烯模拟的物性包
"""

import numpy as np

# 各组分参数：链段数m、链段直径s、能量参数e、偏心因子w、临界温度Tc、临界压力Pc、相对分子质量mw
# 计算理想气体摩尔比热容cp_mol_ig所需的公式形式eqno和参数params
# 参考温度下理想气体摩尔生成焓enth_mol_form_ig_ref, 理想气体摩尔生成熵entr_mol_form_ig_ref, 
# 组分："C3H6", "C3H8", "H2", "N2", "PP", "TICL4", "TEA" 
args = {
    "C2H4": {
        "m": 1.593,   
        "s": 3.445,
        "e": 176.47,
        "w": 0.0862484,
        "Tc": 282.34,
        "Pc": 5041000,
        "mw": 28.05376,
        "cp_mol_ig": {"eqno": None,    # 单位: J/mol/K
                "params": None },
        "enth_mol_form_ig_ref": None,     # 单位: J/mol
        "entr_mol_form_ig_ref": None,  # 单位: J/mol/K
    },

    "C3H6": {
        "m": 1.9597,   
        "s": 3.5356,
        "e": 207.19,
        "w": 0.137588,
        "Tc": 364.85,
        "Pc": 4600000,
        "mw": 42.08064,
        "cp_mol_ig": {"eqno": 107,    # 单位: J/mol/K
                "params": np.array([43.852, 150.600, 1398.8, 74.754, 616.46, 130, 1500])},
        "enth_mol_form_ig_ref": 20230,     # 单位: J/mol
        "entr_mol_form_ig_ref": -142.2438369948013,  # 单位: J/mol/K
    },

    "C3H8": {
        "m": 2.002,   
        "s": 3.6184,
        "e": 208.11,
        "w": 0.152291,
        "Tc": 369.83,
        "Pc": 4248000,
        "mw": 44.09652,
        "cp_mol_ig": {"eqno": 107,
                "params": np.array([59.474, 126.610, 844.31, 86.165, 2482.7, 298.15, 1500])},
        "enth_mol_form_ig_ref": -104680,     # 单位: J/mol
        "entr_mol_form_ig_ref": -269.29397954049974,  # 单位: J/mol/K
    },
    # # H2数据来源于POLYPCSAFT数据库
    # "H2": {
    #     "m": 0.8285,   
    #     "s": 2.9729,
    #     "e": 12.53,
    #     "w": -0.215993,
    #     "Tc": 33.19,
    #     "Pc": 1313000,
    #     "mw": 2.01588,
    #     "cp_mol_ig": {"eqno": 107,
    #             "params": np.array([27.617, 9.560, 2466, 3.760, 567.6, 250, 1500])},
    #     "enth_mol_form_ig_ref": 0,     # 单位: J/mol
    #     "entr_mol_form_ig_ref": 0,  # 单位: J/mol/K
    # },

    # H2数据来源于POLYPCSAFT数据库
    "H2": {
        "m": 0.8285,   
        "s": 2.9729,
        "e": 12.53,
        "w": -0.215993,
        "Tc": 33.19,
        "Pc": 1313000,
        "mw": 2.01588,
        "cp_mol_ig": {"eqno": 107,
                "params": np.array([27.617, 9.560, 2466, 3.760, 567.6, 250, 1500])},
        "enth_mol_form_ig_ref": 0,     # 单位: J/mol
        "entr_mol_form_ig_ref": 0,  # 单位: J/mol/K
    },

    "N2": {
        "m": 1.2053,   
        "s": 3.313,
        "e": 90.96,
        "w": 0.0377215,
        "Tc": 126.2,
        "Pc": 3400000,
        "mw": 28.01348,
        "cp_mol_ig": {"eqno": 107,    # 单位: J/mol/K
                "params": np.array([29.105, 8.6149, 1701.6, 0.10347, 909.79, 50, 1500])},
        "enth_mol_form_ig_ref": 0,     # 单位: J/mol
        "entr_mol_form_ig_ref": 0,  # 单位: J/mol/K
    },

    "PP": {
        "m": None,   
        "s": None,
        "e": None,
        "w": 0,
        "Tc": 2000,
        "Pc": 5000000,
        "mw": 42.08064,
        "cp_mol_ig": None,
        "enth_mol_form_ig_ref": None,  # 单位: J/mol
        "entr_mol_form_ig_ref": None,  # 单位: J/mol/K
    },

    "HDPE": {
        "m": None,   
        "s": None,
        "e": None,
        "w": 0,
        "Tc": 2000,
        "Pc": 5000000,
        "mw": 28.05376,
        "cp_mol_ig": {"eqno": None,    # 单位: J/mol/K
                "params": None },
        "enth_mol_form_ig_ref": None,     # 单位: J/mol
        "entr_mol_form_ig_ref": None,  # 单位: J/mol/K
    },



    # 数据来源于案例Aspen文件：pp
    "TiCl4": {
        "m": 20,   
        "s": 3.5356,
        "e": 207.19,
        "w": 0.283732,
        "Tc": 638,
        "Pc": 4660950,
        "mw": 189.6908,
        "cp_mol_ig": {"eqno": 107,    # 单位: J/mol/K
            "params": np.array([53.717, 54.368, 342.05, 28.474, 159.71, 100, 1500])},
        "enth_mol_form_ig_ref": -761660,     # 单位: J/mol
        "entr_mol_form_ig_ref": -117.25505803719125,  # 单位: J/mol/K
    },

    # 数据来源于案例Aspen文件：pp
    "TEA": {
        "m": 20,   
        "s": 3.5356,
        "e": 207.19,
        "w": 0.841783,
        "Tc": 678.15,
        "Pc": 8930000,
        "mw": 114.16664,
        "cp_mol_ig": {"eqno": 107,    # 单位: J/mol/K
            "params": np.array([112.060, 384.330, 1533, 242.740, 707.5, 298.15, 1500])},
        "enth_mol_form_ig_ref": -163600,     # 单位: J/mol
        "entr_mol_form_ig_ref": -548.7115250378416,  # 单位: J/mol/K
    },

    "C3H6-R": {
        "r": 0.02528,   
        "s": 4.1473,
        "e": 298.6392,
        "w": None,
        "Tc": None,
        "Pc": None,
        "mw": 42.08064,
        "cp_mol_ig": {"eqno": 0,
            "params": np.array([-42.339, 0.50092, -0.0005574, 0.0000002412, 0, 0, 280, 1000, 36.0292, 0.00000161263, 2.927166])},
        "enth_mol_form_ig_ref": -70700,     # 单位: J/mol
        "entr_mol_form_ig_ref": -317,  # 单位: J/mol/K
    },

    "C2H4-R": {
        "r": 0.04132,   
        "s": 3.4751,
        "e": 267.1854,
        "w": None,
        "Tc": None,
        "Pc": None,
        "mw": 28.05376,
        "cp_mol_ig": {"eqno": None,    # 单位: J/mol/K
                "params": None },
        "enth_mol_form_ig_ref": None,     # 单位: J/mol
        "entr_mol_form_ig_ref": None,  # 单位: J/mol/K
    },


    "Hexane": {
        "m": 3.0576,   
        "s": 3.7983,
        "e": 236.77,
        "w": 0.301261,
        "Tc": 507.6,
        "Pc": 3025000,
        "mw": 86.17716,
        "cp_mol_ig": {"eqno": None,    # 单位: J/mol/K
                "params": None },
        "enth_mol_form_ig_ref": None,     # 单位: J/mol
        "entr_mol_form_ig_ref": None,  # 单位: J/mol/K
    },

    "H2O": {
        "m": 1.0656,   
        "s": 3.0007,
        "e": 366.51,
        "e_assoc": 2500.7,
        "vol_a": 0.034868,     
        "w": 0.344861,
        "Tc": 647.096,
        "Pc": 22064000,
        "mw": 18.01528,
        "cp_mol_ig": {"eqno": 107,    # 单位: J/mol/K
            "params": np.array([33.363, 26.790, 2610.5, 8.896, 1169, 100, 2273.15])},
        "enth_mol_form_ig_ref": -241814,     # 单位: J/mol
        "entr_mol_form_ig_ref": -44.353513332215334,  # 单位: J/mol/K
    }
}