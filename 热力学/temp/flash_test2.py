from polymer_model.thermodynamics.flash_algorithm4 import *
import numpy as np

# # Flash test 2: "C3H6", "C3H8", "H2", "N2", "PP", "TICL4", "TEA"
# T = 338.15
# P = 500000
# PP_DPN = 535.4335 
# PP_MWN = 42.0806*PP_DPN
# z = np.array([0.2098286, 0.00408862, 0.0005843, 0.000137371, 0.7851512/PP_DPN,
#                 3.22413E-05, 0.000177694])
# # z = np.array([0.0251157, 0.000495333, 1.33715E-06, 7.41334E-07, 0.9741264/PP_DPN, 4.00013E-05, 0.000220463])

# z = z/np.sum(z)
# print("z", z)
# args = {"Tc": np.array([364.85, 369.83, 33.19, 126.2, 2000, 638, 678.15]), 
#         "Pc": np.array([4600000, 4248000, 1313000, 3400000, 5000000, 4660950, 8930000]), 
#         "w": np.array([0.137588, 0.152291, -0.215993, 0.0377215, 0, 0.283732, 0.841783]),
#         "m": np.array([1.9597, 2.002, 0.8285, 1.2053, 0.02528*PP_MWN, 20, 20]),
#         "s": np.array([3.5356, 3.6184, 2.9729, 3.313, 4.1473, 3.5356, 3.5356]),
#         "e": np.array([207.19, 208.11, 12.53, 90.96, 298.6392, 207.19, 207.19])}

# beta, x, y, iter_total = VLE_TPflash(T, P, z, args)
# print("Vapor fraction: ", beta)
# print("Liquid composition: ", x)
# x[4] = x[4]*PP_DPN
# beta = beta/(beta+(1-beta)*sum(x))
# print("新气相分数：", beta) 
# x = x/sum(x)
# print("新液相组成：", x)

# print("Vapor composition: ", y)
# print("iter_total: ", iter_total)
# # print("iter_inner: ", iter_inner)


# # Flash test 2: "C3H6", "C3H8", "H2", "N2", "PP", "TICL4", "TEA"
# T = 333.15
# P = 500000
# PP_DPN = 622.6436
# PP_MWN = 42.0806*PP_DPN
# z = np.array([0.2483033, 0.00487137, 0.000664964, 0.000158333, 0.7458073/PP_DPN,
#                 2.98812E-05, 0.000164897])
# # z = np.array([0.0251157, 0.000495333, 1.33715E-06, 7.41334E-07, 0.9741264/PP_DPN, 4.00013E-05, 0.000220463])

# z = z/np.sum(z)
# print("z", z)
# args = {"Tc": np.array([364.85, 369.83, 33.19, 126.2, 2000, 638, 678.15]), 
#         "Pc": np.array([4600000, 4248000, 1313000, 3400000, 5000000, 4660950, 8930000]), 
#         "w": np.array([0.137588, 0.152291, -0.215993, 0.0377215, 0, 0.283732, 0.841783]),
#         "m": np.array([1.9597, 2.002, 0.8285, 1.2053, 0.02528*PP_MWN, 20, 20]),
#         "s": np.array([3.5356, 3.6184, 2.9729, 3.313, 4.1473, 3.5356, 3.5356]),
#         "e": np.array([207.19, 208.11, 12.53, 90.96, 298.6392, 207.19, 207.19])}

# beta, x, y, iter_total = VLE_TPflash(T, P, z, args)
# print("Vapor fraction: ", beta)
# print("Liquid composition: ", x)
# x[4] = x[4]*PP_DPN
# beta = beta/(beta+(1-beta)*sum(x))
# print("新气相分数：", beta) 
# x = x/sum(x)
# print("新液相组成：", x)

# print("Vapor composition: ", y)
# print("iter_total: ", iter_total)





# # Pb = bubble_pressure(T, z, args)  # 泡点压力
# # print(Pb)

# x = np.array([0.0251184, 0.00049539, 1.04542E-06, 7.41444E-07, 0.9741239/PP_DPN, 4.00012E-05, 0.000220462])
# x = x/np.sum(x)
# print("x:",x) 

# y = np.array([0.9772695, 0.0190178, 0.00300763, 0.000705043, 0, 4.3295E-15, 2.3862E-14])
# print("y:",y)
# rho_liq = pcsaft.pcsaft_den(T, P, x, args, "liq")
# fug_liq = pcsaft.pcsaft_fugcoef(T, rho_liq, x, args)
# rho_vap = pcsaft.pcsaft_den(T, P, y, args, "vap")
# fug_vap = pcsaft.pcsaft_fugcoef(T, rho_vap, y, args)
# K = np.zeros(len(fug_liq))
# for i in range(len(fug_liq)):
#         if fug_vap[i] != 0:
#                 K[i] = fug_liq[i]/fug_vap[i]
# print("fug_liq", fug_liq)
# print("fug_vap", fug_vap)
# print("K:", K)
# print(np.sum(K*z))
# print(fug_liq*x-fug_vap*y)

# beta, x, y = Rachford_Rice_solver(K,z)
# print(beta)
# print(Rachford_Rice(0, z, K))
# print(np.sum(z*(K-1)/(1-beta+beta*K)))


# x = np.array([0.0251184, 0.00049539, 1.04542E-06, 7.41444E-07, 0.9741239/PP_DPN, 4.00012E-05, 0.000220462])
# x = x/np.sum(x)
# x[4] = x[4]*PP_DPN
# print(beta/(beta+(1-beta)*sum(x)))

# Flash test 1: "C3H6", "C3H8", "H2", "N2"
T = 341.15
P = 3000000
z = np.array([0.82, 0.16, 0.005, 0.015])
args = {"Tc": np.array([364.85, 369.83, 33.19, 126.2]), 
        "Pc": np.array([4600000, 4248000, 1313000, 3400000]), 
        "w": np.array([0.137588, 0.152291, -0.215993, 0.0377215]),
        "m": np.array([1.9597, 2.002, 0.9846, 1.2053]),
        "s": np.array([3.5356, 3.6184, 2.8263, 3.313]),
        "e": np.array([207.19, 208.11, 20.893, 90.96])}

beta, x, y, iter_total = VLE_TPflash(T, P, z, args)
print("Vapor fraction: ", beta)
print("Liquid composition: ", x)
print("Vapor composition: ", y)
print("iter_total: ", iter_total)

# Flash test 2: "C3H6", "C3H8", "H2", "N2", "PP", "TICL4", "TEA"
T = 223.15
P = 100000
PP_MWN = 119735
PP_DPN = 2845.37
z = np.array([0.7614716, 0.1435221, 0.0784851, 0.00720122, 0.00931595/PP_DPN, 6.25573E-07, 3.44657E-06])
z = z/sum(z)
args = {"Tc": np.array([364.85, 369.83, 33.19, 126.2, 2000, 638, 678.15]), 
        "Pc": np.array([4600000, 4248000, 1313000, 3400000, 5000000, 4660950, 8930000]), 
        "w": np.array([0.137588, 0.152291, -0.215993, 0.0377215, 0, 0.283732, 0.841783]),
        "m": np.array([1.9597, 2.002, 0.9846, 1.2053, 0.02528*PP_MWN, 20, 20]),
        "s": np.array([3.5356, 3.6184, 2.8263, 3.313, 4.1473, 3.54, 3.54]),
        "e": np.array([207.19, 208.11, 20.893, 90.96, 298.6392, 210, 210])}

beta, x, y, iter_total = VLE_TPflash(T, P, z, args)
print("Vapor fraction: ", beta)
print("Liquid composition: ", x)
x[4] = x[4]*PP_DPN
beta = beta/(beta+(1-beta)*sum(x))
print("新气相分数：", beta)
x = x/sum(x)
print("新液相组成：", x)

# Flash test 2: "C3H6", "C3H8", "H2", "N2", "PP", "TICL4", "TEA"
T = 333.15
P = 3000000
PP_MWN = 76386.77
PP_DPN = 1815.247
# z = np.array([0.7682945, 0.143522, 0.0784859, 0.00720122, 0.00249226/PP_DPN, 6.25573E-07, 3.46464E-06])
z = np.array([7.682945e-01, 1.435220e-01, 7.848590e-02, 7.201220e-03, 2.492260e-03/PP_DPN,
                6.255730e-07, 3.464640e-06])

z = z/np.sum(z)
print("z", z)
args = {"Tc": np.array([364.85, 369.83, 33.19, 126.2, 2000, 638, 678.15]), 
        "Pc": np.array([4600000, 4248000, 1313000, 3400000, 5000000, 4660950, 8930000]), 
        "w": np.array([0.137588, 0.152291, -0.215993, 0.0377215, 0, 0.283732, 0.841783]),
        "m": np.array([1.9597, 2.002, 0.9846, 1.2053, 0.02528*PP_MWN, 20, 20]),
        "s": np.array([3.5356, 3.6184, 2.8263, 3.313, 4.1473, 3.54, 3.54]),
        "e": np.array([207.19, 208.11, 20.893, 90.96, 298.6392, 210, 210])}

beta, x, y, iter_total = VLE_TPflash(T, P, z, args)
print("Vapor fraction: ", beta)
print("Liquid composition: ", x)
x[4] = x[4]*PP_DPN
beta = beta/(beta+(1-beta)*sum(x))
print("新气相分数：", beta)
x = x/sum(x)
print("新液相组成：", x)

print("Vapor composition: ", y)
print("iter_total: ", iter_total)