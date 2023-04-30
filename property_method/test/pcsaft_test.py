import polymer_model.ploymer_eos.pcsaft_elec as pcsaft
import numpy as np

T = 303.15
P = 3e6
m = np.array([20])
s = np.array([3.5356])
e = np.array([207.19])
pyargs= {"m":np.array([20]), "s":np.array([3.5356]), "e":np.array([207.19])}
x = np.array([1.0])

rho = pcsaft.pcsaft_den(x, m, s, e, T, P, phase='liq')
print(rho)