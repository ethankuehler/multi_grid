import numpy as np
import matplotlib.pyplot as pp


params = np.loadtxt('cmake-build-release/params.txt')
N = int(params[0])
L = params[1]
dx = params[2]
R = params[3]
dens = params[4]
a = params[5]
b = params[6]
c = params[7]
N2 = int(params[8])

