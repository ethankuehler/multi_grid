import numpy as np
import matplotlib.pyplot as pp

# getting params
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

m = (4.0 / 3.0) * np.pi * a * b * c * dens
# load in data

data = np.loadtxt('cmake-build-release/data.txt')
data = data.reshape(N, N, N)

plane = data[int(N / 2)]
line = plane[int(N / 2)]
del data
del plane

data2 = np.loadtxt('cmake-build-release/data2.txt')
data2 = data2.reshape(N, N, N)

plane2 = data2[int(N2 / 2)]
line2 = plane2[int(N2 / 2)]
del data2


data3 = np.loadtxt('cmake-build-release/data3.txt')
data3 = data3.reshape(N2, N2, N2)

plane3 = data3[int(N2 / 2)]
line3 = plane3[int(N2 / 2)]
del data3

# load in density
# dens_plot = np.loadtxt('cmake-build-release/f.txt')
# dens_plot = dens_plot.reshape(N, N, N)

# plot plane


x = np.linspace(-L / 2, L / 2, N)
x2 = np.linspace(-L / 2, L / 2, N2)
y = -m * (3 * R ** 2 - x ** 2) / (2 * R ** 3)
y[x < -a] = m / x[x < -a]
y[x > a] = -m / x[x > a]
# pp.plot(x, y, color='g')
pp.plot(x2, line, color='b')
pp.plot(x2, line2, color='r')
pp.plot(x2, line3, color='g')
pp.show()
'''
pp.pcolor(plane)
pp.show()
'''
pp.pcolor(plane2)
pp.show()
pp.pcolor(plane3)
pp.show()

# plot density
# dens_plane = dens_plot[int(N/2)]
# dens_line = dens_plane[int(N/2)]
# x = np.linspace(-L/2, L/2, N)
# pp.plot(x, dens_line)
# pp.show()
# pp.pcolor(dens_plane)
# pp.show()
