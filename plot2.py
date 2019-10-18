import numpy as np
import matplotlib.pyplot as pp


x = np.loadtxt('cmake-build-release/line.txt')
x3 = np.loadtxt('cmake-build-release/line3.txt')
pp.plot(x, 'b')
pp.plot(x3, 'g')
pp.show()

