import numpy as np
import matplotlib.pyplot as plt

filename = input()

plt.figure()
n = int(input())
x1 = 0
xN = 0
for i in range(n):
    x_temp = float(input())
    y_temp = float(input())
    if i == 0:
        x1 = x_temp 
    if i == n-1:
        xN = x_temp
    plt.scatter([x_temp], [y_temp], color = 'red', zorder = 5)

x = np.linspace(x1, xN, 1000)
y = []
for i in range(1000):
    y.append(float(input()))
plt.plot(x, y, label = 'The spline')
plt.legend()
plt.title(filename)
plt.xlabel('x')
plt.ylabel('y')
plt.savefig("../figure/" + filename + ".png")