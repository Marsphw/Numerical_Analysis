import numpy as np
import matplotlib.pyplot as plt

N = int(input())
for j in range(N):

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
        plt.scatter([x_temp], [y_temp], color = 'red', zorder = 2, s = 10)

    x = np.linspace(x1, xN, 1000)
    y1 = []
    y2 = 1 / (1 + 25 * x ** 2)
    for i in range(1000):
        y1.append(float(input()))
    plt.plot(x, y1, label = 'The spline')
    plt.plot(x, y2, label = 'The true function')
    plt.legend()
    plt.title(filename)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.savefig("../figure/" + filename + ".png")