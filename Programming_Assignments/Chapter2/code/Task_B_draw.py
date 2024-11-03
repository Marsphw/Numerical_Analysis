import numpy as np
import matplotlib.pyplot as plt

for i in range(4):
    n = int(input())
    x = np.linspace(-5, 5, 1000)
    y2 = 1 / (1 + pow(x, 2))
    y1 = []
    for j in range(1000):
        y_temp = float(input())
        y1.append(y_temp)
    plt.figure()
    plt.plot(x, y1, label = 'Interpolation Points')
    plt.plot(x, y2, label = 'True Function')
    for j in range(n + 1):
        plt.scatter([-5 + 10*j/n], [1 / (1 + pow(-5 + 10*j/n, 2))], color = 'red', zorder = 5)
    plt.legend()
    plt.title('Interpolation Plot for n = '+str(n))
    plt.xlabel('x')
    plt.ylabel('y')
    plt.savefig('../results/Task_B/Interpolation_Plot_n_'+str(n)+'.png')