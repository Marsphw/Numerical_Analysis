import numpy as np
import matplotlib.pyplot as plt

Day = [0, 6, 10, 13, 17, 20, 28]
Sp = [[6.67, 17.3, 42.7, 37.3, 30.1, 29.3, 28.7],
      [6.67, 16.1, 18.9, 15.0, 10.6, 9.44, 8.89]]
x = np.linspace(0, 28, 1000)

for i in range(2):
    y1 = []
    for j in range(1000):
        y_temp = float(input())
        y1.append(y_temp)
    plt.figure()
    plt.plot(x, y1, label = 'Interpolation Points')
    for j in range(7):
        plt.plot([Day[j]], [Sp[i][j]], 'ro', zorder = 5)
    plt.legend()
    plt.title('Interpolation Plot for Sp'+str(i + 1))
    plt.xlabel('x')
    plt.ylabel('y')
    plt.savefig('../results/Task_E/Interpolation_Plot_Sp'+str(i + 1)+'.png')
    plt.close()