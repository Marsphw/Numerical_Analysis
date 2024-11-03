import numpy as np
import matplotlib.pyplot as plt
import math

for i in range(4):
    n = int(input())
    x = np.linspace(-1, 1, 1000)
    y2 = 1 / (1 + 25*pow(x, 2))
    y1 = []
    for j in range(1000):
        y_temp = float(input())
        y1.append(y_temp)
        
    plt.figure()
    plt.plot(x, y1, label = 'Interpolation Points')
    plt.plot(x, y2, label = 'True Function')
    for j in range(1, n + 1):
        plt.scatter([math.cos((2*j - 1)*math.pi/(2*n))], [1 / (1 + 25*pow(math.cos((2*j - 1)*math.pi/(2*n)), 2))], color = 'red', zorder = 5)
    plt.legend()
    plt.title('Interpolation Plot for n = '+str(n))
    plt.xlabel('x')
    plt.ylabel('y')
    plt.savefig('../results/Task_C/Interpolation_Plot_n_'+str(n)+'.png')
    plt.close()
    
    plt.figure()
    plt.plot(x, np.abs(y1 - y2), linestyle = '--', color = 'orange', label = 'Absolute Error')
    plt.legend()
    plt.title('Absolute Error Plot for n = '+str(n))
    plt.xlabel('x')
    plt.ylabel('Absolute Error')
    plt.savefig('../results/Task_C/Absolute_Error_Plot_n_'+str(n)+'.png')
    plt.close()