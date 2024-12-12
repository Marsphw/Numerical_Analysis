import numpy as np
import matplotlib.pyplot as plt

epsilon = 2**(-24)
x = np.linspace(0.001, 1, 1000)
y1 = abs(x*np.exp(-x)/(1 - np.exp(-x)))
y2 = abs(1 / epsilon * np.log(1 + epsilon - epsilon*np.exp(x)) / x)

plt.plot(x, y1, label='condition number of function', color='red', linestyle = '--')
plt.plot(x, y2, label='condition number of algorithm', color='blue')
title = 'Condition number of function and algorithm'
plt.title(title)
plt.xlabel('x')
plt.ylabel('Condition number')
plt.legend()
plt.savefig('../images/9_draw.png')
plt.show()