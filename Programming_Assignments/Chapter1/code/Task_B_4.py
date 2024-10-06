import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(0.1, 0.15, 1000)
y = (x**3 + 4*x**2 + 3*x + 5)/(2*x**3 - 9*x**2 + 18*x - 2)

y = np.where((2*x**3 - 9*x**2 + 18*x - 2) == 0, np.nan, y)

plt.plot(x, y)
plt.axhline(0, color='black',linewidth = 0.5)
plt.xlabel("x")
plt.ylabel("y")
plt.title("Task B.4")
plt.savefig('../images/task_b4.png')