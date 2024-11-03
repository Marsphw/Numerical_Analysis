import numpy as np
import matplotlib.pyplot as plt

n = int(input())
n = n * 101

x = []
y = []

for i in range(n):
    x_temp = float(input())
    y_temp = float(input())
    x.append(x_temp)
    y.append(y_temp)
    
plt.figure()
plt.plot(x, y)
plt.title("Bezier Curve for m = " + str(n//101))
plt.xlabel('x')
plt.ylabel('y')
plt.grid(True)
plt.savefig("../results/Task_F/Task_F_m_"+str(n//101)+".png")