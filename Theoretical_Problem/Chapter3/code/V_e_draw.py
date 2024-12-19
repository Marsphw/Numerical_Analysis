import numpy as np
from matplotlib.pyplot import *

x1 = np.linspace(-1, 0, 3000)
y1 = ((x1 + 1)**2 / 2)
x2 = np.linspace(0, 1, 3000)
y2 = ((x2 + 1)*(1 - x2)/2 + (2 - x2)*x2/2)
x3 = np.linspace(1, 2, 3000)
y3 = ((2 - x3)**2 / 2)

figure
plot(x1, y1, 'r')
plot(x2, y2, 'r')
plot(x3, y3, 'r')
xlabel('x')
ylabel('B(x)')
title('B(x) for different x')
savefig('../images/B(x).png')
