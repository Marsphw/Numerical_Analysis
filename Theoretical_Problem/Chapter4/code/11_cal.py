import numpy as np

p = [20, 30, 40]
for i in range(len(p)):
    prod = 1
    mul = 1
    for j in range(1, p[i] + 1):
        prod *= (1 + j * p[i])
        mul *= j
    result = (prod - np.power(p[i], p[i])) / mul
    print("p =", p[i], " result =", result)
