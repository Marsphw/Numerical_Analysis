result1 = 0
result2 = 0

num = [1, 1e-9, 9e-10, 8e-10, 7e-10, 6e-10, 5e-10, 4e-10, 3e-10, 2e-10, 1e-10]

for i in range(0, len(num)):
    result1 += num[i]
    result2 += num[len(num)-1-i]

print(result1)
print(result2)
