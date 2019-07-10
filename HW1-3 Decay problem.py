# coding=UTF-8
# 此程序用来求衰变规律
import math
import matplotlib.pyplot as plt


#Na
i=0
for step in [0.4, 0.2, 0.1, 0.05]:  # 有四种时间步长，选定Taub=10, Taua=1
    Na = 1
    Nb = 1
    t = 0
    a = [0]
    b = [1]  # 添加初始值
    c = [1]
    while t <= 5:  # 添加元素
        Na = Na-step*Na
        Nb = Nb + (Na - Nb/10) * step
        t = t + step
        a.append(t)
        b.append(Na)
        c.append(Nb)
    plt.figure(1)
    i += 1
    if i == 1:  # 对于四种不同的步长使用不同标签标出
        plt.plot(a, b, label='step=0.4, Na')
    if i == 2:
        plt.plot(a, b, label='step=0.2, Na')
    if i == 3:
        plt.plot(a, b, label='step=0.1, Na')
    if i == 4:
        plt.plot(a, b, label='step=0.05, Na')

plt.figure(1)
x = [0]
for i in range(1, 101):  # 求出准确结果
    x.append(i * 0.05)
f2 = [1]
for i in range(1, 101):
    f2.append(1 * math.exp(-x[i]))


plt.plot(x, f2, label='accurate Na')
plt.title('decay curve')
plt.xlabel('Time/s')
plt.ylabel('N')
plt.legend()
plt.show()


#Nb
i = 0

for step in [0.4, 0.2, 0.1, 0.05]:  # 有四种时间步长，选定Taub=10, Taua=1
    Na = 1
    Nb = 1
    t = 0
    a = [0]
    b = [1]  # 添加初始值
    c = [1]
    while t <= 5:  # 添加元素
        Na = Na - step * Na
        Nb = Nb + (Na - Nb / 10) * step
        t = t + step
        a.append(t)
        b.append(Na)
        c.append(Nb)
    plt.figure(2)
    i += 1
    if i == 1:  # 对于四种不同的步长使用不同标签标出
        plt.plot(a, c, label='step=0.4, Nb')
    if i == 2:
        plt.plot(a, c, label='step=0.2, Nb')
    if i == 3:
        plt.plot(a, c, label='step=0.1, Nb')
    if i == 4:
        plt.plot(a, c, label='step=0.05, Nb')

plt.figure(2)
x = [0]
for i in range(1, 101):  # 求出准确结果
    x.append(i * 0.05)
f3 = [1]
for i in range(1, 101):
    f3.append((-10 / 9) * (math.exp(-x[i]) - math.exp((-x[i]) / 10)) + math.exp((-x[i]) / 10))

plt.plot(x, f3, label='accurate Nb')
plt.title('decay curve')
plt.xlabel('Time/s')
plt.ylabel('N')
plt.legend()
plt.show()


#Nb
k=0

for taub in [0.1, 1, 10]: #对不同的衰变特征时间, Taua=1
    Na = 1
    Nb = 1
    t = 0
    a = [0]
    b = [1]  # 添加初始值
    c = [1]
    step = 0.05
    while t <= 5:  # 添加元素
        Na = Na - step * Na
        Nb = Nb + (Na - Nb / taub) * step
        t = t + step
        a.append(t)
        c.append(Nb)
    plt.figure(3)
    k += 1
    if k==1:
        plt.plot(a, c, label='taub=0.1, Nb')
    if k==2:
        plt.plot(a, c, label='taub= 1, Nb')
    if k==3:
        plt.plot(a, c, label='taub= 10, Nb')

plt.figure(3)
x = [0]
for i in range(1, 101):  # 求出准确结果
    x.append(i * 0.05)
f3 = [1]
f4 = [1]
f5 = [1]
for i in range(1, 101):
    f3.append((1/9) * (math.exp(-x[i])-math.exp((-x[i])*10)) + math.exp((-x[i])*10))
plt.plot(x, f3, label='accurate Nb, taub = 0.1')
for i in range(1, 101):
    f4.append((1 + x[i])* math.exp(-x[i]))
plt.plot(x, f4, label='accurate Nb , taub = 1')
for i in range(1, 101):
    f5.append((-10/ 9) * (math.exp(-x[i]) - math.exp((-x[i]) / 10)) + math.exp((-x[i]) / 10))
plt.plot(x, f5, label='accurate Nb , taub = 10')

plt.title('decay curve')
plt.xlabel('Time/s')
plt.ylabel('N')
plt.legend()
plt.show()
