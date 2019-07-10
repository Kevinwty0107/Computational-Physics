import math
import matplotlib.pyplot as plt
from copy import deepcopy


def Q(k,q,n):
    res = 0
    for i in range(n+2):
        res += math.sin(math.pi*k*i/(n+1))*q[i]
    res = res * math.sqrt(2 / n)
    return res


def Q_dot(p, k, n):
    res = 0
    for i in range(n + 2):
        res += math.sin(math.pi * k * i / (n + 1)) * p[i]
    res = res * math.sqrt(2 / n)
    return res


def E(n, k, q, p):  # 求能量
    E_1 = 0.5 * Q_dot(p, k, n) ** 2
    E_2 = 0.5 * Q(k, q, n) ** 2 * (2 * math.sin(math.pi * k / (2*(n+1)))) ** 2
    return E_1 + E_2


n = 32
p = [0] * (n + 2)
Q_0 = [20] + [0] * (n-1)
q = []
for t in range(n+2):
    q.append(20*math.sin(t*math.pi/(n+1))/math.sqrt(n/2) *n/(n+1) )
q_real = q  # 真实的q
p_origin = p  # 记录初始p
E_list = [[] for i in range(32)]
for i in range(32):
    E_list[i].append(E(n, i+1, q_real, p))  # 添加t=0时能量

length = 4000
t_list = [i for i in range(length+1)]
resolution = 400000
dt = 2*math.pi/(2*math.sin(math.pi/2/(n+1)))*length/resolution
ddt = dt / 1000  # 第一次计算的时间分辨率
for i in range(500):  # 初始化计算
    q_old = q
    for i in range(n):
        q[i + 1] += ddt * p[i + 1]
    for i in range(n):  # 计算p和q的变化,假设得到的q是k+0.5的
        p[i + 1] += - ddt * (q_old[i + 1] - q_old[i + 2] + 0.25 * (q_old[i + 1] - q_old[i + 2]) ** 2
                            - q_old[i] + q_old[i + 1] - 0.25 * (q_old[i] - q_old[i + 1]) ** 2)
for j in range(resolution):
    print('完成', j/resolution*100, "%")
    q_old = q
    p = p_origin
    for i in range(n):  # 计算p和q的变化,得到的q是j+0.5的
        p[i + 1] += - dt * (q_old[i + 1] - q_old[i + 2] + 0.25 * (q_old[i + 1] - q_old[i + 2]) ** 2
                               - q_old[i] + q_old[i + 1] - 0.25 * (q_old[i] - q_old[i + 1]) ** 2)
    for i in range(n):
        q[i + 1] += dt * p[i + 1]
    for u in range(len(q)):  # 插值计算q
        q_real[u] = (q[u]+q_old[u])/2
    if (j+1)%(resolution/length) == 0:  # 在一定的位置再计算E
        for y in range(32):
            E_list[y].append(E(n, y+1, q_real, p))
plt.figure(1)
plt.plot(t_list, E_list[0], label='E1')
plt.plot(t_list, E_list[1], label='E2')
plt.plot(t_list, E_list[2], label='E3')
plt.plot(t_list, E_list[3], label='E4')
plt.xlabel('t')
plt.ylabel('E')
plt.title('Chaos for a long period of time')
plt.legend()
plt.show()

E_ave = [[] for i in range(32)]

for i in range(4):
    for j in range(len(E_list[i])):
        E_ave[i].append(0)
        for k in range(0,j):
            E_ave[i][j] += E_list[i][k]


for i in range(4):
    E_ave[i][0] = E_list[i][0]
    for j in range(1,len(E_list[i])):
        E_ave[i][j] = E_ave[i][j]/(j)


plt.figure(2)
plt.plot(t_list, E_ave[0], label='E1')
plt.plot(t_list, E_ave[1], label='E2')
plt.plot(t_list, E_ave[2], label='E3')
plt.plot(t_list, E_ave[3], label='E4')
plt.legend()
plt.xlabel('t')
plt.ylabel('average engergy')
plt.show()

E_ave = [0 for i in range(32)]
for i in range(32):
    for j in range(len(E_list[i])):
        E_ave[i] += E_list[i][j]
for i in range(32):
    E_ave[i] = E_ave[i]/(len(E_list[i]))
k_list = [i for i in range(1, 33)]
plt.plot(k_list, E_ave)
plt.xlabel('k')
plt.ylabel('E_k')
plt.title('average energy of different modes')
plt.show()