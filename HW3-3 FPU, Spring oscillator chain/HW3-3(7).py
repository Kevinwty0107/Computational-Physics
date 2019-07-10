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


n = 16
p = [0] * (n + 2)
q = []
Q_0 =[0]*10 + [1] + [0] * 5

for t in range(n+2):
    q.append(1*math.sin(11*t*math.pi/(n+1))*math.sqrt(2/n) *n/(n+1) )
print(Q(11,q,n))
print(Q(10,q,n))
print(Q(9,q,n))
q_real = q  # 真实的q
p_origin = p  # 记录初始p
E_1_list = []
E_2_list = []
E_3_list = []
E_4_list = []
E_5_list = []
E_1_list.append(math.log10(E(n, 9, q_real, p)))
E_2_list.append(math.log10(E(n, 10, q_real, p)))
E_3_list.append(math.log10(E(n, 11, q_real, p)))
E_4_list.append(math.log10(E(n, 12, q_real, p)))
E_5_list.append(math.log10(E(n, 13, q_real, p)))# 将t=0的能量添加
t_list = [i for i in range(161)]
resolution = 160000  # 分辨率要是160的倍数
dt = 2*math.pi/(2*math.sin(math.pi/2/(n+1)))*160/resolution
ddt = dt / 1000  # 第一次计算的时间分辨率
for i in range(500):  # 初始化计算
    q_old = q
    for i in range(n):
        q[i + 1] += ddt * p[i + 1]
    for i in range(n):  # 计算p和q的变化,假设得到的q是k+0.5的
        p[i + 1] += - ddt * (q_old[i + 1] - q_old[i + 2] + (q_old[i + 1] - q_old[i + 2]) ** 3
                            - q_old[i] + q_old[i + 1] - (q_old[i] - q_old[i + 1]) ** 3)
for j in range(resolution):
    print('完成', j/resolution*100, "%")
    q_old = q
    p = p_origin
    for i in range(n):  # 计算p和q的变化,得到的q是j+0.5的
        p[i + 1] += - dt * (q_old[i + 1] - q_old[i + 2] + (q_old[i + 1] - q_old[i + 2]) ** 3
                               - q_old[i] + q_old[i + 1] - (q_old[i] - q_old[i + 1]) ** 3)
    for i in range(n):
        q[i + 1] += dt * p[i + 1]
    for u in range(len(q)):  # 插值计算q
        q_real[u] = (q[u]+q_old[u])/2
    if (j+1)%(resolution/160) == 0:  # 在一定的位置再计算E
        E_1_list.append(math.log10(E(n, 9, q_real, p)))
        E_2_list.append(math.log10(E(n, 10, q_real, p)))
        E_3_list.append(math.log10(E(n, 11, q_real, p)))
        E_4_list.append(math.log10(E(n, 12, q_real, p)))
        E_5_list.append(math.log10(E(n, 13, q_real, p)))
plt.figure(1)
plt.plot(t_list,E_3_list,label='E11')
plt.plot(t_list,E_1_list,label='E9')
plt.plot(t_list,E_2_list,label='E10')
plt.plot(t_list,E_4_list,label='E12')
plt.plot(t_list,E_5_list,label='E13')
plt.xlabel('t')
plt.ylabel('lgE')
plt.title('energy change over time')
plt.legend()
plt.show()