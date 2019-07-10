import math
import matplotlib.pyplot as plt


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
    q.append(20*math.sin(t*math.pi/(n+1))/math.sqrt(n/2) *n/(n+1))
q_real = q  # 真实的q
p_origin = p  # 记录初始p
E_1_list = []
E_2_list = []
E_3_list = []
E_4_list = []
E_1_list.append(E(n, 1, q_real, p))
E_2_list.append(E(n, 2, q_real, p))
E_3_list.append(E(n, 3, q_real, p))
E_4_list.append(E(n, 4, q_real, p))  # 将t=0的能量添加
t_list = [i for i in range(2001)]
resolution = 4000000  # 分辨率要是160的倍数
dt = 2*math.pi/(2*math.sin(math.pi/2/(n+1)))*2000/resolution
ddt = dt / 1000  # 第一次计算的时间分辨率
for i in range(500):  # 初始化计算
    q_old = q
    for i in range(n):
        q[i + 1] += ddt * p[i + 1]
    for i in range(n):  # 计算p和q的变化,假设得到的q是k+0.5的
        p[i + 1] += - ddt * (q_old[i + 1] - q_old[i + 2] + 0.6 * (q_old[i + 1] - q_old[i + 2]) ** 3
                            - q_old[i] + q_old[i + 1] - 0.6 * (q_old[i] - q_old[i + 1]) ** 3)
for j in range(resolution):
    print('完成', j/resolution*100, "%")
    q_old = q
    p = p_origin
    for i in range(n):  # 计算p和q的变化,得到的q是j+0.5的
        p[i + 1] += - dt * (q_old[i + 1] - q_old[i + 2] + 0.6 * (q_old[i + 1] - q_old[i + 2]) ** 3
                               - q_old[i] + q_old[i + 1] - 0.6 * (q_old[i] - q_old[i + 1]) ** 3)
    for i in range(n):
        q[i + 1] += dt * p[i + 1]
    for u in range(len(q)):  # 插值计算q
        q_real[u] = (q[u]+q_old[u])/2
    if (j+1)%(resolution/2000) == 0:  # 在一定的位置再计算E
        E_1_list.append(E(n, 1, q_real, p))
        E_2_list.append(E(n, 2, q_real, p))
        E_3_list.append(E(n, 3, q_real, p))
        E_4_list.append(E(n, 4, q_real, p))
plt.figure(1)
plt.plot(t_list,E_1_list,label='E1')
plt.plot(t_list,E_2_list,label='E2')
plt.plot(t_list,E_3_list,label='E3')
plt.plot(t_list,E_4_list,label='E4')
plt.xlabel('t')
plt.ylabel('E')
plt.title('energy change over time')
plt.legend()
plt.show()