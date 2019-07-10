import math
import matplotlib.pyplot as plt
from copy import deepcopy


def draw(beta, t_max, t_inter, t_run, t_start, t_stop):
    n = 128
    p = [0]
    q = [0]
    B = 0.5
    k = 11
    ome = 2*math.sin(math.pi*k/2/(n+1))
    for i in range(1, n+1):
        q.append(B * math.cos(math.pi * k * (i - n / 2) / (n + 1)) / math.cosh(math.sqrt(3 / 2) * B * ome * (i - n / 2)))
        p.append(B / math.cosh(math.sqrt(3 / 2) * B * ome * (i - n / 2))
               * (ome * (1 + 3 * (ome ** 2) * (B ** 2) / 16)
                  * math.sin(math.pi * k * (i - n / 2) / (n + 1)) + math.sqrt(3 / 2)
                  * B * math.cos(math.pi * k * (i - n / 2) / (n + 1))
                  * math.sin(math.pi * k / (n + 1)) * math.tanh(math.sqrt(3 / 2) * B * ome * (i - n / 2))))
    p.append(0)
    q.append(0)
    q_real = q  # 真实的q
    p_origin = p  # 记录初始p
    q_record = []
    q_record.append(q)
    t_list = [t_start+i*t_inter for i in range(int((t_stop-t_start)/t_inter)+1)]
    x_list = [i for i in range(130)]
    dt = 2*math.pi/(2*math.sin(math.pi/2/(n+1)))*t_run
    ddt = dt / 1000  # 第一次计算的时间分辨率
    for t in range(500):  # 初始化计算
        q_old = q
        for i in range(n):
            q[i + 1] += ddt * p[i + 1]
        for i in range(n):  # 计算p和q的变化,假设得到的q是k+0.5的
            p[i + 1] += - ddt * (q_old[i + 1] - q_old[i + 2] + beta * (q_old[i + 1] - q_old[i + 2]) ** 3
                                - q_old[i] + q_old[i + 1] - beta * (q_old[i] - q_old[i + 1]) ** 3)
    p = p_origin
    for j in range(int(t_max/t_run)):
        print('完成', j / int(t_max/t_run) * 100, "%")
        q_old = q
        for i in range(n):
            p[i + 1] += - dt * (q_old[i + 1] - q_old[i + 2] + beta * (q_old[i + 1] - q_old[i + 2]) ** 3
                                - q_old[i] + q_old[i + 1] - beta * (q_old[i] - q_old[i + 1]) ** 3)
        for i in range(n):
            q[i + 1] += dt * p[i + 1]
        for u in range(len(q)):
            q_real[u] = (q[u] + q_old[u]) / 2
        if j/int(t_max/t_run) > t_start/t_max and j/int(t_max/t_run) < t_stop/t_max and (j + 1) % (int(t_inter/t_run)) == 0:
            q_u = deepcopy(q_real)
            q_record.append(q_u)
    fig1 = plt.figure(1)
    surf1 = plt.contourf(x_list, t_list, q_record)
    fig1.colorbar(surf1)
    if beta == 0:
        plt.title('beta=0')
    else:
        plt.title('beta=1')
    plt.xlabel('x')
    plt.ylabel('t')
    plt.show()


draw(0, 1, 0.005, 0.0001, 0, 1)
draw(0, 10, 0.005, 0.0001, 9, 10)
draw(1, 1, 0.005, 0.0001, 0, 1)
draw(1, 10, 0.005, 0.0001, 9, 10)
