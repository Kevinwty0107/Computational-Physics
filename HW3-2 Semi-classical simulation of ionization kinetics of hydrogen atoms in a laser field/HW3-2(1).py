import random
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['font.sans-serif'] = ['KaiTi']
mpl.rcParams['font.serif'] = ['KaiTi']
mpl.rcParams['axes.unicode_minus'] = False # 解决保存图像是负号'-'显示为方块的问题,或者转换负号为字符串



def E(t):  # 线偏振电场函数
    omega = 0.057
    E = -1.325 * omega * (
            0.5 * math.cos(omega * t) + 0.5 * math.cos(omega * t) * math.cos(omega * t / 4) - 1 / 8 * math.sin(
        omega * t) * math.sin(omega * t / 4))
    return E


def E_ellipse(t):  # 椭圆偏振电场函数
    omega = 0.057
    E_x = -1/math.sqrt(1.25) * 1.325 * omega * (
            0.5 * math.cos(omega * t) + 0.5 * math.cos(omega * t) * math.cos(omega * t / 4) - 1 / 8 * math.sin(
        omega * t) * math.sin(omega * t / 4))
    E_y = -0.5 / math.sqrt(1.25) * 1.325 * omega * (
            -0.5 * math.sin(omega * t) - 0.5 * math.sin(omega * t) * math.cos(omega * t / 4) - 1 / 8 * math.cos(
        omega * t) * math.sin(omega * t / 4))
    return E_x,E_y


def species(E):  # 对于固定的E选取一些v的样本
    max = 1 / 0.076 ** 2 * math.exp(-(2 / 3) / 0.076)  # max是全局的
    v_list = []  # 初始化列表
    for iter_a in range(1000):
        Chisi = random.random()  # 随机生成的参考向量
        v = random.uniform(-1, 1)  # 随机选取v，在-1到1之间
        W = 1 / E ** 2 * math.exp(-(2 / 3 + v ** 2) / E)  # 隧穿几率
        if W / max > Chisi:  # 判定是否选取
            v_list.append(v)
    return v_list  # 返回一系列v的样本



t1_list = []
t2_list = []
v1_list = []
v2_list = []
omega = 0.057
resolution = 10000
t_0 = -4*math.pi/omega
for n in range(resolution):  # 对于时间t做循环，计算末状态
    t = t_0 + 8*math.pi/omega * n / resolution
    print('已完成计算：', n/resolution*100, '%')
    E_1 = abs(E(t))
    E_2x, E_2y = E_ellipse(t)
    E_2 = math.sqrt(E_2x**2+E_2y**2)
    if E_1 > 0.03:
        v_ver1 = species(E_1)
    else:
        v_ver1 = []
    if E_2 > 0.03:
        v_ver2 = species(E_2)
    else:
        v_ver2 = []
    for i in range(len(v_ver1)):
        t1_list.append(t)
    v1_list = v1_list + v_ver1
    for i in range(len(v_ver2)):
        t2_list.append(t)
    v2_list = v2_list + v_ver2
plt.figure(1)
plt.scatter(t1_list, v1_list, s=1)
plt.xlabel('t')
plt.ylabel('v')
plt.title('line')
plt.show()
plt.figure(2)
plt.scatter(t2_list, v2_list, s=1)
plt.xlabel('t')
plt.ylabel('v')
plt.title('ellipse')
plt.show()

max = 1/0.076**2*math.exp(-(2/3)/0.076)  # 最大概率
E_list = []
v_list = []  # 初始化列表
for iter_a in range(500000):
    Chisi = random.random()  # 随机生成的参考向量
    E = random.uniform(0, 0.076)
    v = random.uniform(0, 1)
    W = 1/E**2*math.exp(-(2/3+v**2)/E)  # 隧穿几率
    if W/max > Chisi:  # 判定是否选取
        E_list.append(E)
        v_list.append(v)
plt.figure(3)
plt.scatter(E_list, v_list, s=1)
plt.xlabel('E')
plt.ylabel('v')
plt.title('电子随电场和垂直电场速度分布图')
plt.show()