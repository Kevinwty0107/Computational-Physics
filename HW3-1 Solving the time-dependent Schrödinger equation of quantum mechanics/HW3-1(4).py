import random
import math
import matplotlib.pyplot as plt
from copy import deepcopy
import numpy as np
import matplotlib as mpl
mpl.rcParams['font.sans-serif'] = ['KaiTi']
mpl.rcParams['font.serif'] = ['KaiTi']
mpl.rcParams['axes.unicode_minus'] = False # 解决保存图像是负号'-'显示为方块的问题,或者转换负号为字符串


def V(x):
    return -1/math.sqrt(2+x**2)


def E2(t, N, omega):  # 第三题所用E
    return math.sqrt(1E14/3.5094448314E16) * (math.sin(omega*t/2/N))**2 * math.sin(omega*t)


def listnum_multiply(n, list):  # 向量数乘
    vec = []
    for i in range(0, len(list)):
        vec.append(n * list[i])
    return vec


def listminus_function(list1, list2):  # 向量相减
    list_return = [0] * len(list2)
    if len(list1) != len(list2):
        print("List the length does not match")
        exit()
    else:
        for n in range(len(list1)):
            list_return[n] = list1[n] - list2[n]
    return list_return


def max_num(v):  # 向量最大值
    n = len(v)
    res = 0
    num = 0
    for i in range(0, n):
        if math.fabs(v[i]) > res:
            res = math.fabs(v[i])
            num = i
    return v[num]


def chasesolvingfunction(n, a, b, c, f):  # 定义追逐法解方程的函数
    a1 = deepcopy(a)
    b1 = deepcopy(b)
    c1 = deepcopy(c)
    f1 = deepcopy(f)
    for i in range(1, n):
        m = a1[i] / b1[i-1]
        b1[i] = b1[i] - m * c1[i-1]
        f1[i] = f1[i] - m * f1[i-1]
    x = [0] * n
    x[n-1] = f1[n-1]/b1[n-1]
    for i in range(n-2, -1, -1):
        x[i] = (f1[i] - c1[i] * x[i + 1])/b1[i]
    return x


h = 0.1
boundary = 600  # 最大值
n = int(2*boundary/h)+1
a = []
b = []
c = []
v = []
x = []
for i in range(n):
    b.append(1/h**2 + V(-boundary+i*h)+0.48)
    a.append(-0.5/h**2)
    c.append(-0.5/h**2)
for i in range(n):
    v.append(random.random())
    x.append(-boundary+i*h)
u = deepcopy(v)
iteration = 1
while iteration == 1 or math.fabs(max_num(listminus_function(u1, u))) > 1E-14:
    u1 = u
    v = chasesolvingfunction(n, a, b, c, u)
    lam = 1 / max_num(v)
    u = listnum_multiply(lam, v)
    iteration += 1
# 得到u,特征向量
sum_u = 0
for i in range(len(u)):
    sum_u += u[i]*u[i]*h
for i in range(len(u)):
    u[i] = u[i]/math.sqrt(sum_u)  # 得到基态波函数


def a_function(t, i):  # 求a（t）的函数
    global h
    global boundary
    global dt
    global N
    global omega_light
    x = (-boundary + i * h)
    t_real = t*dt
    return -x/(2+x*x)**(1.5)+E2(t_real, N, omega_light)


def d_function(t, i):  # 求d（t）的函数
    global h
    global boundary
    return -boundary + i * h


def d_or_a(y, function, N, dt):  # 本函数用来求解d或者a
    global omega_light
    d_list = []  # 初始化，d是t的函数
    d = 0  # d是最后的求和
    for t in range(int(2 * N * math.pi /omega_light / dt)):
        for i in range(n):
            d += y[t][i].conjugate() * function(t, i) * y[t][i] * h
        d_list.append(d)  # 得到t关于时间的函数
    return d_list


def gauss_absorb(u):  # 高斯吸收序列
    global h
    global boundary
    gauss = []
    for i in range(len(u)):
        x = -boundary + i * h
        if abs(x) > 0.75 * boundary:
            gauss.append(math.exp(-(abs(x)-0.75*boundary)**2/0.04))
        else:
            gauss.append(1)
    return gauss


def psi_xt(u, dt, N, omega_light, gauss):  # 本函数用来求解波函数,输出的长度是int(2 * N * math.pi /omega_light/ dt)
    global boundary
    y = [[] for t in range(int(2 * N * math.pi /omega_light/ dt))]
    y[0] = deepcopy(u)  # 引入第1问中得到的波函数
    const = complex(0, 1) * dt / 2  # 后面用到的一个量i delta t/2
    a1 = []
    c1 = []
    a2 = []
    c2 = []
    for i in range(n):
        a1.append(-0.5 / h ** 2 * const)
        c1.append(-0.5 / h ** 2 * const)
        a2.append(0.5 / h ** 2 * const)
        c2.append(0.5 / h ** 2 * const)  # 初始化三对角矩阵，1是更新的时间下标,2是原时间下矩阵
    for t in range(int(2 * N * math.pi /omega_light/ dt) - 1):  # t不是真实时间
        b1 = []
        b2 = []
        for i in range(n):
            b1.append((1 / h ** 2 + V(-boundary + i * h) + (-boundary + i * h) * E2((t + 1) * dt, N,
                                                                                  omega_light)) * const + 1)  # 现时刻的矩阵
            b2.append(-(1 / h ** 2 + V(-boundary + i * h) + (-boundary + i * h) * E2(t * dt, N,
                                                                                   omega_light)) * const + 1)  # 原时刻的矩阵
        y_old = [b2[0] * y[t][0] + c2[0] * y[t][1]]
        for i in range(1, n - 1):
            y_old.append(a2[i] * y[t][i - 1] + b2[i] * y[t][i] + c2[i] * y[t][i + 1])
        y_old.append(a2[n - 1] * y[t][n - 2] + b2[n - 1] * y[t][n - 1])  # 计算演化时右侧的矩阵乘向量
        y[t+1] = chasesolvingfunction(n, a1, b1, c1, y_old)  # 追逐法求y(时间t时)
        for i in range(n):  # 计算吸收
            y[t+1][i] = y[t+1][i] * gauss[i]
        print('波函数求解已完成', t/(2*N*math.pi/omega_light/dt)*100, '%')
    return y


omega_300 = 45.5633525316/300
N = 4
lam = 400
omega_light = 45.5633525316/lam
dt = 0.1  # 时间步长
gauss = gauss_absorb(u)
y = psi_xt(u, dt, N, omega_light, gauss)  # 得到的y是xt的函数，即各个时刻的波函数
a_list = d_or_a(y, a_function, N, dt)
t_0_list = []
d_ome = 0.01  # ome的步长
dt_0 = 2
A_list = [[] for t_0 in range(int(2 * N * math.pi /omega_light/dt_0))]
ome_list = []
for ome in range(1, int(20*omega_300/d_ome)):
    ome_real = ome * d_ome
    ome_list.append(ome_real / omega_300)
for t_0 in range(int(2 * N * math.pi / omega_light / dt_0)):  # 对不同时刻循环
    t_0_list.append(t_0*dt_0)
    for ome in range(1, int(20*omega_300/d_ome)):  # 对不同频率循环
        ome_real = ome*d_ome
        A = 0
        for t in range(int(2*N*math.pi/omega_light/dt)):
            A += complex(math.cos(ome_real * t * dt),-math.sin(ome_real*t * dt))*a_list[t]*dt*math.exp(-(t*dt-t_0*dt_0)**2/2/15**2)
        A_list[t_0].append(math.log10(abs(A)**2))
    print('功率谱已完成', t_0 / (2 * N * math.pi /omega_light/dt_0) * 100, '%')
fig2 = plt.figure(2)
surf2 = plt.contourf(ome_list, t_0_list, A_list)
fig2.colorbar(surf2)
plt.title('具有时间分辨率的谐波功率谱分布')
plt.xlabel('omega/omega300')
plt.ylabel('t')
plt.show()
