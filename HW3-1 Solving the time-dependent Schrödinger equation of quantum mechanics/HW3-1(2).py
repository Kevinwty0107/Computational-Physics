import random
import math
import matplotlib.pyplot as plt
from copy import deepcopy
import matplotlib as mpl
mpl.rcParams['font.sans-serif'] = ['KaiTi']
mpl.rcParams['font.serif'] = ['KaiTi']
mpl.rcParams['axes.unicode_minus'] = False # 解决保存图像是负号'-'显示为方块的问题,或者转换负号为字符串


def listnum_multiply(n, list):  # 向量数乘
    vec = []
    for i in range(0, len(list)):
        vec.append(n * list[i])
    return vec


def listminusfunction(list1, list2):  # 向量相减
    list_return = [0] * len(list2)
    if len(list1) != len(list2):
        print("length does not match")
        exit()
    else:
        for n in range(len(list1)):
            list_return[n] = list1[n] - list2[n]
    return list_return


def V(x):
    return -1/math.sqrt(2+x**2)


def E(t, N):
    return math.sqrt(1E16/3.5094448314E16) * (math.sin(t/2/N))**2 * math.sin(t)



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


h = 0.1  # 步长
boundary = 200  # 求解的x_max
n = int(2*boundary/h)+1  # 区间数
a = []
b = []
c = []
for i in range(n):  # 初始化追逐法用到的三列向量
    b.append(1/h**2 + V(-boundary+i*h)+0.48)
    a.append(-0.5/h**2)
    c.append(-0.5/h**2)
v = []
x = []
for i in range(n):  # 初始化反幂法用到的解向量v以及画图时要用到的坐标向量x
    v.append(random.random())
    x.append(-boundary+i*h)
u = deepcopy(v)
iteration = 1
while iteration == 1 or math.fabs(max_num(listminusfunction(u1, u))) > 1E-14:  # 运用反幂法求解，其中LU分解解方程的部分已经化为追逐法
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

# 下面为第二问
a1 = []
c1 = []
a2 = []
c2 = []
dt = 0.05  # 时间步长
N = 18
t_list = []  # 初始化时间列表
p_list = []  # 初始化
y = deepcopy(u)  # y的初始态是u
dot = complex(0, 1)*dt/2  # 后面用到的一个量i delta t/2
y_origin = deepcopy(u)  # 基态的波函数

for i in range(n):
    a1.append(-0.5 / h ** 2 * dot)
    c1.append(-0.5 / h ** 2 * dot)
    a2.append(0.5 / h ** 2 * dot)
    c2.append(0.5 / h ** 2 * dot)  # 初始化三对角矩阵，1是更新的时间下标,2是原时间下矩阵，
for t in range(int(2*N*math.pi/dt)):  # t不是真实时间
    b1 = []
    b2 = []
    t_list.append(t * dt)
    for i in range(n):
        b1.append((1 / h ** 2 + V(-boundary + i * h) + (-boundary + i * h)*E((t+1)*dt, N))*dot + 1)  # 现时刻的矩阵
        b2.append(-(1 / h ** 2 + V(-boundary + i * h) + (-boundary + i * h)*E(t*dt, N))*dot + 1)  # 原时刻的矩阵
    y_old = [b2[0] * y[0] + c2[0] * y[1]]
    for i in range(1, n-1):
        y_old.append(a2[i] * y[i-1] + b2[i] * y[i] + c2[i] * y[i+1])
    y_old.append(a2[n - 1] * y[n - 2] + b2[n - 1] * y[n - 1])  # 计算演化时右侧的矩阵乘向量
    y = chasesolvingfunction(n, a1, b1, c1, y_old)  # 追逐法求y
    p = 0  # p是最后的求和
    for i in range(n):
        p += y[i] * y_origin[i] * h
    p_list.append(abs(p)**2)
    print('布居数已完成：', t/(2*N*math.pi/dt)*100, '%')
plt.plot(t_list, p_list)
plt.ylabel('P')
plt.xlabel('t')
plt.title('布居数随时间变化图')
plt.show()  # 画布居数的图

psi_f = []
for i in range(n):
    psi_f.append(y[i]-p*y_origin[i])
P_k = []
dk = 0.01
k_list = []
for k in range(int(5/dk)+1):
    sum_k = 0
    k_list.append(k*dk-2.5)
    for i in range(n):
        sum_k += (math.cos((k*dk-2.5)*(-boundary + i * h))-complex(0,1)*math.sin((k*dk-2.5)*(-boundary + i * h)))*psi_f[i]*h
    P_k.append(abs(sum_k)**2)
    print('动量谱已完成：', k / (5/dk) * 100, '%')
plt.plot(k_list, P_k)
plt.ylabel('P(K)')
plt.xlabel('k')
plt.title('电离电子动量谱')
plt.show()  # 画动量谱的图

for i in range(0,len(P_k)):
    P_k[i] = math.log10(P_k[i])

plt.plot(k_list, P_k)
plt.ylabel('logP(K)')
plt.xlabel('k')
plt.title('电离电子动量谱（对数坐标）')
plt.show()  # 画动量谱的图

