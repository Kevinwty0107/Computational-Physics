import random
import math
import matplotlib.pyplot as plt
from copy import deepcopy
import matplotlib as mpl
mpl.rcParams['font.sans-serif'] = ['KaiTi']
mpl.rcParams['font.serif'] = ['KaiTi']
mpl.rcParams['axes.unicode_minus'] = False # 解决保存图像是负号'-'显示为方块的问题,或者转换负号为字符串



def listminusfunction(list1, list2):  # 向量相减
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

def V(x):
    return -1/math.sqrt(2+x**2)


def E2(t, N, omega):  # 第三题所用E
    return math.sqrt(2E14/3.5094448314E16) * (math.sin(omega*t/2/N))**2 * math.sin(omega*t)


def listnum_multiply(n, list):  # 向量数乘
    vec = []
    for i in range(0, len(list)):
        vec.append(n * list[i])
    return vec


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
boundary = 10  # 最大值
n = int(2*boundary/h)+1
a = []
b = []
c = []
for i in range(n):
    b.append(1/h**2 + V(-boundary+i * h)+0.48)
    a.append(-0.5/h**2)
    c.append(-0.5/h**2)
v = []
x = []
for i in range(n):
    v.append(random.random())
    x.append(-boundary+i * h)
u = deepcopy(v)
iteration = 1
while iteration == 1 or math.fabs(max_num(listminusfunction(u1, u))) > 1E-14:
    u1 = u
    v = chasesolvingfunction(n, a, b, c, u)
    lam = 1 / max_num(v)
    u = listnum_multiply(lam, v)
    iteration += 1
# 得到u,特征向量
sum_u = 0
for i in range(len(u)):
    sum_u += u[i]*u[i] * h
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
    return (-x/((2+x**2)**(1.5)))+E2(t_real, N, omega_light)


def d_function(t, i):  # 求d（t）的函数
    global h
    global boundary
    return -boundary + i * h

def cons_vector(lam, v):
    re = [lam * v[i] for i in range(len(v))]
    return re

def normalization(v, dx):  # 归一
    sum = 0
    for i in range(len(v)):
        sum += v[i] ** 2
    con = 1/math.sqrt(sum * dx)
    return cons_vector(con, v)


def d_or_a(y, function, N, dt):  # 本函数用来求解d或者a
    global omega_light
    global h
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
    dot = complex(0, 1) * dt / 2  # 后面用到的一个量i delta t/2
    a1 = []
    c1 = []
    a2 = []
    c2 = []
    for i in range(n):
        a1.append(-0.5 / h ** 2 * dot)
        c1.append(-0.5 / h ** 2 * dot)
        a2.append(0.5 / h ** 2 * dot)
        c2.append(0.5 / h ** 2 * dot)  # 初始化三对角矩阵，1是更新的时间下标,2是原时间下矩阵
    for t in range(int(2 * N * math.pi /omega_light/ dt) - 1):  # t不是真实时间
        sum_y=0
        b1 = []
        b2 = []
        for i in range(n):
            b1.append((1 / h ** 2 + V(-boundary + i * h) + (-boundary + i * h) * E2((t + 1) * dt, N,
                                                                                  omega_light)) * dot + 1)  # 现时刻的矩阵
            b2.append(-(1 / h ** 2 + V(-boundary + i * h) + (-boundary + i * h) * E2(t * dt, N,
                                                                                   omega_light)) * dot + 1)  # 原时刻的矩阵
        y_old = [b2[0] * y[t][0] + c2[0] * y[t][1]]
        for i in range(1, n - 1):
            y_old.append(a2[i] * y[t][i - 1] + b2[i] * y[t][i] + c2[i] * y[t][i + 1])
        y_old.append(a2[n - 1] * y[t][n - 2] + b2[n - 1] * y[t][n - 1])  # 计算演化时右侧的矩阵乘向量
        y[t+1] = chasesolvingfunction(n, a1, b1, c1, y_old)  # 追逐法求y(时间t时)
        for i in range(n):  # 计算吸收
            y[t+1][i] =( y[t+1][i] * gauss[i] )
        print('波函数求解已完成', t/(2*N*math.pi/omega_light/dt)*100, '%')
    return y



omega_300 = 45.5633525316/300
N = 48
lam = 300
omega_light = 45.5633525316/lam
dt = 0.05  # 时间步长
gauss = gauss_absorb(u)
y = psi_xt(u, dt, N, omega_light, gauss)# 得到的y是xt的函数，即各个时刻的波函数

# 下面部分是画偶极矩下的谐波功率谱
d_list = d_or_a(y, d_function, N, dt)
d_ome = 0.005  # ome的步长
ome_list = []
A_list = []

for ome in range(1, int(15*omega_300/d_ome)):  # 对不同频率积分
    ome_real = ome*d_ome
    ome_list.append(ome_real/omega_300)
    A = 0
    for t in range(int(2*N*math.pi/omega_light/dt)):
        A += complex(math.cos(ome_real * t * dt), -math.sin(ome_real*t * dt))*d_list[t]*dt
    A_list.append(math.log10((A.real**2+A.imag**2)*abs((1/math.sqrt(2*math.pi))*(-ome_real**2))**2))
    print('偶极矩已完成', ome/(15*omega_300/d_ome)*100, '%')
plt.plot(ome_list, A_list, label='偶极矩')
plt.ylabel('lg|A(omega)|^2')
plt.xlabel('omega/omega300')
plt.title('谐波辐射频谱')
# 下面部分是画电子加速度的谐波功率谱
a_list = d_or_a(y, a_function, N, dt)
print(a_list[2])
d_ome = 0.005  # ome的步长
ome_list = []
A_list = []
for ome in range(1, int(15*omega_300/d_ome)):  # 对不同频率积分
    ome_real = ome*d_ome
    ome_list.append(ome_real/omega_300)
    A = 0
    for t in range(int(2*N*math.pi/omega_light/dt)):
        A += (math.cos(ome_real * t * dt)-complex(0,1)*math.sin(ome_real*t * dt))*a_list[t]*dt
    A_list.append(math.log10((A.real**2+A.imag**2)/(2 * math.pi)))
    print('电子加速已完成', ome / (15 * omega_300 / d_ome) * 100, '%')
plt.plot(ome_list, A_list, label='电子加速')
plt.legend()
plt.yticks([])
plt.show()  # 画两个的图