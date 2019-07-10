import random
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['font.sans-serif'] = ['KaiTi']
mpl.rcParams['font.serif'] = ['KaiTi']
mpl.rcParams['axes.unicode_minus'] = False # 解决保存图像是负号'-'显示为方块的问题,或者转换负号为字符串

from copy import deepcopy


def listminusfunction(list1, list2):  # 向量相减
    list_return = [0] * len(list2)
    if len(list1) != len(list2):
        print("length does not match")
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


h = 0.1  # 步长
boundary = 2000  # 求解的x_max
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
# 得到u,特征向量，lam是特征值
sum_u = 0
for i in range(len(u)):
    sum_u += u[i]*u[i]
for i in range(len(u)):  # 对得到的基态波函数进行归一化
    u[i] = u[i]*u[i]/sum_u/h
print('得到的氢原子基态能量为：', lam-0.48)
plt.plot(x, u)
plt.xlabel('x')
plt.ylabel('概率')
plt.title('氢原子基态时在各处出现的概率分布图')
#plt.legend()
plt.show()
