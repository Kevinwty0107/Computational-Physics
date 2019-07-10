import random
import math
import matplotlib.pyplot as plt

I0 = 5 * 10 ** (13)
E0 = math.sqrt(I0 / (3.5094448314 * 10 ** (16)))
w = 45.5633525316 / 3200
Ip = 0.5
A0 = E0 / w
jie = [[0 for i in range(6)] for i in range(200)]  # 存放解的矩阵
p = [(i + 1) * 0.01 for i in range(200)]  # 可调参数p
Ek = [0 for i in range(200)]
for i in range(0, 200, 1):
    Ek[i] = p[i] ** 2 / 2  # 横坐标
tf = 4 * math.pi / w
Nx = 200
ti = 150
Ny = 60
MSPM = [0.0 for i in range(200)]
MDI = [0.0 for i in range(200)]


# SS=[[0 for i in range(10000)]for i in range(200)]


def exp(x):
    return complex(math.cos(x.real), math.sin(x.real)) * float(math.exp(-x.imag))  # 数学上exp(ix)


def csin(x):  # sinx,x为复数
    return (exp(x) - exp(-x)) / 2j


def ccos(x):  # 复数的余弦函数
    return (exp(x) + exp(-x)) / 2


def f(t, i):  # 就是本题中的x
    y = A0 * csin(w * t) * csin(w * t / (4 + 0j)) ** 2
    y = p[i] + y  # p可变
    y = 0.5 * y ** 2 + 0.5
    return y


def jiansuo(x, l, n):  # 用来检索,如果两个解相同，就返回1，不相同，返回1
    k = 0
    for i in range(0, n, 1):
        a = x.real - jie[l][i].real  # 分别对比实部和虚部
        b = x.imag - jie[l][i].imag
        if a < 10 ** (-3) and b < 10 ** (-3):
            k = 1
            break
    return k


def fjifen(t, p):  # 用来存要积分 的方程，也就是第一问为0的方程
    return (p + A(t)) ** 2 / 2 + Ip  # 为时间t，和参数pz的二元函数


def S(t, p):  # 题目中的S函数，是第一问函数，也即f的积分，是t和pz的二元函数
    dt = t / 1000  # 分成1000份
    a = 0
    for i in range(1, 1001, 1):
        a += 4 * fjifen((i - 0.5) * dt, p) + 2 * fjifen(i * dt, p)  # 采用复化抛物求积
    a = (a + fjifen(0, p) - fjifen(t, p)) * dt / 6
    return a  # 返回积分值


def A(t):  # 为了方便，把电磁四矢量存为一个函数
    return A0 * csin(w * t) * csin(w * t / (4 + 0j)) ** 2  # 返回题目中A(t)


def Sdd(t, p):  # S函数的二阶导函数,二阶导的解析表示就是下面方程右边
    a = (p + A(t)) * w * A0 * (ccos(w * t) * csin(w * t / 4) ** 2 + 1 / 4 * csin(w * t) * csin(w * t / 2))
    return a


def m(t, p):  # 返回鞍点近似下每一个M求和中的一项
    a = exp(S(t, p)) / Sdd(t, p)
    return a


def mdi(t, p):  # 直接求积分算MDI的被积函数
    q11 = p + A(t)
    a = (-1) * q11 * w * A0 * (ccos(w * t) * csin(w * t / 4) ** 2 + 1 / 4 * csin(w * t) * csin(w * t / 2))
    b = math.pi * (q11 ** 2 + 2 * Ip) ** 3
    d = a * exp(S(t, p)) / b
    return d


def mjifen(t, p):  # 复化抛物求积求MDI
    dt = t / 5000  # 分成10000份
    a = 0
    for i in range(1, 5001, 1):
        a += 4 * mdi((i - 0.5) * dt, p) + 2 * mdi(i * dt, p)
    a = (a + mdi(0, p) - mdi(t, p)) * dt / 6
    return a * 2 ** (7 / 2)  # 返回MDI积分的值


nj = 0  # 解的个数
for i in range(0, Nx, 1):
    for j in range(0, Ny, 1):
        a = complex(i * tf / Nx, j * ti / Ny)  # 初始两个端点随机数产生
        b = complex((i + 1) * tf / Nx, (j + 1) * ti / Ny)
        x0 = complex(random.uniform(i * tf / Nx, (i + 1) * tf / Nx), random.uniform(j * ti / Ny, (j + 1) * ti / Ny))
        sd = 0  # 存放while循环步数
        k = 0  # 存放判定数k，在while循环里，一旦出现判定结果为break的条件，k=1
        while sd == 0 or math.sqrt(f(x0, 0).real ** 2 + f(x0, 0).imag ** 2) > 10 ** (-12):
            x0 = x0 - (b - a) * f(x0, 0) / (f(b, 0) - f(a, 0))
            sd += 1
            try:
                ans = f(x0, 0)
            except OverflowError:  # 防止溢出，出现溢出，break，k=1
                k = 1
                break
            try:
                ans = math.sqrt(f(x0, 0).real ** 2 + f(x0, 0).imag ** 2)
            except OverflowError:
                k = 1
                break
            if sd > 100 or abs(x0.imag) > 100 * math.pi:  # 循环次数控制在100以内，一般能收敛，100次就足够用
                k = 1
                break
        if x0.imag > 0 and 0 < x0.real < 4 * math.pi / w and k == 0 and jiansuo(x0, 0, nj) == 0:  # 满足t的各个条件才能把x赋在解向量里
            print(sd)  # 输出步数，一般可以证明20步以内即可
            jie[0][nj] = x0  # 存入解向量之中
            nj += 1
            print(x0)

for l in range(1, 200, 1):
    nj = 0
    while nj < 6:
        for i in range(6):
            if jie[l - 1][i].real < 1:
                a = jie[l - 1][i]  # 直接取上一步的值
            else:
                a = jie[l - 1][i] - 1 - 1j
            b = jie[l - 1][i] + 1 + 1j  # 另一个端点取稍大的值
            x0 = complex(random.uniform(a.real, b.real), random.uniform(a.imag, b.imag))
            # x0=jie[l-1][i]
            sd = 0
            k = 0
            while sd == 0 or math.sqrt(f(x0, l).real ** 2 + f(x0, l).imag ** 2) > 10 ** (-12):
                x0 = x0 - (b - a) * f(x0, l) / (f(b, l) - f(a, l))
                sd += 1
                try:
                    ans = f(x0, l)
                except OverflowError:  # 防止溢出，出现溢出，break，k=1
                    k = 1
                    break
                try:
                    ans = math.sqrt(f(x0, l).real ** 2 + f(x0, l).imag ** 2)
                except OverflowError:
                    k = 1
                    break
                if sd > 20 or abs(x0.imag) > 100 * math.pi:  # 循环次数控制在20以内，一般能收敛，20次就足够用
                    k = 1
                    break
            if x0.imag > 0 and 0 < x0.real < 4 * math.pi / w and k == 0 and jiansuo(x0, l,
                                                                                    nj) == 0:  # 满足t的各个条件才能把x赋在解向量里
                print(sd)  # 输出步数，一般可以证明20步以内即可
                jie[l][nj] = x0  # 存入解向量之中
                nj += 1
                print(x0)

for i in range(0, 200, 1):
    a = 0.0
    print(i)
    for j in range(0, 6, 1):
        a += m(jie[i][j], p[i])  # 每个pz，6个解
    MSPM[i] = (a.real ** 2 + a.imag ** 2) / 2  # 得到的MPSM模方, 前面系数根号2在这里统一除掉
    b = mjifen(4 * math.pi / w, p[i])  # 得到MDI模方
    MDI[i] = b.real ** 2 + b.imag ** 2

plt.plot(Ek, MSPM, label="MSPM")
plt.plot(Ek, MDI, label="MDI")
plt.xlabel('Ek')
plt.ylabel('|M|2')
plt.legend()
plt.show()





























