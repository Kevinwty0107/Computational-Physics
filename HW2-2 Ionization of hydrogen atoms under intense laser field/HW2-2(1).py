import random
import math

I0 = 5 * 10 ** (13)
E0 = math.sqrt(I0 / (3.5094448314 * 10 ** (16)))
w = 45.5633525316 / 3200
Ip = 0.5
A0 = E0 / w
jie = [0 for i in range(6)]  # 存放解的矩阵
p = 1  # 可调参数p
tf = 4 * math.pi / w
Nx = 200
ti = 150
Ny = 60


def exp(x):
    return complex(math.cos(x.real), math.sin(x.real)) * float(math.exp(-x.imag))  # 数学上exp(ix)


def csin(x):  # sinx,x为复数
    return (exp(x) - exp(-x)) / 2j


def f(t):  # 就是本题中的x
    y = A0 * csin(w * t) * csin(w * t / complex(4,0)) ** 2
    y = p + y  # p可变
    y = 0.5 * y ** 2 + 0.5
    return y


def jiansuo(x, n, jie):  # 用来检索,如果两个解相同，就返回1，不相同，返回0
    k = 0
    for i in range(0, n, 1):
        a = x.real - jie[i].real  # 分别对比实部和虚部
        b = x.imag - jie[i].imag
        if a < 10 ** (-3) and b < 10 ** (-3):
            k = 1
            break
    return k


nj = 0  # 解的个数
for i in range(0, Nx, 1):
    for j in range(0, Ny, 1):  # 在解空间切割小空间，产生端点
        a = complex(i * tf / Nx, j * ti / Ny)  # 初始两个端点随机数产生
        b = complex((i + 1) * tf / Nx, (j + 1) * ti / Ny)
        x0 = complex(random.uniform(i * tf / Nx, (i + 1) * tf / Nx), random.uniform(j * ti / Ny, (j + 1) * ti / Ny))
        sd = 0  # 存放while循环步数
        k = 0  # 存放判定数k，在while循环里，一旦出现判定结果为break的条件，k=1
        while sd == 0 or math.sqrt(f(x0).real ** 2 + f(x0).imag ** 2) > 10 ** (-12):
            x0 = x0 - (b - a) * f(x0) / (f(b) - f(a))
            sd += 1
            try:
                ans = f(x0)
            except OverflowError:  # 防止溢出，出现溢出，break，k=1
                k = 1
                break
            try:
                ans = math.sqrt(f(x0).real ** 2 + f(x0).imag ** 2)
            except OverflowError:
                k = 1
                break
            if sd > 100 or abs(x0.imag) > 100 * math.pi:  # 循环次数控制在100以内，一般能收敛，100次就足够用
                k = 1
                break
        if x0.imag > 0 and x0.real>0 and  x0.real < 4 * math.pi / w and k == 0 and jiansuo(x0, nj, jie) == 0:  # 满足t的各个条件才能把x赋在解向量里
            #print(sd)  # 输出步数，一般可以证明100步以内即可
            jie[nj] = x0  # 存入解向量之中
            nj += 1
            print(x0)

print(jie)































