# 本程序用来计算束缚态的波函数
import math


def f(r):
    return (math.exp(-r / 3) * (2 * r / 3) * (4 - 2 * r / 3))


def differential1_low(h1, r):  # 定义了左端一阶导数
    return (f(r + h1) - f(r)) / h1


def differential1_high(h1, r):  # 定义了右端一阶导数
    return (f(r) - f(r - h1)) / h1


def differential2(h1, r):  # 定义了三点一阶导数
    return (f(r + h1) - f(r - h1)) / (2 * h1)


def differential3(h1, r):  # 定义了五点一阶导数
    return (f(r - 2 * h1) - 8 * f(r - h1) + 8 * f(r + h1) - f(r + 2 * h1)) / (12 * h1)


def differential1_two_low(h2, r):  # 定义了左端二阶导数
    return ((r + h2) ** 2 * differential2(h2 / 3, r + h2) - (r) ** 2 * differential1_low(h2 / 3, r)) / h2


def differential1_two_high(h2, r):  # 定义了右端二阶导数
    return ((r) ** 2 * differential1_high(h2 / 3, r) - (r - h2) ** 2 * differential2(h2 / 3, r - h2)) / h2


def differential3_two(h2, r):  # 定义了中间二阶导数
    return (8 * (r + h2) ** 2 * differential3(h2 / 3, r + h2) - 8 * (r - h2) ** 2 * differential3(h2 / 3, r - h2) +
            (r - 2 * h2) ** 2 * differential3(h2 / 3, r - 2 * h2) - (r + 2 * h2) ** 2 * differential3(h2 / 3,
                                                                                                      r + 2 * h2)) / (
           12 * h2)


def G1_low(r, h2):  # 定义了左端积分函数值
    G = f(r) * ((1 - r) * f(r) - 0.5 * differential1_two_low(h2, r))
    return G


def G1_high(r, h2):  # 定义了右端积分函数值
    G = f(r) * ((1 - r) * f(r) - 0.5 * differential1_two_high(h2, r))
    return G


def G3(r, h2):  # 定义了中间积分函数值
    G = f(r) * ((1 - r) * f(r) - 0.5 * differential3_two(h2, r))
    return G


def integral(h3):  # 复化牛顿-科斯特公式
    A = (7 * G1_high(60, h3 / 4) + 7 * G1_low(h3, h3 / 4))  # A是求和的左右端点值
    B=0
    C=0
    D=0
    E=0
    for i in range(int(60 / h3) ):
        B = B+ 32 * G3((i + 0.25) * h3, h3 / 4)
        C = C+12 * G3((i + 0.5) * h3, h3 / 4)
        D = D +32 * G3((i + 0.75) * h3, h3 / 4)


    for i in range(1, int(60 / h3) - 1):
        E = E+ 14 * G3((i) * h3, h3 / 4)

    intsum = h3 / 90 * (A + B + C + D + E)  # 对积分值进行加和
    return intsum


n = 3
l = 1
m = 1
h = 0.01  # 设置步长
C = (2 / n) ** 3 * math.factorial(n - l - 1) / (2 * n * math.factorial(n + l))
E = C * integral(h)  # 得到能量
print('l=1,m=1,n=3时的能量为：', E)

