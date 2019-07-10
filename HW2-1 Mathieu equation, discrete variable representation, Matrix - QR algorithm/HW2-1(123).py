#!/usr/bin/python
# -*- coding: UTF-8 -*-


import math
import cmath
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from pylab import *  # 导入pylab，构建类似于MATLAB的环境,用于画图

temp = input("请输入M值： ")
M = int(temp)
n = 2 * M + 1
phi = [0 for i in range(2 * M + 2)]  # 就是基函数里的phik格点，第0个设为0，然后从第一个按讲义上写
for i in range(1, 2 * M + 2, 1):
    phi[i] = (2 * math.pi * i) / (2 * M + 1)  # 赋值
ql = 0  # ql、qr分别是参数q的左右选值
qr = 21
pxl = [[0 for i in range(qr - ql)] for i in range(n)]  # pxl存储q下前几个偶函数的本征值
sxl = [[0 for i in range(qr - ql)] for i in range(n)]  # pxl存储q下本征值
bzxl = [[0 for i in range(qr - ql)] for i in range(n)]  # baxl用来存参数q下从小到大第几个本征值的本征矢量在本征值矩阵中（下面Q）的列数，通过访问Q就知道本征向量的具体数值


def sign1(a):  # sign1就是返回正负的函数
    if a > 0: return 1
    if a < 0: return -1


def fs2(a, m):  # 定义2范数函数，返回模长
    c = 0
    for i in range(0, m, 1):
        c += a[i] ** 2
    c = math.sqrt(c)
    return c


def jdz(a):  # 返回函数的绝对值
    if a >= 0:
        return a
    else:
        return -a


for q in range(ql, qr, 1):  # 从ql到qr（qr取不到），循环求解矩阵
    A = [[0 for i in range(2 * M + 2)] for i in range(2 * M + 2)]
    for j in range(1, 2 * M + 2, 1):  # 赋值
        for k in range(1, 2 * M + 2, 1):
            if j == k:
                A[j][j] = M * (M + 1) / 3 + 2 * q * math.cos(2 * phi[j])
            else:
                A[j][k] = ((-1) ** (j - k)) * math.cos(math.pi * (j - k) / (2 * M + 1)) / (
                2 * math.sin(math.pi * (j - k) / (2 * M + 1)) ** 2)
    A = np.asarray(A)

    Uk = [[0 for i in range(2 * M + 1)] for i in range(2 * M + 1)]  # 用来存储每一步的H反射向量，累乘即可得到本征向量
    for i in range(0, 2 * M + 1, 1):
        Uk[i][i] = 1.0
    Uk = np.asarray(Uk)

    for k in range(1, 2 * M, 1):
        x = A[k + 1:2 * M + 2, k]
        e1 = [0 for i in range(2 * M - k + 1)]
        e1[0] = 1.0
        e1 = np.asarray(e1)
        v = sign1(x[0]) * fs2(x, 2 * M - k + 1) * e1 + x  # 计算向量v
        v = v / fs2(v, 2 * M - k + 1)
        vvT = [[0 for i in range(2 * M - k + 1)] for i in range(2 * M - k + 1)]  # vvT矩阵
        Uk0 = [[0 for i in range(2 * M + 2)] for i in range(2 * M + 2)]  # 用来存储每一步的H反射矩阵
        for i in range(1, 2 * M + 2, 1):
            Uk0[i][i] = 1.0
        Uk0 = np.asarray(Uk0)
        for i in range(0, 2 * M - k + 1, 1):
            for j in range(0, 2 * M - k + 1, 1):
                vvT[i][j] = v[i] * v[j]  # vvT矩阵
        vvT = np.asarray(vvT)
        Uk0[k + 1:2 * M + 2, k + 1:2 * M + 2] = Uk0[k + 1:2 * M + 2, k + 1:2 * M + 2] - 2 * vvT  # UK0部分赋予反射部分
        Uk0 = Uk0[1:2 * M + 2, 1:2 * M + 2]
        Uk = np.matmul(Uk, Uk0)  # 每次累乘得到最终的H反射矩阵
        A[k + 1:2 * M + 2, k:2 * M + 2] = A[k + 1:2 * M + 2, k:2 * M + 2] - 2 * np.matmul(vvT, A[k + 1:2 * M + 2,
                                                                                               k:2 * M + 2])  # H反射矩阵对A的分块矩阵左乘
        A[1:2 * M + 2, k + 1:2 * M + 2] = A[1:2 * M + 2, k + 1:2 * M + 2] - 2 * np.matmul(
            A[1:2 * M + 2, k + 1:2 * M + 2], vvT)  # H反射矩阵对A的分块矩阵右乘
    A = np.asarray(A[1:2 * M + 2, 1:2 * M + 2])
    for i in range(0, 2 * M + 1, 1):
        for j in range(0, 2 * M + 1, 1):
            if jdz(A[i][j]) < 10 ** (-10):
                A[i][j] = 0  # 原则上此时A应该是三对角矩阵，但因为计算原因，会保留非常小的数，此处直接把其赋为0

    GT = [[0 for i in range(n)] for i in range(n)]  # Givens矩阵，下面累乘就能获得最终Givens矩阵
    for i in range(0, n, 1):
        GT[i][i] = 1.0
    GT = np.asarray(GT)
    k = 2 * M
    while k > 0:
        if jdz(A[k][k - 1]) < 10 ** (-12):
            k = k - 1
        s = A[k][k]
        for j in range(0, k + 1, 1):
            A[j][j] -= s  # 原点位移
        for i in range(1, k + 1, 1):
            a = math.sqrt(A[i - 1][i - 1] ** 2 + A[i][i - 1] ** 2)
            c = A[i - 1][i - 1] / a
            s1 = A[i][i - 1] / a
            G = [[0 for i in range(k + 1)] for i in range(k + 1)]
            for j in range(0, k + 1, 1):
                G[j][j] = 1.0
            G[i][i] = G[i - 1][i - 1] = A[i - 1][i - 1] / a
            G[i][i - 1] = -s1
            G[i - 1][i] = s1  # 单个Givens矩阵赋值
            G = np.asarray(G)
            A[0:k + 1, 0:k + 1] = np.matmul(G, A[0:k + 1, 0:k + 1])  # Givens矩阵左乘A的分块矩阵
            Gni = G
            Gni[i][i - 1] = s1
            Gni[i - 1][i] = -s1
            Gni = np.asarray(Gni)  # 单个Givens矩阵的逆
            A[0:k + 1, 0:k + 1] = np.matmul(A[0:k + 1, 0:k + 1], Gni)  # Givens矩阵的逆右乘A的分块矩阵
            GT0 = [[0 for i in range(n)] for i in range(n)]  # 每步的全块的Givens矩阵
            for i in range(0, n, 1):
                GT0[i][i] = 1.0
            GT0 = np.asarray(GT0)
            GT0[0:k + 1, 0:k + 1] = Gni
            GT = np.matmul(GT, GT0)  # 累乘每步的Givens矩阵，获得最终的Givens矩阵

        for j in range(0, k + 1, 1):
            A[j][j] += s  # 原点位移复原
    Q = np.matmul(Uk, GT)  # 获得最终Q矩阵，Q矩阵每列对应A每列本征值的本征向量

    B = [0 for i in range(2 * M + 1)]  # c存储按照从小到大排的本征值
    pp = [i for i in range(2 * M + 1)]  # 存储本征值变换的顺序
    for i in range(0, 2 * M + 1, 1):
        B[i] = A[i][i]

    for i in range(0, 2 * M + 1, 1):
        for j in range(i + 1, 2 * M + 1, 1):
            if B[j] < B[i]:  # 按照从小到大排序
                sxs = pp[i]
                pp[i] = pp[j]
                pp[j] = sxs
                cc = B[j]
                B[j] = B[i]
                B[i] = cc



    sd = 5
    for l in range(0, sd, 1):  # 按顺序对所有本征值遍历
        pxl[l][q - ql] = B[2 * l]  # 存q下前5个偶函数本征值

    for l in range(0, 11, 1):  # 按顺序对所有本征值遍历
        sxl[l][q - ql] = B[l]  # 存q下11个本征值

ta = [0 for i in range(qr - ql)]
for i in range(0, qr - ql, 1):
    ta[i] = i  # 图像的横坐标
for i in range(0, sd, 1):
    P = [0 for i in range(qr - ql)]

    for j in range(0, qr - ql, 1):
        P[j] = pxl[i][j]  # 把本征值赋予P，然后画图
        print(P[j])

    print('\n')
    plt.plot(ta, P, label=i)
    plt.legend()

show()


for i in range(0, 11, 1):
    P = [0 for i in range(qr - ql)]

    for j in range(0, qr - ql, 1):
        P[j] = sxl[i][j]  # 把本征值赋予P，然后画图
        print(P[j])

    plt.plot(ta, P, label=i)
    plt.legend()

show()










































