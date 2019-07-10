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
phi = [0 for i in range(2 * M + 2)]
for i in range(1, 2 * M + 2, 1):
    phi[i] = (2 * math.pi * i) / (2 * M + 1)
ql = 10
qr = 11
mr = 11
pxl = [[0 for i in range(qr - ql)] for i in range(n)]
bzxl = [[0 for i in range(qr - ql)] for i in range(n)]


def sign1(a):
    if a > 0: return 1
    if a < 0: return -1


def fs2(a, m):  # 定义2范数函数
    c = 0
    for i in range(0, m, 1):
        c += a[i] ** 2
    c = math.sqrt(c)
    return c


def jdz(a):
    if a >= 0:
        return a
    else:
        return -a


for q in range(ql, qr, 1):
    A = [[0 for i in range(2 * M + 2)] for i in range(2 * M + 2)]
    for j in range(1, 2 * M + 2, 1):
        for k in range(1, 2 * M + 2, 1):
            if j == k:
                A[j][j] = M * (M + 1) / 3 + 2 * q * math.cos(2 * phi[j])
            else:
                A[j][k] = ((-1) ** (j - k)) * math.cos(math.pi * (j - k) / (2 * M + 1)) / (
                2 * math.sin(math.pi * (j - k) / (2 * M + 1)) ** 2)
    A = np.asarray(A)

    Uk = [[0 for i in range(2 * M + 1)] for i in range(2 * M + 1)]
    for i in range(0, 2 * M + 1, 1):
        Uk[i][i] = 1.0
    Uk = np.asarray(Uk)

    for k in range(1, 2 * M, 1):
        x = A[k + 1:2 * M + 2, k]
        e1 = [0 for i in range(2 * M - k + 1)]
        e1[0] = 1.0
        e1 = np.asarray(e1)
        v = sign1(x[0]) * fs2(x, 2 * M - k + 1) * e1 + x
        v = v / fs2(v, 2 * M - k + 1)
        vvT = [[0 for i in range(2 * M - k + 1)] for i in range(2 * M - k + 1)]
        Uk0 = [[0 for i in range(2 * M + 2)] for i in range(2 * M + 2)]
        for i in range(1, 2 * M + 2, 1):
            Uk0[i][i] = 1.0
        Uk0 = np.asarray(Uk0)
        for i in range(0, 2 * M - k + 1, 1):
            for j in range(0, 2 * M - k + 1, 1):
                vvT[i][j] = v[i] * v[j]
        vvT = np.asarray(vvT)
        Uk0[k + 1:2 * M + 2, k + 1:2 * M + 2] = Uk0[k + 1:2 * M + 2, k + 1:2 * M + 2] - 2 * vvT
        Uk0 = Uk0[1:2 * M + 2, 1:2 * M + 2]
        Uk = np.matmul(Uk, Uk0)
        A[k + 1:2 * M + 2, k:2 * M + 2] = A[k + 1:2 * M + 2, k:2 * M + 2] - 2 * np.matmul(vvT, A[k + 1:2 * M + 2,
                                                                                               k:2 * M + 2])
        A[1:2 * M + 2, k + 1:2 * M + 2] = A[1:2 * M + 2, k + 1:2 * M + 2] - 2 * np.matmul(
            A[1:2 * M + 2, k + 1:2 * M + 2], vvT)
    A = np.asarray(A[1:2 * M + 2, 1:2 * M + 2])
    for i in range(0, 2 * M + 1, 1):
        for j in range(0, 2 * M + 1, 1):
            if jdz(A[i][j]) < 10 ** (-10):
                A[i][j] = 0


    GT = [[0 for i in range(n)] for i in range(n)]
    for i in range(0, n, 1):
        GT[i][i] = 1.0
    GT = np.asarray(GT)

    k = 2 * M
    while k > 0:
        if jdz(A[k][k - 1]) < 10 ** (-12):
            # A[k][k-1]=0
            k = k - 1
            # print(k)
        s = A[k][k]
        for j in range(0, k + 1, 1):
            A[j][j] -= s
        for i in range(1, k + 1, 1):
            a = math.sqrt(A[i - 1][i - 1] ** 2 + A[i][i - 1] ** 2)
            c = A[i - 1][i - 1] / a
            s1 = A[i][i - 1] / a
            G = [[0 for i in range(k + 1)] for i in range(k + 1)]
            for j in range(0, k + 1, 1):
                G[j][j] = 1.0
            G[i][i] = G[i - 1][i - 1] = A[i - 1][i - 1] / a
            G[i][i - 1] = -s1
            G[i - 1][i] = s1
            G = np.asarray(G)
            A[0:k + 1, 0:k + 1] = np.matmul(G, A[0:k + 1, 0:k + 1])
            Gni = G
            Gni[i][i - 1] = s1
            Gni[i - 1][i] = -s1
            Gni = np.asarray(Gni)
            A[0:k + 1, 0:k + 1] = np.matmul(A[0:k + 1, 0:k + 1], Gni)
            GT0 = [[0 for i in range(n)] for i in range(n)]
            for i in range(0, n, 1):
                GT0[i][i] = 1.0
            GT0 = np.asarray(GT0)
            GT0[0:k + 1, 0:k + 1] = Gni
            GT = np.matmul(GT, GT0)
            # print(GT)
        for j in range(0, k + 1, 1):
            A[j][j] += s
    Q = np.matmul(Uk, GT)



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





    sd = 0
    for l in range(0, n, 1):  # 按顺序对所有本征值遍历

        maxdf = 0
        maxddf = 0
        #print(B[l])
        #print(l)
        for i in range(1, M + 1, 1):
            df = 0
            ddf = 0
            for k in range(0, n, 1):
                df += Q[k][pp[l]] * math.sin(i * phi[k + 1])  #
                ddf += Q[k][pp[l]] * math.cos(i * phi[k + 1])

            if maxdf < jdz(df):
                maxdf = jdz(df)
            if maxddf < jdz(ddf):
                maxddf = jdz(ddf)

        if maxddf == 0:
            continue

        ggg = maxdf / maxddf
        if jdz(ggg) < 0.01:  # 对任何n，此时sin函数前面系数均为0，即为偶函数
            pxl[sd][q - ql] = B[l]  # 存本征值
            bzxl[sd][q - ql] = pp[l]
            sd += 1
            #print("haha")
        if sd == 6:
            break
            # print(sd)




jiaodu = [0.0 for i in range(101)]  # 0到90度取100个数
for i in range(0, 101, 1):
    jiaodu[i] = i * math.pi / 200  # 每个数赋值
jiaodu = np.asarray(jiaodu)
for i in range(0, 6, 1):
    wz = bzxl[i][0]  # wz就是本征向量在Q的列数
    y = [0.0 for j in range(101)]
    y = np.asarray(y)  # y就是每个本征值下的解函数
    for nn in range(1, M + 1, 1):
        osx = 0.0
        jsx = 0.0
        for k in range(0, n, 1):
            osx += Q[k][wz] * math.cos(nn * phi[k + 1])  # 用来计算每个cos函数前的系数
            jsx += Q[k][wz] * math.sin(nn * phi[k + 1])  # 用来计算每个sin函数前的系数
        y += osx * np.cos(nn * jiaodu) + jsx * np.sin(nn * jiaodu)  # 角度相关赋予
        print(osx)
        print(jsx)
    y *= 2
    for k in range(0, n, 1):
        y += Q[k][wz]  # 加上常数项
    y=y/(2*math.pi*n)
    if y[1] < 0: y = -y
    plt.plot(jiaodu, y, label=pxl[i][0])
    plt.legend()
show()


for q in range(ql, qr, 1):
    A = [[0 for i in range(2 * M + 2)] for i in range(2 * M + 2)]
    for j in range(1, 2 * M + 2, 1):
        for k in range(1, 2 * M + 2, 1):
            if j == k:
                A[j][j] = M * (M + 1) / 3 + 2 * q * math.cos(2 * phi[j])
            else:
                A[j][k] = ((-1) ** (j - k)) * math.cos(math.pi * (j - k) / (2 * M + 1)) / (
                2 * math.sin(math.pi * (j - k) / (2 * M + 1)) ** 2)
    A = np.asarray(A)

for i in range(n):
    print(np.matmul(A[1:n+1,1:n+1],Q[:,pp[i]])-B[i]*Q[:,pp[i]])









