# coding=UTF-8
# Guass求积公式
import math
import matplotlib.pyplot as plt
from pylab import *
mpl.rcParams['font.sans-serif'] = ['SimHei']


def GEM(A, b, x):  # 高斯消元法。Ax=b
    n = len(A)  # A的行数
    m = len(A[0])  # A的列数
    B = [[0 for i in range(m)] for j in range(n)]  # B为增广矩阵
    for i in range(n):
        for j in range(n):
            # print(B)
            B[i][j] = A[i][j]  # 为增广矩阵赋值
        B[i].append(b[i])

    for i in range(n - 1):
        # 找主元
        PivotColumn = []
        for iPivot in range(i, n):
            PivotColumn.append(abs(B[iPivot][i]))
        if not max(PivotColumn):
            print('矩阵奇异')
            return

        PivotIndex = PivotColumn.index(max(PivotColumn))
        temp1 = B[i]
        B[i] = B[i + PivotIndex]
        B[i + PivotIndex] = temp1
        # 消元
        for iothers in range(i + 1, n):
            Coeff = B[iothers][i] / B[i][i]
            B[iothers][i] = 0
            for j in range(i + 1, m + 1):
                B[iothers][j] = B[iothers][j] - Coeff * B[i][j]

    # 反代
    for i in range(m - 1, -1, -1):
        x[i] = B[i][m]
        for j in range(m - 1, i, -1):
            x[i] = x[i] - B[i][j] * x[j]
        x[i] = x[i] / B[i][i]
    return x


def multiple(list1, list2, list_result): # 向量内积
    for i in range(len(list_result)):
        for x in range(i + 1):
            list[i] += list[x] * list[i - x]


def f1(n):   # 定义与根号x在0到1区间的积分结果
    return 1/(n + 3/2)


def f2(n):     # 定义与(1+x^2.py)在-1到1区间的积分结果
    return (1/(n + 1) + 1/(n + 3)) * (1 ** (n + 1) - (-1) ** (n + 1))


def equation_2(list, list_result):  # 解一元二次方程
    a = list[2]
    b = list[1]
    c = list[0]
    delta = b * b - 4 * a * c
    if delta < 0:
        print('No result.')
        exit()
    x1 = (-b - math.sqrt(delta)) / (2 * a)
    x2 = (-b + math.sqrt(delta)) / (2 * a)
    list_result[0] = x1
    list_result[1] = x2


def det(a11, a12, a21, a22):       # 行列式
    return a11 * a22 - a21 * a12


def equa(a11, a12, a21, a22, b1, b2, list_result):   # 解二元一次方程组
    s = det(a11, a12, a21, a22)
    x1 = det(b1, a12, b2, a22) / s
    x2 = det(a11, b1, a21, b2) / s
    list_result[0] = x1
    list_result[1] = x2


# 1
a1 = -f1(1) / f1(0)
b1 = -f1(2)
a11 = f1(1)
a12 = f1(0)
b2 = - f1(3)
a21 = f1(2)
a22 = f1(1)
list_1 = [0, 0]
equa(a11, a12, a21, a22, b1, b2, list_1)
list_X = [0, 0]
list = [list_1[1], list_1[0], 1]
equation_2(list, list_X)
list_A = [0, 0]
equa(1, 1, list_X[0], list_X[1], f1(0), f1(1), list_A)
print('方程1：\n','[x1,x2] =',list_X,'\n','[A1,A2] =',list_A)
# 2.py
#a2 = -f2(1) / f2(0)
b1 = -f2(2)
a11 = f2(1)
a12 = f2(0)
b2 = -f2(3)
a21 = f2(2)
a22 = f2(1)
list_2 = [0, 0]
equa(a11, a12, a21, a22, b1, b2, list_2)
list_X = [0, 0]
list = [list_2[1], list_2[0], 1]
equation_2(list, list_X)
list_A = [0, 0]
equa(1, 1, list_X[0], list_X[1], f2(0), f2(1), list_A)
print('方程2：\n','[x1,x2] =', list_X, '\n', '[A1,A2] =', list_A)
# 3
a = a1
b = list_1[0]
c = list_1[1]
A1 = [f1(2), f1(1), f1(0)]
A2 = [f1(3), f1(2), f1(1)]
A3 = [f1(4), f1(3), f1(2)]
A = [A1, A2, A3]
b = [-f1(3), -f1(4), -f1(5)]
x = [0, 0, 0]
GEM(A, b, x)
d = x[0]
e = x[1]
f = x[2]
p = e - d ** 2 / 3
q = f - d ** 3 / 27 - d / 3 * p
r = math.sqrt(-(p/3)**3)
theta = 1/3 * math.acos(-q /2/r)
u1 = 2 * math.pow(r, 1/3) * math.cos(theta)
u2 = 2 * math.pow(r, 1/3) * math.cos(theta + 2*math.pi/3)
u3 = 2 * math.pow(r, 1/3) * math.cos(theta + 4*math.pi/3)
x = [u1 - d/3, u2 - d/3, u3 - d/3]
A_1 = [1, 1, 1]
A_2 = [x[0], x[1], x[2]]
A_3 = [x[0]**2, x[1]**2, x[2]**2]
A = [A_1, A_2, A_3]
b = [f1(0), f1(1), f1(2)]
x_new =[0, 0, 0]
GEM(A, b, x_new)
print('方程3：\n','[x1,x2,x3] =',x,'\n','[A1,A2,A3] =',x_new)
