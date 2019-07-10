# 利用共轭梯度法求解方程组
import math
import matplotlib.pyplot as plt

N = 100  # 阶数


class matrix_value:  # 这一段定义了稀疏矩阵中的每一个值的结构
    def __init__(self, x, y, value):
        self.row = x
        self.column = y
        self.number = value


class sparse_matrix:  # 矩阵从1开始，这一段定义了矩阵的结构
    def __init__(self, x, y):
        self.row = x
        self.column = y
        self.value = []

    def assignment(self, x, y, value):  # 定义了对矩阵进行添加值的操作
        self.value.append(matrix_value(x, y, value))


def matrix_multiply(matrix, list):  # 定义了矩阵与向量相乘的函数
    if matrix.column != len(list):
        print("矩阵维度不匹配")
        exit()
    else:
        list_return = [0 for n in range(matrix.row)]
        for n in range(len(matrix.value)):
            list_return[matrix.value[n].row - 1] = list_return[matrix.value[n].row - 1] + matrix.value[n].number * list[
                matrix.value[n].column - 1]
    return list_return


def list_minus(list1, list2):  # 定义向量减法
    list_return = [0] * len(list2)
    if len(list1) != len(list2):
        print("列表长度不匹配")
        exit()
    else:
        for n in range(len(list1)):
            list_return[n] = list1[n] - list2[n]
    return list_return


def list_plus(list1, list2):  # 定义向量加法
    list_return = [0] * len(list2)
    if len(list1) != len(list2):
        print("列表长度不匹配")
        exit()
    else:
        for n in range(len(list1)):
            list_return[n] = list1[n] + list2[n]
    return list_return


def list_distance(list1, list2):  # 定义向量距离
    distance = 0
    if len(list1) != len(list2):
        print("列表长度不匹配")
        exit()
    else:
        for n in range(len(list1)):
            distance += (list1[n] - list2[n]) ** 2
        distance = distance ** 0.5
        return distance


def list_multiply(list1, list2):  # 定义向量乘法
    multiply = 0
    if len(list1) != len(list2):
        print("列表长度不匹配")
        exit()
    else:
        for n in range(len(list1)): multiply += list1[n] * list2[n]
        return multiply


def list_nummultiply(num, list2): #定义常数向量乘法
    multiply = []
    for n in range(len(list2)):
        multiply.append(num * list2[n])
    return multiply


A = sparse_matrix(N, N)  # 对A进行初始化,此稀疏矩阵不能更新值
for n in range(int(N / 2 - 1)):
    A.assignment(n + 1, N - n - 1, 0.5)
for n in range(int(N / 2 + 1), N):
    A.assignment(n + 1, N - n - 1, 0.5)
for n in range(N):
    A.assignment(n + 1, n + 1, 3)
for n in range(N - 1):
    A.assignment(n + 1, n + 2, -1)
for n in range(N - 1):
    A.assignment(n + 2, n + 1, -1)

b = [2.5]  # 初始化b
for n in range(int((N - 4) / 2)):
    b.append(1.5)
b.append(1)
b.append(1)
for n in range(int((N - 4) / 2)):
    b.append(1.5)
b.append(2.5)

x = [3 for n in range(N)]  # 初始化x
r = list_minus(b, matrix_multiply(A, x))
p = r
k = 0
listplot=[]

while k == 0 or list_distance(x, x_old) >= 1E-6:
    alpha = list_multiply(r, r) / list_multiply(p, matrix_multiply(A, p))
    x_old = x
    x = list_plus(x, list_nummultiply(alpha, p))
    r_old = r
    r = list_minus(r, list_nummultiply(alpha, matrix_multiply(A, p)))
    beta = list_multiply(r, r) / list_multiply(r_old, r_old)
    p = list_plus(r, list_nummultiply(beta, p))
    #print("迭代误差为：", math.log(list_distance(x, x_old), 10))
    listplot.append(math.log(list_distance(x, x_old), 10))
    k = k + 1

print(len(listplot))
print(listplot[0:30])
print("得到的解是：",x)
