# 本程序用来解带状对称正定矩阵
import matplotlib.pyplot as plt
N = 100
global m
m = 2


class band_matrix:  # 设立结构体
    def __init__(self, x, y):
        self.rank = x
        self.band = y
        self.value = [([0] * self.band) for n in range(self.rank)]

    def assignment(self, x, y, value):  # x,y都从0开始，赋值为value
        self.value[x][y] = value


def list_sum(i, j, start, end):  # 为后面的函数写一个a求和函数
    sum = 0
    for n in range(start, end+1):  # start 和 end都包括在内
        sum += get_a(i, n) * get_a(j, n) / get_a(n, n)
    return sum


def list_sumb(start, end, i):  # 为a和b写一个求和函数
    sum = 0
    for j in range(start, end+1):
        sum += get_a(i, j)*b[j-1]/get_a(j, j)
    return sum


def list_sumx(start, end, i):  # 为a和x写一个求和函数
    sum = 0
    for j in range(start, end + 1):
        sum += get_a(j, i) * x[j-1]
    return sum


def get_a(i, j):  # 规定i《j,从带状数组的存储空间得到存储值
    if i <= m + 1:
        return A.value[i - 1][j-i+m]
    else:
        if j < i - m:
            return 0
        else:
            return A.value[i - 1][j - i + m]


def put_a(i, j, value):  # 此写的范围仅限于在c范围内，将得到的值写入a
    A.value[i - 1][j - i + m] = value





def Cholesky():
    for i in range(2, A.rank+1):
        if i <= m + 1:# i在前面几行时
            for j in range(2, i + 1):
                r=1
                put_a(i, j, (get_a(i, j) - list_sum(i, j, r, j - 1)))
        else:
            for j in range(i-m,i+1):
                r = i - m  # 求和时在i-m开始求和
                put_a(i, j, (get_a(i, j) - list_sum(i, j, r, j - 1)))  # 将得到的a值写入





def Choleskyb():  # 对b进行操作的函数
    for i in range(2, len(b)+1):
        if i <= m + 1:
            r = 1
        else:
            r = i - m
        b[i - 1] = b[i - 1] - list_sumb(r, i-1, i)


def solve():  # 对x进行操作的函数
    x[N-1] = b[N-1]/get_a(N, N)
    for i in range(N-1, 0, -1):
        if i > N-m-1:
            t = N
        else:
            t = i + m
        x[i-1] = (b[i - 1]-list_sumx(i+1, t, i))/get_a(i, i)


global A  # 定义A为全局变量
A = band_matrix(N, 3)
A.assignment(0, 2, 5)
A.assignment(N - 1, 2, 5)
for n in range(1, N - 1):
    A.assignment(n, 2, 6)
for n in range(1, N):
    A.assignment(n, 1, 4)
for n in range(2, N):
    A.assignment(n, 0, 1)  # 初始化A矩阵
global b
b = [60]
b += (N - 2) * [120]
b += [60]  # 初始化b
global x
x = [0] * A.rank
Cholesky()  # 依次进行三个函数的操作
Choleskyb()
solve()
print(x)
