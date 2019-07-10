import matplotlib.pyplot as plt
  # 以二阶导数表示三次样条插值函数


def S(x, xlist, h, ylist, M):  # 定义三次样条函数，输入坐标，得到y值
    for i in range(len(xlist)-1):
        if x<=xlist[i+1] and x>=xlist[i]:
            j = i
    s = (xlist[j + 1] - x) ** 3 * M[j] / 6 / h[j] \
    + (x - xlist[j]) ** 3 * M[j + 1] / 6 / h[j]\
    + (ylist[j] - M[j] * h[j] ** 2 /6) *(xlist[j+1] - x) / h[j]\
    + (ylist[j + 1] - M[j + 1] * h[j] ** 2 / 6) *(x - xlist[j]) / h[j]
    return s


def chase(n, a, b, c, f):  # 定义追逐法解方程的函数
    for i in range(1, n):
        m = a[i]/b[i-1]
        b[i] = b[i] - m * c[i-1]
        f[i] = f[i] - m * f[i-1]
    x = [0] * n
    x[n-1] = f[n-1]/b[n-1]
    for i in range(n-2, -1, -1):
        x[i] = (f[i] - c[i] * x[i + 1])/b[i]
    return x


x = [0, 3, 5, 7, 9, 11, 12, 13, 14, 15]  # 初始化后面要用到的变量
y = [0, 1.2, 1.7, 2.0, 2.1, 2.0, 1.8, 1.2, 1.0, 1.6]
h = [x[i+1] - x[i] for i in range(len(x)-1)]
miu = [(h[j])/(h[j] + h[j + 1]) for j in range(len(h) - 1)]
lam = [(h[j + 1])/(h[j] + h[j + 1]) for j in range(len(h) - 1)]
d = [6 / (h[j] + h[j + 1]) * ((y[j] - y[j + 1]) / h[j] + (y[j + 2] - y[j + 1]) / h[j + 1]) for j in range(len(h)-1)]
M = chase(len(h)-1, miu, [2]*(len(h)-1), lam, d)
M = [0] + M + [0]  # M的初末添加0
aa = [0.1 * i for i in range(151)]
bb = []
for i in range(len(aa)):
    bb.append(S(aa[i], x, h, y, M))
plt.plot(aa, bb)  # 绘图
plt.xlabel("x")
plt.ylabel("S(x)")
plt.title('Wing outline')
plt.show()
print(bb)
