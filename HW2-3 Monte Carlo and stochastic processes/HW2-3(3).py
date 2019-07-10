import math
import random
import matplotlib.pyplot as plt


def walk(x, y):
    r = []
    for i in range(1, 101):
        theta = random.uniform(0, 2 * math.pi)  # 随机输出角度
        x = x + math.cos(theta)
        y = y + math.sin(theta)
        R = math.sqrt(x * x + y * y)  # 求距离
        r.append(R)  # r是储存距离的list
        print('第', i, '步行走后：', 'x坐标为：', x, 'y坐标为：', y, '与原点距离：', R, '\n')
    return r


def linefit(x , y):  # 最小二乘法函数
    N = float(len(x))
    sx,sy,sxx,syy,sxy=0,0,0,0,0
    for i in range(0,int(N)):
        sx += x[i]
        sy += y[i]
        sxx += x[i]*x[i]
        syy += y[i]*y[i]
        sxy += x[i]*y[i]
    a = (sy*sx/N -sxy)/( sx*sx/N -sxx)
    b = (sy - a*sx)/N
    r = abs(sy*sx/N-sxy)/math.sqrt((sxx-sx*sx/N)*(syy-sy*sy/N))
    return a,b,r


x = 0
y = 0  # 初始在原点
r = [0] * 100
step = 5000
for j in range(step):
    a = walk(x, y)
    r = [r[i] + a[i] for i in range(min(len(r), len(a)))]
r = [r[i] / step for i in range(len(r))]
n = [math.sqrt(i) for i in range(1, 101)]
N = [i for i in range(1, 101)]
print(r)
plt.scatter(N, r)
plt.xlabel('N')
plt.ylabel('E(r)')
plt.show()
aa = linefit(n, r)
print('回归系数为：',aa[0])