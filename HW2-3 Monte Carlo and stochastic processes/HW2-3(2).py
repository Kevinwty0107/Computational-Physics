# 本程序用来解决第三题第二问
import math
import random

x = 0
y = 0  #初始在原点
for i in range(1, 51):
    a = random.randint(1, 4)
    if a==1:  # 1代表向左
        x = x - 1
    elif a==2:  # 2代表向右
        x = x + 1
    elif a==3:  # 3代表向上
        y = y + 1
    elif a==4:  # 4代表向下
        y = y - 1
    R = math.sqrt(x*x+y*y)  # 求距离
    print('第', i,'步行走后：', 'x坐标为：', x, 'y坐标为：', y, '与原点距离：', R, '\n')
