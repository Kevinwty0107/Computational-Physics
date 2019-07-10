import random
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['font.sans-serif'] = ['KaiTi']
mpl.rcParams['font.serif'] = ['KaiTi']
mpl.rcParams['axes.unicode_minus'] = False # 解决保存图像是负号'-'显示为方块的问题,或者转换负号为字符串


def E(t):  # 线偏振电场函数
    omega = 0.057
    E = -1.325 * omega * (
            0.5 * math.cos(omega * t) + 0.5 * math.cos(omega * t) * math.cos(omega * t / 4) - 1 / 8 * math.sin(
        omega * t) * math.sin(omega * t / 4))
    return E


def E_ellipse(t):  # 椭圆偏振电场函数
    omega = 0.057
    E_x = -1/math.sqrt(1.25) * 1.325 * omega * (
            0.5 * math.cos(omega * t) + 0.5 * math.cos(omega * t) * math.cos(omega * t / 4) - 1 / 8 * math.sin(
        omega * t) * math.sin(omega * t / 4))
    E_y = -0.5 / math.sqrt(1.25) * 1.325 * omega * (
            -0.5 * math.sin(omega * t) - 0.5 * math.sin(omega * t) * math.cos(omega * t / 4) - 1 / 8 * math.cos(
        omega * t) * math.sin(omega * t / 4))
    return E_x,E_y


def species(E):  # 对于固定的E选取一些v的样本
    max = 1 / 0.076 ** 2 * math.exp(-(2 / 3) / 0.076)  # max是全局的
    v_list = []  # 初始化列表
    for iter_a in range(1000):
        kesai = random.random()  # 随机生成的参考向量
        v = random.uniform(-1, 1)  # 随机选取v，在-1到1之间
        W = 1 / E ** 2 * math.exp(-(2 / 3 + v ** 2) / E)  # 隧穿几率
        if W / max > kesai:  # 判定是否选取
            v_list.append(v)
    return v_list  # 返回一系列v的样本


def line(t_0, v_y):  # 对于线偏振，给定t和v计算末状态
    v_x = 0  # x方向速度为0
    E_0 = E(t_0)
    omega = 0.057
    #  得到初始时的电场
    r_0 = (0.5/E_0)
    x = -r_0
    y = 0
    dt = 1
    t = t_0
    x_list = []
    y_list = []
    while t < 4*math.pi/omega:
        x = x + v_x * dt
        y = y + v_y * dt
        x_list.append(x)
        y_list.append(y)
        v_x = v_x + (-E(t) - x/(x**2+y**2+0.2**2)**(1.5))*dt
        v_y = v_y - y/(x**2+y**2+0.2**2)**(1.5)*dt
        t = t + dt
    return [x, y, v_x, v_y]


def ellipse(t_0, v_ver):  # 对于椭圆偏振，给定t和v计算末状态
    E_x,E_y = E_ellipse(t_0)
    E_0 = math.sqrt(E_x**2+E_y**2)
    omega = 0.057
    #  得到初始时的电场
    r_0 = (0.5 / E_0)
    x = -r_0*E_x/E_0  # 初始位置
    y = -r_0*E_y/E_0
    v_y = v_ver * E_x / E_0
    v_x = -v_ver * E_y / E_0
    dt = 1
    t = t_0
    x_list = []
    y_list = []
    while t < 4 * math.pi / omega:
        x = x + v_x * dt
        y = y + v_y * dt
        x_list.append(x)
        y_list.append(y)
        E_x, E_y = E_ellipse(t)
        v_x = v_x + (-E_x - x / (x ** 2 + y ** 2 + 0.2 ** 2) ** (1.5)) * dt
        v_y = v_y + (-E_y - y / (x ** 2 + y ** 2 + 0.2 ** 2) ** (1.5)) * dt
        t = t + dt
    return [x, y, v_x, v_y]


def p_infinite(x, y, v_x, v_y):  # 对于给定末状态，计算无穷远处动量
    p_inf = math.sqrt(v_x**2 + v_y**2 - 2/math.sqrt(x**2+y**2))
    L = x * v_y - y * v_x
    a_x = L * v_y - x / math.sqrt(x**2+y**2)
    a_y = -L * v_x - y / math.sqrt(x**2+y**2)
    p_infx = p_inf/(1+p_inf**2*L**2)*(-p_inf*L*a_y-a_x)
    p_infy = p_inf/(1+p_inf**2*L**2)*(p_inf*L*a_x-a_y)
    return p_infx, p_infy


def cac(E_1, t, p_inf1_num, function):  # function是用来计算到运动的函数
    if E_1 > 0.03:   # E很小时不会发生隧穿
        v_ver1 = species(E_1)  # 生成样品
        for i in range(len(v_ver1)):
            result_1 = function(t, v_ver1[i])
            if result_1[2] ** 2 + result_1[3] ** 2 - 2 / math.sqrt(result_1[0] ** 2 + result_1[1] ** 2) > 0:
                p_infx1, p_infy1 = p_infinite(result_1[0], result_1[1], result_1[2], result_1[3])
                if abs(p_infx1) < 1.5 and abs(p_infy1) < 1.5:
                    p_inf1_num[round((p_infy1 + 1.49) / 0.02)][round((p_infx1 + 1.49) / 0.02)] += 1
    return p_inf1_num


omega = 0.057
resolution = 10000
t_0 = -4*math.pi/omega
p_inf1_num = [[0 for i in range(150)]for j in range(150)]
p_inf2_num = [[0 for i in range(150)]for j in range(150)]
for n in range(resolution):  # 对于时间t做循环，计算末状态
    t = t_0 + 8*math.pi/omega * n / resolution
    print('已完成计算：', n/resolution*100, '%')
    E_1 = E(t)
    E_2x, E_2y = E_ellipse(t)
    E_2 = math.sqrt(E_2x**2+E_2y**2)  # 计算E
    p_inf1_num = cac(E_1, t, p_inf1_num, line)
    p_inf2_num = cac(E_2, t, p_inf2_num, ellipse)

p_x_list = []
p_y_list = []
for t in range(150):  # 初始化x和y坐标
    p_x_list.append(-1.49 + t * 0.02)
    p_y_list.append(-1.49 + t * 0.02)

p_inf1_num_x = [0] * len(p_x_list)  # 分到得到x和y分布的频数
p_inf1_num_y = [0] * len(p_x_list)
p_inf2_num_x = [0] * len(p_x_list)
p_inf2_num_y = [0] * len(p_x_list)
for i in range(len(p_x_list)):
    for j in range(len(p_y_list)):
        p_inf1_num_x[j] += p_inf1_num[i][j]
        p_inf1_num_y[i] += p_inf1_num[i][j]
        p_inf2_num_x[j] += p_inf2_num[i][j]
        p_inf2_num_y[i] += p_inf2_num[i][j]
line_x_max1 = p_inf1_num_x.index(max(p_inf1_num_x[0:int(len(p_inf1_num_x)/2)]))
line_x_max2 = p_inf1_num_x.index(max(p_inf1_num_x[int(len(p_inf1_num_x)/2):len(p_inf1_num_x)]))
line_y_max = p_inf1_num_y.index(max(p_inf1_num_y))
ellipse_x_max = p_inf2_num_x.index(max(p_inf2_num_x))
ellipse_y_max1 = p_inf2_num_y.index(max(p_inf2_num_y[0:int(len(p_inf2_num_y)/2)]))
ellipse_y_max2 = p_inf2_num_y.index(max(p_inf2_num_y[int(len(p_inf2_num_y)/2):len(p_inf2_num_y)]))
print('线偏振下x方向动量峰值为：','p_x1=',line_x_max1*0.02-1.49,'p_x2=',line_x_max2*0.02-1.49)
print('线偏振下y方向动量峰值为：','p_y=',line_y_max*0.02-1.49)
print('椭圆偏振下x方向动量峰值为：','p_x=',ellipse_x_max*0.02-1.49)
print('椭圆偏振下y方向动量峰值为：','p_y1=',ellipse_y_max1*0.02-1.49,'p_y2=',ellipse_y_max2*0.02-1.49)

fig1 = plt.figure(1)
surf1 = plt.contourf(p_x_list, p_y_list, p_inf1_num)
fig1.colorbar(surf1)
plt.title('线偏振光下动量分布')
plt.xlabel('px')
plt.ylabel('py')
plt.show()

fig2 = plt.figure(2)
surf2 = plt.contourf(p_x_list, p_y_list, p_inf2_num)
fig2.colorbar(surf2)
plt.title('椭圆偏振光下动量分布')
plt.xlabel('px')
plt.ylabel('py')
plt.show()

fig3 = plt.figure(3)
plt.plot(p_x_list, p_inf1_num_x, label = 'px')
plt.plot(p_y_list, p_inf1_num_y, label = 'py')
plt.xlabel('p')
plt.ylabel('frequency')
plt.title('line_p')
plt.legend()
plt.show()


fig4 = plt.figure(4)
plt.plot(p_x_list, p_inf2_num_x,label = 'px')
plt.plot(p_y_list, p_inf2_num_y,label = 'py')
plt.xlabel('p')
plt.ylabel('frequency')
plt.title('ellipse_px')
plt.legend()
plt.show()





