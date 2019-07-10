import math
import matplotlib.pyplot as plt


def E(t):
    omega = 0.057
    E = -1.325 * omega * (
            0.5 * math.cos(omega * t) + 0.5 * math.cos(omega * t) * math.cos(omega * t / 4) - 1 / 8 * math.sin(
        omega * t) * math.sin(omega * t / 4))
    return E


def E_ellipse(t):
    omega = 0.057
    E_x = -1/math.sqrt(1.25) * 1.325 * omega * (
            0.5 * math.cos(omega * t) + 0.5 * math.cos(omega * t) * math.cos(omega * t / 4) - 1 / 8 * math.sin(
        omega * t) * math.sin(omega * t / 4))
    E_y = -0.5 / math.sqrt(1.25) * 1.325 * omega * (
            -0.5 * math.sin(omega * t) - 0.5 * math.sin(omega * t) * math.cos(omega * t / 4) - 1 / 8 * math.sin(
        omega * t) * math.cos(omega * t / 4))
    return E_x,E_y


def line(t_0, v_y):
    v_x = 0  # x方向速度为0
    E_0 = E(t_0)
    omega = 0.057
    #  得到初始时的电场
    r_0 = abs(0.5/E_0)
    x = -r_0
    y = 0
    dt = 0.001
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
    r_final = math.sqrt(x**2+y**2)
    v_final = math.sqrt(v_x**2+v_y**2)
    print('线偏振光作用下的末态','x=', x, 'y=', y, 'v_x=', v_x, 'v_y=', v_y)
    plt.plot(x_list, y_list)
    plt.title('line')
    #plt.legend()
    plt.show()


def ellipse(t_0, v_ver):
    E_x,E_y = E_ellipse(t_0)
    E_0 = math.sqrt(E_x**2+E_y**2)
    omega = 0.057
    #  得到初始时的电场
    r_0 = abs(0.5 / E_0)
    x = -r_0*E_x/E_0  # 初始位置
    y = -r_0*E_y/E_0
    v_y = v_ver * E_x/E_0
    v_x = -v_ver * E_y/E_0
    dt = 0.001
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
    r_final = math.sqrt(x ** 2 + y ** 2)
    v_final = math.sqrt(v_x ** 2 + v_y ** 2)
    print('椭圆偏振光作用下的末态','x=', x, 'y=', y, 'v_x=', v_x, 'v_y=', v_y)
    plt.plot(x_list, y_list)
    plt.title('ellipse')
    #plt.legend()
    plt.show()


t_1 = float(input('输入线偏振光初始时间: '))
v_1 = float(input('输入线偏振光垂直速度: '))
t_2 = float(input('输入椭圆偏振光初始时间: '))
v_2 = float(input('输入椭圆偏振光垂直速度: '))
line(t_1, v_1)
ellipse(t_2, v_2)