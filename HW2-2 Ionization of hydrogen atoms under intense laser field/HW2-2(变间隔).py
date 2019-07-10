import math
import random
import matplotlib.pyplot as plt


I_0 = 5E13  # 首先是输入变量
lam = 3200

omega = 45.5633525316/lam
E_0 = math.sqrt(I_0/3.5094448314E16)
A_0 = E_0/omega


def csin(com):  # 复数三角函数
    return (complex(math.exp(-com.imag),0)*complex(math.cos(com.real),math.sin(com.real))-
            complex(math.exp(com.imag), 0) * complex(math.cos(-com.real), math.sin(-com.real)))\
            / complex(0, 2)


def ccos(com):  # 复数三角函数
    return (complex(math.exp(-com.imag), 0) * complex(math.cos(com.real), math.sin(com.real)) +
            complex(math.exp(com.imag), 0) * complex(math.cos(-com.real), math.sin(-com.real))) \
           / complex(2, 0)


def function(com, success, p):  # 抛除已有结果的方程
    global A_0
    global omega
    if len(success) != 0:
        result = 0.5*(A_0*csin(omega*com)*(csin(omega*com/4))**2+p)**2 + 0.5 #方程
        for t in range(len(success)):  # 除掉已有的结果
            result = result/(com - success[t])
        return result
    else:
        return 0.5*(A_0*csin(omega*com)*(csin(omega*com/4))**2+p)**2 + 0.5 #方程


def function_origin(com, p):  # 不抛除已有结果的方程
    global A_0
    global omega
    return 0.5*(A_0 * csin(omega * com) * (csin(omega * com / 4)) ** 2 + p) ** 2 + 0.5


def solve_chase(com_a, com_b, success, error, p):  # 抛除已有结果的解步骤
    iter_a = 1  # 迭代次数
    fun = 0
    while (abs(fun) > error) or iter_a == 1:
        change = com_b
        fun = function(com_b, success, p)
        com_b = com_a - (com_a - com_b) / (function(com_a, success, p) - function(com_b, success, p)) * function(com_a,
                                                                                                  success, p)
        com_a = change
        iter_a += 1
        if iter_a > 1000 or abs(com_b.imag) > 100 or (function(com_a, success, p) - function(com_b, success, p)) == 0:
            break  # 防止除数为0与溢出
    return com_b


def solve_chase_origin(com_b0, error, p):  # 不抛除已有结果的解步骤
    iter_a = 1  # 迭代次数
    fun = 0
    com_b = com_b0
    com_a = com_b + complex(random.uniform(0, 1), random.uniform(0, 1))
    while (abs(fun) > error) or iter_a == 1:
        change = com_b
        while (function_origin(com_a, p) - function_origin(com_b, p)) == 0 or iter_a > 1000:  # 防止出现为0值
            com_b = com_b0
            com_a = com_b0 + complex(random.uniform(0, 1), random.uniform(0, 1))
            iter_a = 1
        com_b = com_a - (com_a - com_b) / (function_origin(com_a, p) - function_origin(com_b, p)) * (function_origin(com_a, p))
        com_a = change
        fun = function_origin(com_b, p)
        iter_a += 1
    return com_b


def solve_t(p):  # 解解向量的函数
    goal = 0
    success = []
    print('p为', p, '时解为：')
    while goal < 6:
        error = 1E-14  # 误差
        com_a = complex(random.uniform(0, 4*math.pi/omega), random.uniform(0, 100))  # 初始化
        com_b = complex(random.uniform(0, 4*math.pi/omega), random.uniform(0, 100))
        com_b = solve_chase(com_a, com_b, success, error, p)  # 解出x值
        if com_b.imag > 0 and com_b.real < 4*math.pi/omega and com_b.real>0 and com_b.imag < 100 and\
                                abs(function(com_b, success, p)) < error:
            com_b = solve_chase_origin(com_b, error, p)  #在原有解上加一微扰
            print('解为：', com_b)
            success.append(com_b)  # 添加到序列里
            goal = goal+1
    print('\n')
    return success


def S_two_derivative(com, p):  # 二阶导数
    global A_0
    global omega
    return (A_0 * csin(omega * com) * (csin(omega * com / 4)) ** 2 + p)*(omega * A_0 * ccos(omega * com) *
            (csin(omega * com / 4)) ** 2 + omega/4*A_0*csin(omega * com)*csin(omega * com / 2))


def S_sum(com, p, n):  # 求和计算S积分
    global A_0
    global omega
    h1 = com.real/n
    h2 = com.imag/n
    sum_real = 0
    sum_imag = 0
    for i in range(n):
        sum_real = sum_real + function_origin(i*h1, p) * h1
        sum_imag = sum_imag + function_origin(complex(com.real, i*h2), p) * complex(0, h2)
    sum = sum_real + sum_imag
    return sum


def q(t, p):  # 积分用到的q函数
    global omega
    global A_0
    return A_0*csin(omega*t)*(csin(omega*t/4))**2+p


def E(t):   # 积分用到的E函数
    global omega
    global A_0
    return -(omega * A_0 * ccos(omega * t) *
            (csin(omega * t / 4)) ** 2 + omega/4*A_0*csin(omega * t)*csin(omega * t / 2))


E_k = []  # 能量向量
M_p_list1 = []  # 求和的M向量
M_p_list2 = []  # 积分的M向量
E_k_now = 0  # 现在的能量
n = 10000  # 分辨率，即积分区间份数
while E_k_now < 2:
    if E_k_now < 0.2:
        dE = 0.001
    elif E_k_now> 0.2 and E_k_now <1:
        dE = 0.003
    else:
        dE = 0.015
    E_k_now = E_k_now + dE  # 现在的能量
    print('E=', E_k_now)
    E_k.append(E_k_now)
    p = math.sqrt(2 * E_k_now)  # 现在的动量
    # 首先进行求和的积分
    t_list = solve_t(p)  # 解所有的t
    M_p = 0
    for j in range(len(t_list)):
        S_num1 = S_sum(t_list[j], p, n)
        M_p = M_p + math.exp(-S_num1.imag)*complex(math.cos(S_num1.real), math.sin(S_num1.real)) / S_two_derivative(t_list[j], p)
    M_p = (-M_p) * (2*0.5) ** (5/4) / math.sqrt(2)
    M_p_list1.append(abs(M_p)**2)
    # 然后进行积分的求解
    M_p = 0
    S_num2 = 0
    for j in range(1, 1000):
        t = 4*math.pi/omega/1000 * j
        S_num2 += (function_origin(t, p) * 4*math.pi/omega/1000).real  # S_num2 是实数了
        M_p = M_p + (q(t, p) * E(t)) / (q(t, p)**2 + 1)**3 * (math.cos(S_num2)+complex(0, math.sin(S_num2))) * 4*math.pi/omega/1000
    M_p = M_p * 2**(7/2) * (2*0.5)**(5/4) / math.pi
    M_p_list2.append(abs(M_p) ** 2)
plt.plot(E_k, M_p_list1, label='求和的结果')
plt.plot(E_k, M_p_list2, label='积分的结果')
plt.xlabel('E_k')
plt.ylabel('电离几率')
plt.legend()
plt.show()



