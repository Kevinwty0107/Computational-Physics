# 本程序用来计算球谐函数
import math
import sys
import numpy as np

import pandas as pd
sys.setrecursionlimit(1000000)

def scinum(a):  # a是小数项，b是指数项
    b = 0
    if abs(a) >= 10:
        while abs(a) >= 10:
            a = a/10
            b += 1
    elif abs(a) < 1:
        while abs(a) < 1:
            a = a*10
            b -= 1
    return [a, b]


def scinum_adjust(a):  # 用以稍作调整为规范科学计数法
    if abs(a[0]) >= 10:
        while abs(a[0]) >= 10:
            a[0] = a[0] / 10
            a[1] += 1
    elif abs(a[0]) < 1:
        while abs(a[0]) < 1:
            a[0] = a[0] * 10
            a[1] -= 1
    return a

def scinummultiply(a, b):  # 科学计数法的乘法
    c = [0] * 2
    c[0] = a[0] * b[0]
    c[1] = a[1] + b[1]
    return scinum_adjust(c)


def scinumdivide(a, b):  # 科学计数法的除法
    c = [0] * 2
    c[0] = a[0] / b[0]
    c[1] = a[1] - b[1]
    return scinum_adjust(c)


def scinumplus(a, b):  # 科学计数法的加法
    c = [0] * 2
    if a[1] > b[1]:
        c[1] = b[1]
        c[0] = a[0] + b[0] / (10 ** (a[1]-b[1]))
    else:
        c[1] = a[1]
        c[0] = b[0] + a[0] / (10 ** (b[1] - a[1]))
    return scinum_adjust(c)


def scinumminus(a, b):  # 科学计数法的减法
    c = [0] * 2
    if a[1] > b[1]:
        c[1] = a[1]
        c[0] = a[0] - b[0] / (10 ** (a[1]-b[1]))
    else:
        c[1] = b[1]
        c[0] = a[0] / 10 ** (b[1] - a[1]) - b[0]
    return scinum_adjust(c)


def scinumsqrt(a):  #科学计数法的开方
    if (a[1] % 2) == 0:
        a[1] = a[1]/2
        a[0] = math.sqrt(a[0])
    else:
        a[1] = (a[1] - 1)/2
        a[0] = math.sqrt((a[0]*10))
    return a


def scinumfactorial(t):  # 输入是整数，输出科学计数法的乘方
    a = [1, 0]
    if t==0:
        return a
    else:
        while t != 1:
            a = scinummultiply(a, scinum(t))
            t -= 1
        return a


def Lnax(n,alpha,x):
    t = n
    if n == 0:
        return 1
    else:
        while t > 0:
            a = (t+alpha)*Lnax(t-1,alpha,x) - x * Lnax(t-1,alpha+1,x)
            a = a / t
            t -= 1
            #print(t)
            return a




def pll(l, xx):  # 定义ll阶勒让德函数，返回值是科学计数法,输入值是浮点数
    if l == 0:
        return [1, 0]
    else:
        a = scinummultiply(scinum(-(2*l-1)),
                           scinummultiply(scinumsqrt(scinum(1 - xx ** 2)), pll(l-1, xx)))
        return a


def legendre(l, m, x):  # 定义勒让德函数
    if m < 0:  # 如果m为负值
        m = -m
        a = scinummultiply(scinum(math.pow(-1, m)),
                        scinummultiply(scinumdivide(scinumfactorial(l - m), scinumfactorial(l + m)), legendre(l, m, x)))
        return a
    elif l == m:  # 如果l=m返回pll值
        a = pll(l, x)
        return a
    else:  # 否则利用递推公式
        t = m + 1
        mm = pll(m, x)
        mm2 = scinummultiply(scinum(x * (2 * m + 1)), pll(m, x))
        while l > t:
            t += 1
            mm_exchange = mm2
            mm2 = scinumdivide(scinumminus(scinummultiply(scinum((2 * t - 1) * x), mm2),
                                           scinummultiply(scinum(m + t - 1), mm)),
                               scinum(t - m))
            mm = mm_exchange
        return mm2

# Lnax
for n in [3, 10 ,30]:
    for alpha in [2, 20, 40]:
        for x in [0.001, 1, 100]:
            print("n=", n, "alpha=", alpha, "x=", x, "Lnax=", Lnax(n,alpha,x))


#Ylm,theta,phi

#input parameters
#print("L和M值自动取整，请输入大于0的L和M值，且L>=M")  # 输入L和M的值
#L = int(input("Enter the input SphericalHarmonics parameter L: "))
#while L < 0:
#    print("您输入的数不在范围内，请重新输入：\n")
#    L = int(input("Enter the input SphericalHarmonics parameter L: "))
#M = int(input("Enter the input SphericalHarmonics parameter M: "))
#while abs(M) > L:
#    print("您输入的数不在范围内，请重新输入：\n")
#    L = int(input("Enter the input SphericalHarmonics parameter M: "))
#theta = float(input("Enter the input SphericalHarmonics parameter theta(角度制): ")) * math.pi / 180  # 实际的储存用弧度制
#phi = float(input("Enter the input SphericalHarmonics parameter phi(角度制): ")) * math.pi / 180  # 实际的储存用弧度制


#pre-defined parameters




for theta in [math.pi /1000, 3*math.pi/10 , 501* math.pi/1000 ]:
    for L in [100, 500 ,1000]:
        for M in [1, L/100, L/10, L-1 ]:
            phi = math.pi / 5

            x = math.cos(theta)
            Ya = scinumsqrt(scinumdivide(scinummultiply(scinum(2 * L + 1), scinumfactorial(L - M)),
                                         scinummultiply(scinum(4 * math.pi), scinumfactorial(L + M))))  # 计算第二部分
            Yb_real = scinum(math.cos(M * phi))  # 计算第三部分的实部和虚部
            Yb_imag = scinum(math.sin(M * phi))
            Yc = legendre(L, M, x)  # 计算第一部分
            SH_real = scinummultiply(Ya, scinummultiply(Yb_real, Yc))  # 计算球谐函数的实部和虚部
            SH_imag = scinummultiply(Ya, scinummultiply(Yb_imag, Yc))
            print('L=', str(L),'M=', str(M), 'theta=', str(theta), 'phi=', str(phi),'球谐函数值为',SH_real[0],'* 10 ^',
                  SH_real[1],'+',SH_imag[0],'* 10 ^',SH_imag[1],'i\n')









# PHI n=2,l=1,m r, theta , phi
#M = +1

L = 1
M = 1
lists = [[] for i in range(100000)]
phi = math.pi / 5
listnum = 0
for theta in np.arange(math.pi/2, 3*math.pi/2, 1*math.pi/40 ):
    for r in np.arange(0,15,0.01):
        if theta > math.pi:
            x= math.cos(theta-math.pi)
            Yc = legendre(L, M, x)  # 计算第一部分
        elif theta == math.pi:
            x= math.cos(1*math.pi/80) #做一个小量近似
            Yc = legendre(L, M, x)  # 计算第一部分
        else:
            x = math.cos(theta)
            Yc = legendre(L, M, x)  # 计算第一部分

        Ya = scinumsqrt(scinumdivide(scinummultiply(scinum(2 * L + 1), scinumfactorial(L - M)),
                                     scinummultiply(scinum(4 * math.pi), scinumfactorial(L + M))))  # 计算第二部分
        Yb_real = scinum(math.cos(M * phi))  # 计算第三部分的实部和虚部
        Yb_imag = scinum(math.sin(M * phi))

        SH_real = scinummultiply(Ya, scinummultiply(Yb_real, Yc))  # 计算球谐函数的实部和虚部
        SH_imag = scinummultiply(Ya, scinummultiply(Yb_imag, Yc))
        cof1 = math.exp(-r / 2) * r / math.sqrt(24)
        SH_realns = SH_real[0] * (10 ** SH_real[1])  # 非科学计数法
        SH_imagns = SH_imag[0] * (10 ** SH_imag[1])
        R = (cof1  * SH_realns) ** 2 + (cof1 * SH_imagns) ** 2
        lists[listnum].append(r)
        lists[listnum].append(theta)
        lists[listnum].append(R)
        listnum +=1

name = ['r', 'theta', 'Rou']
test=pd.DataFrame(columns=name,data=lists)
test.to_csv('M1data.csv',encoding='gbk')

#M = +1,phi

L = 1
M = 1
lists0 = [[] for t in range(100000)]
theta = math.pi / 2
listnum = 0
for phi in np.arange(math.pi/2, 21*math.pi/10, 1*math.pi/20):
    for r in np.arange(0,15,0.01):

        x = math.cos(theta)
        Ya = scinumsqrt(scinumdivide(scinummultiply(scinum(2 * L + 1), scinumfactorial(L - M)),
                                     scinummultiply(scinum(4 * math.pi), scinumfactorial(L + M))))  # 计算第二部分
        Yb_real = scinum(math.cos(M * phi))  # 计算第三部分的实部和虚部
        Yb_imag = scinum(math.sin(M * phi))
        Yc = legendre(L, M, x)  # 计算第一部分
        SH_real = scinummultiply(Ya, scinummultiply(Yb_real, Yc))  # 计算球谐函数的实部和虚部
        SH_imag = scinummultiply(Ya, scinummultiply(Yb_imag, Yc))
        cof1 = math.exp(-r / 2) * r / math.sqrt(24)
        SH_realns = SH_real[0] * (10 ** SH_real[1])  # 非科学计数法
        SH_imagns = SH_imag[0] * (10 ** SH_imag[1])
        R = (cof1 * SH_realns) ** 2 + (cof1 * SH_imagns) ** 2
        lists0[listnum].append(r)
        lists0[listnum].append(phi)
        lists0[listnum].append(R)
        listnum +=1

name = ['r', 'phi', 'Rou']
test=pd.DataFrame(columns=name,data=lists0)
test.to_csv('M1dataphi.csv',encoding='gbk')


#M = -1
L = 1
M = -1
lists1 = [[] for k in range(100000)]
phi = math.pi / 5
listnum = 0
for theta in np.arange(math.pi/2, 3*math.pi/2, 1*math.pi/40 ):
    for r in np.arange(0,15,0.01):
        if theta > math.pi:
            x= math.cos(theta-math.pi)
            Yc = legendre(L, M, x)  # 计算第一部分
        elif theta == math.pi:
            x = math.cos(1 * math.pi / 80)  # 做一个小量近似
            Yc = legendre(L, M, x)  # 计算第一部分
        else:
            x = math.cos(theta)
            Yc = legendre(L, M, x)  # 计算第一部分
        Ya = scinumsqrt(scinumdivide(scinummultiply(scinum(2 * L + 1), scinumfactorial(L - M)),
                                     scinummultiply(scinum(4 * math.pi), scinumfactorial(L + M))))  # 计算第二部分
        Yb_real = scinum(math.cos(M * phi))  # 计算第三部分的实部和虚部
        Yb_imag = scinum(math.sin(M * phi))

        SH_real = scinummultiply(Ya, scinummultiply(Yb_real, Yc))  # 计算球谐函数的实部和虚部
        SH_imag = scinummultiply(Ya, scinummultiply(Yb_imag, Yc))

        cof1 = math.exp(-r / 2) * r / math.sqrt(24)
        SH_realns = SH_real[0] * (10 ** SH_real[1])  # 非科学计数法
        SH_imagns = SH_imag[0] * (10 ** SH_imag[1])
        R = (cof1 * SH_realns) ** 2 + (cof1  * SH_imagns) ** 2
        lists1[listnum].append(r)
        lists1[listnum].append(theta)
        lists1[listnum].append(R)
        listnum +=1

name1 = ['r', 'theta', 'Rou']
test=pd.DataFrame(columns=name1,data=lists1)
test.to_csv('M-1data.csv',encoding='gbk')



#M = -1,phi

L = 1
M = -1
lists2 = [[] for p in range(100000)]
theta = math.pi / 2
listnum = 0
for phi in np.arange(1*math.pi/10, 21*math.pi/10, 1*math.pi/20 ):
    for r in np.arange(0,15,0.01):

        x = math.cos(theta)
        Ya = scinumsqrt(scinumdivide(scinummultiply(scinum(2 * L + 1), scinumfactorial(L - M)),
                                     scinummultiply(scinum(4 * math.pi), scinumfactorial(L + M))))  # 计算第二部分
        Yb_real = scinum(math.cos(M * phi))  # 计算第三部分的实部和虚部
        Yb_imag = scinum(math.sin(M * phi))
        Yc = legendre(L, M, x)  # 计算第一部分
        SH_real = scinummultiply(Ya, scinummultiply(Yb_real, Yc))  # 计算球谐函数的实部和虚部
        SH_imag = scinummultiply(Ya, scinummultiply(Yb_imag, Yc))
        cof1 = math.exp(-r / 2) * r / math.sqrt(24)
        SH_realns = SH_real[0] * (10 ** SH_real[1])  # 非科学计数法
        SH_imagns = SH_imag[0] * (10 ** SH_imag[1])
        R = (cof1 *SH_realns) ** 2 + (cof1 * SH_imagns) ** 2
        lists2[listnum].append(r)
        lists2[listnum].append(phi)
        lists2[listnum].append(R)
        listnum +=1

name = ['r', 'phi', 'Rou']
test=pd.DataFrame(columns=name,data=lists2)
test.to_csv('M-1dataphi.csv',encoding='gbk')


#M = 0

L = 1
M = 0
Rou = []

lists3 = [[] for s in range(100000)]
phi = math.pi / 5
listnum = 0
for theta in np.arange(1*math.pi/40, math.pi, 1*math.pi/40 ):
    for r in np.arange(0,15,0.01):

        x = math.cos(theta)
        Ya = scinumsqrt(scinumdivide(scinummultiply(scinum(2 * L + 1), scinumfactorial(L - M)),
                                     scinummultiply(scinum(4 * math.pi), scinumfactorial(L + M))))  # 计算第二部分
        Yb_real = scinum(1)  # 计算第三部分的实部和虚部
        Yc = legendre(L, 0, x)  # 计算第一部分
        SH_real = scinummultiply(Ya, scinummultiply(Yb_real, Yc))  # 计算球谐函数的实部和虚部
        cof1 = math.exp(-r / 2) * r / math.sqrt(24)
        SH_realns = SH_real[0] * (10 ** SH_real[1])  # 非科学计数法
        R = (cof1 * SH_realns) ** 2
        lists3[listnum].append(r)
        lists3[listnum].append(theta)
        lists3[listnum].append(R)
        listnum +=1

name2 = ['r', 'theta', 'Rou']
test=pd.DataFrame(columns=name2,data=lists3)
test.to_csv('M0data.csv',encoding='gbk')

#M = 0,phi

L = 1
M = 0
lists4 = [[] for b in range(100000)]
theta = math.pi / 2
listnum = 0
for phi in np.arange(1*math.pi/10, 21*math.pi/10, 1*math.pi/20 ):
    for r in np.arange(0,15,0.01):
        cof1 = math.exp(-r / 2) * r / math.sqrt(24)
        SH_realns = 0 # 非科学计数法
        R = (cof1 *  SH_realns) ** 2
        lists4[listnum].append(r)
        lists4[listnum].append(phi)
        lists4[listnum].append(R)
        listnum +=1

name = ['r', 'phi', 'Rou']
test=pd.DataFrame(columns=name,data=lists4)
test.to_csv('M0dataphi.csv',encoding='gbk')
