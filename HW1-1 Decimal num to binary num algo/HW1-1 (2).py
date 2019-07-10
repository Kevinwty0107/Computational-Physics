# coding=UTF-8
# 此程序用来将二进制浮点数转化为十进制数

for k in range(1,100):
    indexnumber = 0
    decimalnumber = 0
    num = input("Enter the input binary number:")  # 输入二进制数
    while len(num) != 64:  # 检验是否是64位
        print("please enter again:\n")
        num = input("Enter the input binary number:")
    if num == 64 * '0':
        print(0)
        exit()
    index = int(num[1:12], 2) - 1023
    indexnumber = 2 ** index  # 指数部分
    for x in range(1, 53):
        decimalnumber += int(num[11 + x]) * (2 ** (-x))  # 小数部分
    number = indexnumber * (1 + decimalnumber)
    if num[0] == '1':
        number = -number
    print(number)
    cob = input("Continue or Break: ")
    if cob == "C" or cob == "Continue":
        continue
    elif cob == "B" or cob == "Break":
        break
    else:
        print("Invalid command! ")
        break
