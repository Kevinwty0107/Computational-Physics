# coding=UTF-8
# 此程序用来将十进制数转化为双精度浮点数
for k in range(1,100):
    num = input("Enter the input decimal number: ")
    num_decimal = float(num)
    if num_decimal == 0:  # 如果数为0，则直接输出
        print(64 * '0')
        exit()
    while abs(num_decimal) > 2.23E308 or abs(num_decimal) < 1.79E-308:  # 检验是否在可表示的范围内
        print("您输入的数不在双精度浮点数表示范围内，请重新输入：\n")
        num_decimal = input("Enter your input: ")
    if num_decimal < 0:  # 看正负
        pm = str(1)
    else:
        pm = str(0)
    num_decimal = abs(num_decimal)  # 取绝对值
    bit = 52
    integer = int(num_decimal)
    decimal = num_decimal - integer
    #print(bin(integer))
    if integer > 1:  # 整数部分大于1
        num_binary_1 = str(bin(integer)[3:])  # 整数部分
        bit = bit - len(num_binary_1)
        num_binary_2 = ""  # 小数部分
        i = 0
        while i < bit:
            result = int(decimal * 2)
            decimal = decimal * 2 - result
            num_binary_2 += str(result)
            i += 1
        index = str(bin(1023 + len(num_binary_1))[2:])  # 指数部分
        x = len(index)
        num_binary = pm + index + num_binary_1 + num_binary_2
    elif integer == 1:  # 整数部分等于1
        num_binary_2 = ""  # 小数部分
        i = 0
        while i < bit:
            result = int(decimal * 2)
            decimal = decimal * 2 - result
            num_binary_2 += str(result)
            i += 1
        index = "01111111111"  # 指数部分
        x = len(index)
        num_binary = pm + index + num_binary_2
    else:  # 整数部分等于0
        num_binary_2 = ""  # 小数部分
        result = 0
        t = 0
        while result == 0:
            result = int(decimal * 2)
            decimal = decimal * 2 - result
            t += 1
        i = 0
        while i < bit:
            result = int(decimal * 2)
            decimal = decimal * 2 - result
            num_binary_2 += str(result)
            i += 1
        index = str(bin(1023 - t))[2:]  # 指数部分
        x = len(index)
        index_str = (11 - x) * '0' + index
        num_binary = pm + index_str + num_binary_2
    print(num_binary)

    cob = input("Continue or Break: ")
    if cob == "C" or cob == "Continue" :
        continue
    elif cob == "B" or cob == "Break":
        break
    else:
        print ("Invalid command! ")
        break








