#!/usr/bin/python3

# f(x)=x^(1.02)
# output f(10)


def one_x(x):
    return 1


def jiecheng(n):
    if n == 0:
        return 1
    elif n == 1:
        return 1
    y = 1
    for i in range(1, n+1):
        y *= i
    return y


def daoshu(n):
    xishu = 1
    zhishu = 1.02
    if n == 0:
        return 1
    else:
        for _ in range(n-1):
            zhishu = zhishu - 1
            xishu = xishu * zhishu
        return xishu * one_x(zhishu)


def taile(n):
    if n == 0:
        return 1
    return daoshu(n) / jiecheng(n) * (10-1)**n


def main():
    he = 0
    for i in range(100):
        he += taile(i)
        print('No:{} taile:{} he:{}'.format(
            i, taile(i), he))


if __name__ == '__main__':
    main()
