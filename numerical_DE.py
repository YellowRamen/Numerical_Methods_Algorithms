#coded in Python

import numpy as np
import math

def Euler (f, a, b, y0, h):
    N = int((b-a)/h)
    t = a
    w = y0
    for i in range (1, N+1):
        w = w+h*f(t, w)
        t = a + i*h
    return [w, t]

def RK2 (f, a, b, y0, h):
    N = int((b-a)/h)
    t = a
    w = y0
    for i in range (1, N+1):
        w = w+h*f(t+(h/2.0), w+(h/2.0)*f(t, w))
        t = a + i*h
    return [w, t]

def RK4 (f, a, b, y0, h):
    N = int((b-a)/h)
    t = a
    w = y0
    for i in range (1, N+1):
        k1 = h*f(t, w)
        k2 = h*f(t+(h/2.0), w+(k1/2.0))
        k3 = h*f(t+(h/2.0), w+(k2/2.0))
        k4 = h*f(t+h, w+k3)
        w = w+(k1+2.0*k2+2.0*k3+k4) / 6.0
        t = a + i*h
    return [w, t]

def AB4 (f, a, b, y0, h):
    N = int((b-a)/h)
    t = [a+i*h for i in range (0, N+1)]
    w = [0.0 for i in range (0, N+1)]
    w[0] = y0
    for i in range (1,4):
        k1 = h*f(t[i-1], w[i-1])
        k2 = h*f(t[i-1]+(h/2.0), w[i-1]+(k1/2.0))
        k3 = h*f(t[i-1]+(h/2.0), w[i-1]+(k2/2.0))
        k4 = h*f(t[i-1]+h, w[i-1]+k3)
        w[i] = w[i-1]+((k1+2.0*k2+2.0*k3+k4)/6.0)
        t[i] = a+i*h

    for i in range (4, N+1):
        w[i] = w[i-1] + (h*(55.0*f(t[i-1], w[i-1]) - 59.0*f(t[i-2], w[i-2]) + 37.0*f(t[i-3], w[i-3]) - 9.0*f(t[i-4], w[i-4]))) / 24.0
    return [w[N], t[N]]

def AMB4 (f, a, b, y0, h):
    N = int((b-a)/h)
    t = [a+i*h for i in range (0, N+1)]
    w = [0.0 for i in range (0, N+1)]
    w[0] = y0
    for i in range (1,4):
        k1 = h*f(t[i-1], w[i-1])
        k2 = h*f(t[i-1]+(h/2), w[i-1]+(k1/2))
        k3 = h*f(t[i-1]+(h/2), w[i-1]+(k2/2))
        k4 = h*f(t[i-1]+h, w[i-1]+k3)
        w[i] = w[i-1]+((k1+2*k2+2*k3+k4)/6)
        t[i] = a+i*h

    for i in range (4, N+1):
        wstar = w[i-1] +(h*(55.0*f(t[i-1], w[i-1]) - 59.0*f(t[i-2], w[i-2]) + 37.0*f(t[i-3], w[i-3]) - 9.0*f(t[i-4], w[i-4]))) / 24.0
        w[i] = w[i-1] + (h*(9.0*f(t[i], wstar) + 19.0*f(t[i-1], w[i-1]) - 5.0*f(t[i-2],w[i-2]) + f(t[i-3],w[i-3]))) / 24.0
    return [w[N], t[N]]

def main ():
    x = lambda t, y: t/y
    y = lambda t: np.sqrt(1+t**2)
    
    y0 = 1.0
    a = 0.0
    b = 10.0
    function = (Euler, RK2, RK4, AB4, AMB4)
    for f in function:
        beforeError = 0
        error = 0
        n = 50
        print('This is ' + str(f.__name__))
        print('h\t\terror\t\t\tratio')
        while (n <= 1600):
            h = 10.0/n
            result = f(x, a, b, y0, h)[0]
            print(str(h)+str('\t\t')),
            error = abs(result - y(10))
            print (error),
            if (n != 50):
                print (str('\t') + str(error/ beforeError)),
            print('')
            beforeError = error
            n = n*2
main()
