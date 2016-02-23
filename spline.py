
import math
import numpy as np
import matplotlib.pyplot as plt

def natural_spline (x, y):
    n = len(x)
    h = [0.0]*n

    for i in range (0 , n-1):
        h[i] = x[i+1] - x[i]

    alpha = [0.0]*n
    for i in range (1, n-1):
        alpha[i] = ((3/h[i])*(y[i+1] - y[i])) - ((3/h[i-1])*(y[i] - y[i-1]))

    l = [0.0]*n
    l[0] = 1.0
    u = [0.0]*n
    z = [0.0]*n
    for i in range (1, n-1):
        l[i] = 2*(x[i+1] - x[i-1]) - (h[i-1]*u[i-1])
        u[i] = h[i] / l[i]
        z[i] = (alpha[i] - h[i-1]*z[i-1]) / l[i] 

    l[n-1] = 1
    z[n-1] = 0

    b = [0.0]*n
    c = [0.0]*n
    d = [0.0]*n
    for j in range (n-2, -1, -1):
        c[j] = z[j] - (u[j]*c[j+1])
        b[j] = ((y[j+1] - y[j])/h[j]) - (h[j]*(c[j+1]+2*c[j])/3)
        d[j] = (c[j+1] - c[j])/(3*h[j])

    return [y,b,c,d]

def notAknot_spline (x,y):
    n = len(x)
    h = [0.0]*n

    for i in range (0 , n-1):
        h[i] = x[i+1] - x[i]

    alpha = [0.0]*n
    for i in range (1, n-1):
        alpha[i] = ((3/h[i])*(y[i+1] - y[i])) - ((3/h[i-1])*(y[i] - y[i-1]))

    l = [0.0]*n
    l[0] = 1.0
    u = [0.0]*n
    z = [0.0]*n
    for i in range (1, n-1):
        l[i] = 2*(x[i+1] - x[i-1]) - (h[i-1]*u[i-1])
        u[i] = h[i] / l[i]
        z[i] = (alpha[i] - h[i-1]*z[i-1]) / l[i] 

    l[n-1] = 1
    z[n-1] = 0

    b = [0.0]*n
    c = [0.0]*n
    d = [0.0]*n
    for j in range (n-2, -1, -1):
        c[j] = z[j] - (u[j]*c[j+1])
        b[j] = ((y[j+1] - y[j])/h[j]) - (h[j]*(c[j+1]+2*c[j])/3)
        d[j] = (c[j+1] - c[j])/(3*h[j])

    return [y,b,c,d]

def eval_spline(x, C, z):
    results = [0.0]*len(z)
    index = 0
    for point in z:
        for i in range(0, len(x)+1):
            if (point >= x[i] and point <= x[i+1]):
                S = C[0][i] + (C[1][i] * (point-x[i])) + (C[2][i] * (point-x[i])**2) + (C[3][i] *(point-x[i])**3)
                results[index] = S
                index = index+1
                break
    return results

def error(hlist, fun, splinefun, numPoints):
    for h in hlist:
        print(h)
        cut = 1.0/h
        x_n = np.linspace(0,1,cut)
        space = np.linspace(0,1,numPoints)
        C = splinefun(x_n, [fun(num) for num in x_n])
        result = eval_spline(x_n, C, space)
        print 'Max error is: ' + str(max([abs(fun(space[i]) - result[i]) for i in xrange(0, len(space))]))
        print("")
def main():
    L = np.linspace(0,1,10)
    
#question 2
    f = lambda x: math.cos(3.0*math.pi *x)
    h = 0.1
    C = natural_spline(L, [f(num) for num in L])
    W = eval_spline(L , C, np.linspace(0,1,400))
    plt.subplot(2,1,1)
    plt.plot (np.linspace(0,1,400), W, '-')
    plt.xlabel('x')
    plt.ylabel('MaxError')
    plt.title('Question 2a')
    plt.show()
    error([.1, 0.05, 0.025, 0.0125,0.00625], f, natural_spline, 200)
    
main()
