# Homework assign. 6
import numpy as np
import math

def compositeTrapezoid (f, a, b, n):
    h = (b - a)/n
    for j in range (0, n):
        ans = ((h/2) * (f(x) + f(x[j+1])))
    return ans

def compositeSimpson (f, a, b, n):
    if ( n%2 == 0):
        h = (b-a)/n
        XI0 = f(a) + f(b)
        XI1 = 0
        XI2 = 0
        for i in range (1, n-1):
            X = a + i*h
            if ( i%2 == 0):
                XI2 = XI2 + f(X)
            else:
                XI1 = XI1 + f(X)
            XI = h*(XI0 + (2*XI2) + (4*XI1)) / 3

            return XI

def romberg (f, a, b, p):
    R = []
    h = b - a
    R[1][1] = (h/2) * (f(a) + f(b))
    for i in range (2, n+1):
        R[2][1] = (1/2)*(R[1][1] + h*(f(a)+(k-0.5)*h))
        for j in range (2, i):
            R[2][j] = R[2][j-1] + (R[2][j-1] - R[1][j-1])/(4**(j-1)) - 1
        return R[2][j]
    h = h/2
        for j in range (1, i):
            R[1][j] = R[2][j]

def adaptiveSimp (f, a, b, TOL):
    APP = 0
    i = 1
    TOL = 10*TOL
    a[i] = a
    h[i] = (b-a)/2
    FA[i] = f(a)
    FC[i] = f(a+h[i])
    FB[i] = f(b)
    S[i] = (h[i]*(FA[i] + 4*FC[i] + FB[i])) /3

    L[i] = 1
    while (i > 0):
        FD = f(a[i] + (h[i]/2))
        FE = f(a[i] + (3*h[i]/2))
        S1 = ( h[i] * (FA[i] + 4*FD + FC[i]) )/ 6

        S2 = (h[i]*(FC[i] + 4*FE + FB[i])) / 6
        v1 = a[i]
        v2 = FA[i]
        v3 = FC[i]
        v4 = FB[i]
        v5 = h[i]
        v6 = TOL[i]
        v7 = S[i]
        v8 = L[i]

    i = i -1
    if abs(S1 + S2 - v7) < v6:
        APP = APP + (S1 + S2)
    else:
        if (v8 >= N):
            stop
        else:
            i = i + 1
            a[i] = v1 + v5
            FA[i] = v3
            FC[i] = FE
            FB[i] = v4
            h[i] = v5/2
            TOL[i] = v6/2
            S[i] = S2
            L[i] = v8 + 1

            i = i + 1
            a[i] = v1
            FA[i] = v2
            FC[i] = FD
            FB[i] = v3
            h[i] = h[i-1]
            TOL[i] = TOL[i-1]
            S[i] = S1
            L[i] = L[i - 1]
    return APP

def main ():




main()
