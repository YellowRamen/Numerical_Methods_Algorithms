import numpy as np
import math

def LU(A):
    n =  len(A)
    rowSwitch = [ i for i in xrange(0, n) ]
    # k represents the current pivot row. Since GE traverses the matrix in the upper 
    # right triangle, we also use k for indicating the k-th diagonal column index.
    for k in xrange(n-1):
        #Choose largest pivot element below (and including) k
        maxindex = abs(A[k:,k]).argmax() + k
        if A[maxindex, k] == 0:
            raise ValueError("Matrix is singular.")
        #Swap rows
        if maxindex != k:
            A[[k,maxindex]] = A[[maxindex, k]]
            temp = rowSwitch[k]
            rowSwitch[k] = rowSwitch[maxindex]
            rowSwitch[maxindex] = temp
        for row in xrange(k+1, n):
            multiplier = A[row][k]/A[k][k]
            #the only one in this column since the rest are zero
            A[row][k] = multiplier
            for col in xrange(k + 1, n):
                A[row][col] = A[row][col] - multiplier*A[k][col]
    L = [ [0.0]* n for i in xrange (0, n) ]
    U = [ [0.0]*n for i in xrange (0, n) ]
    P = [ [0.0]*n for i in xrange(0, n) ]
    for i in xrange(0, n):
        L[i][i] = 1.0
    for i in xrange(1, n):
        for j in range(0, i):
            L[i][j] = A[i][j]
    for i in xrange(0, n):
        for j in xrange(i, n):
            U[i][j] = A[i][j]
    for i in xrange(0, n):
        P[i][rowSwitch[i]] = 1.0
    return [L,U,P]

def GEPP(A, b):
    n =  len(A)
    if b.size != n:
        raise ValueError("Invalid argument: incompatible sizes between A & b.", b.size, n)
    # k represents the current pivot row. Since GE traverses the matrix in the upper 
    # right triangle, we also use k for indicating the k-th diagonal column index.
    for k in xrange(n-1):
        #Choose largest pivot element below (and including) k
        maxindex = abs(A[k:,k]).argmax() + k
        if A[maxindex, k] == 0:
            raise ValueError("Matrix is singular.")
        #Swap rows
        if maxindex != k:
            A[[k,maxindex]] = A[[maxindex, k]]
            b[[k,maxindex]] = b[[maxindex, k]]
        for row in xrange(k+1, n):
            multiplier = A[row][k]/A[k][k]
            #the only one in this column since the rest are zero
            A[row][k] = multiplier
            for col in xrange(k + 1, n):
                A[row][col] = A[row][col] - multiplier*A[k][col]
            #Equation solution column
            b[row] = b[row] - multiplier*b[k]
    x = np.zeros(n)
    k = n-1
    x[k] = b[k]/A[k,k]
    while k >= 0:
        x[k] = (b[k] - np.dot(A[k,k+1:],x[k+1:]))/A[k,k]
        k = k-1
    return x

def main():
    # Question # 3
    A = np.array( [ [4.0, 1.0, 1.0, 0.0, 0.0, 0.0], [1.0, 4.0, 1.0, 1.0 ,0.0, 0.0], [1.0, 1.0, 4.0, 1.0, 1.0, 0.0], [0.0, 1.0, 1.0, 4.0, 1.0, 1.0], [0.0, 0.0, 1.0, 1.0, 4.0, 1.0], [0.0, 0.0, 0.0, 1.0, 1.0, 4.0] ] )
    b = np.array( [ [1.0], [2.0], [3.0], [4.0], [5.0], [6.0] ] )
    C = np.array( [ [1.0, 1.0, 1.0, 1.0], [1.0, 1.0, 2.0, 3.0], [-1.0, 0.0, 2.0, 1.0], [3.0, 2.0, -1.0, 0.0] ] )
    d = np.array( [ [1.0], [2.0], [1.0], [1.0] ] )
    E = np.array( [ [5.0, 5.0, 4.0, 2.0], [3.0, 3.0, 4.0, -2.0], [-2.0, 0.0, 2.0, 0.6], [-1.0, 2.0, 3.4, -1.0] ] )
    f = np.array( [ [4.0], [3.0], [2.0], [1.0] ] )
    print ("Question 3 Part a is: ")
    print (GEPP(A, b))
    print ("Question 3 Part b is: ")
    print (GEPP(C, d))
    print ("Question 3 Part c is: ")
    print (GEPP(E, f))
    print("")

    # Question # 4
    F = np.array( [ [13.0, 39.0, 2.0, 57.0, 28.0], [-4.0, -12.0, 0.0, -19.0, -9.0,], [3.0, 0.0, -9.0, 2.0, 1.0], [6.0, 17.0, 9.0, 5.0, 7.0], [19.0, 42.0, 17.0, 107.0, 44.0] ] )
    g = np.array( [ [1.0], [1.0], [1.0], [1.0], [1.0] ] )
    result = (LU(F))
    print("Solution to Question 4 part a: ")
    print ("")
    L = np.array(result[0])
    print ("This is L: ")
    print (L)
    print ("")
    U = np.array(result[1])
    print("This is U: ")
    print (U)
    print ("")
    P = np.array(result[2])
    print ("This is P: ")
    print (P)
    print ("")
    Pg = np.dot(P, g)
    y = GEPP(L, Pg)
    x = GEPP(U, y)
    print("Solution to Question 4 part b: ")
    print (x)
main()
'''
Question 3 Part a is: 
[ 0.07937164  0.30219099  0.38032245  0.33154196  0.76560562  1.2257131 ]
Question 3 Part b is: 
[-2.  3. -1.  1.]
Question 3 Part c is: 
[ 0.79558824 -1.53088235  1.61029412  0.61764706]

Solution to Question 4 part a: 

This is L: 
[[ 1.          0.          0.          0.          0.        ]
 [ 0.68421053  1.          0.          0.          0.        ]
 [ 0.15789474 -0.64615385  1.          0.          0.        ]
 [ 0.31578947  0.36410256 -0.39862543  1.          0.        ]
 [-0.21052632 -0.30769231 -0.03436426  0.07070707  1.        ]]

This is U: 
[[  1.90000000e+01   4.20000000e+01   1.70000000e+01   1.07000000e+02
    4.40000000e+01]
 [  0.00000000e+00   1.02631579e+01  -9.63157895e+00  -1.62105263e+01
   -2.10526316e+00]
 [  0.00000000e+00   0.00000000e+00  -1.79076923e+01  -2.53692308e+01
   -7.30769231e+00]
 [  0.00000000e+00   0.00000000e+00   0.00000000e+00  -3.30000000e+01
   -9.04123711e+00]
 [  0.00000000e+00   0.00000000e+00   0.00000000e+00   0.00000000e+00
    3.54056024e-03]]

This is P: 
[[ 0.  0.  0.  0.  1.]
 [ 1.  0.  0.  0.  0.]
 [ 0.  0.  1.  0.  0.]
 [ 0.  0.  0.  1.  0.]
 [ 0.  1.  0.  0.  0.]]

Solution to Question 4 part b: 
[ -75.43627451  -88.66176471   -7.19117647  -98.60784314  359.80392157]
'''
