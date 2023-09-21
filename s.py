# Bisection

from sympy import var
from sympy import sympify

a = float(input("Enter the Value of lower limit of the interval : "))
b = float(input("Enter the Value of upper limit of the interval : "))

#---------------------------------------------------------------------
# Getting the Equation
#---------------------------------------------------------------------

x = var('x')
fx = input("Enter the Homogenous equation : ")
fx = sympify(fx)

#---------------------------------------------------------------------

while a >= b:

    print("\n\nWrong Interval")    

    a = float(input("Enter the Value of lower limit of the interval : "))
    b = float(input("Enter the Value of upper limit of the interval : "))


midpoint = float(a + b) / 2

while ( (str(a))[0:6] != (str(b))[0:6] ) :
    
    if ( (fx.subs(x,a)*fx.subs(x,midpoint)) < 0) :
        b = midpoint
    
    else :
        a = midpoint
    
    midpoint = (a + b) / 2
    print(a, " - ", b)    

print("\n\nThe root of the Equation is : ", a)

# Newton Raphson

import sympy as sp

def solve(x):
    return x**3 + (2*x)**2 + 10*x - 20


def newtonRaphson(a0,b0,der,x):
    x0=a0
    condition=True
   
    while condition==True:
        h=-(solve(x0))/der.subs(x,x0)
        x1=x0+h
        if solve(x0)==0 or abs(x1-x0)<0.001:
            return x0
        elif solve(a0)*solve(x1)<0:
            b0=x1
            x0=x1
        elif solve(b0)*solve(x1)<0:
             a0=x1
             x0=x1
       
   
x = sp.symbols('x')
eqn=x**3 + (2*x)**2 + 10*x - 20
der=sp.diff(eqn,x)

a0=float(input("Enter the starting interval : "))
b0=float(input("Enter the ending interval : "))

if solve(a0) * solve(b0) <0:
    final_val = newtonRaphson(a0,b0,der,x)
    print("The final result :",round(final_val,4))
else:
    print("\nVerification step conditions failed")

# Regula Falsi

from sympy import var
from sympy import sympify

a = float(input("Enter the Value of lower limit of the interval : "))
b = float(input("Enter the Value of upper limit of the interval : "))

#---------------------------------------------------------------------
# Getting the Equation
#---------------------------------------------------------------------

x = var('x')
fx = input("Enter the Homogenous equation : ")
fx = sympify(fx)

#---------------------------------------------------------------------

x1 = 1
xprev = 0

while ( (str(xprev))[0:6] != (str(x1))[0:6] ) :
    
    if(xprev == x1):
        break
    
    xprev = x1
    
    print(a, " - ", b)
    
    x1 = a + (( fx.subs(x,a)) * abs(b-a)) / ((fx.subs(x,a)) + (fx.subs(x,b)) )
    
    if( fx.subs(x,a) * fx.subs(x,b) < 0 ):
        b = x1
    
    else:
        a = x1

print("\n\nThe root of the Equation is : ", x1)

# Fixed Point Iteration

import sympy as sp

x = sp.symbols('x')
tolerance = 0.0001
equation_str = input("Enter the equation : ")
equation = sp.sympify(equation_str)
initial_guess = 0
max_iterations = 20
f = sp.lambdify(x, equation, 'numpy')

def fixed_point_iteration(f):
    x_n = initial_guess
    iteration = 0
    print("\nIteration table:")
    print("Iteration\t x_n")
    
    while iteration < max_iterations:
        x_n_plus_1 = f(x_n)
        print(f"{iteration}\t\t {x_n}")
        
        if abs(x_n_plus_1 - x_n) < tolerance:
            print("\nConverged to the desired tolerance.")
            return x_n_plus_1
        
        x_n = x_n_plus_1
        iteration += 1

    print("\nMaximum number of iterations reached.")
    return None

result = fixed_point_iteration(f)

if result is not None:
    print(f"\nApproximate solution: x = {result}")
    sp.plot(equation, (x, result - 2, result + 2), line_color='blue', title='Graph of the Equation', xlabel='x', ylabel='f(x)', show=True)
else:
    print("\nFixed-point iteration did not converge within the maximum number of iterations.")
 
# Linear Interpolation

import matplotlib.pyplot as plt


x1,y1=map(float,input("Enter (x1,y1):").split())
x2,y2=map(float,input("Enter (x2,y2):").split())
point=float(input("Enter x for prediction:"))
ya = y1 + (point - x1)*(y2 - y1)/(x2 - x1)

print(f"f({point}) = {ya}")

x=range(int(x1),int(x2)+1)
m=(y2-y1)/(x2-x1)
y=[]
c=y1-m*x1
for i in range(len(x)):
    y.append(int(m*x[i]+c))
    
print(f"y = {m}*x + {c}")
plt.plot(x,y)
plt.scatter([x1,x2],[y1,y2])
plt.scatter(point,ya,color='red')
 
# Quadratic Interpolation

import matplotlib.pyplot as plt
import numpy as np

def quadratic(x, b0, b1, b2, x1, x2):
    return b0 + b1*(x-x1) + b2*(x-x1)*(x-x2)

x1, y1 = map(float, input("Enter point 1 : ").split())
x2, y2 = map(float, input("Enter point 2 : ").split())
x3, y3 = map(float, input("Enter point 3 : ").split())
x = float(input("Enter x coordinate : "))

b0 = y1
b1 = (y2 - y1) / (x2 - x1)
b2 = (((y3 - y2) / (x3 - x2)) - b1) / (x3 - x1)

print("Quadratic polynomial : y = ", b0, " + ", b1, "(x - ", x1, ") + ", b2, "(x - ", x2, ")")

y = quadratic(x, b0, b1, b2, x1, x2)

print("f(", x, ") = ", y)

x_min = min(x1, x2, x3)
x_max = max(x1, x2, x3)
x_range = np.linspace(x_min, x_max , 100)
y_range = quadratic(x_range, b0, b1, b2, x1, x2)

plt.plot(x_range, y_range)
plt.scatter([x1, x2, x3], [y1, y2, y3])
plt.scatter(x, y, color="red")
plt.xlabel('x')
plt.ylabel('y')
plt.title('Quadratic Polynomial')
plt.show()

# Quadratic Interpolation with plot

import numpy as np
import matplotlib.pyplot as plt

def quadratic_interpolation(x, x0, x1, x2, y0, y1, y2):
    L0 = ((x - x1) * (x - x2)) / ((x0 - x1) * (x0 - x2))
    L1 = ((x - x0) * (x - x2)) / ((x1 - x0) * (x1 - x2))
    L2 = ((x - x0) * (x - x1)) / ((x2 - x0) * (x2 - x1))
    return y0 * L0 + y1 * L1 + y2 * L2

x_values = [0, 1, 2]
y_values = [0, 1, 20]

a0 = y_values[0]
a1 = (y_values[1] - y_values[0]) / (x_values[1] - x_values[0])
a2 = ((y_values[2] - y_values[1]) / (x_values[2] - x_values[1]) - (y_values[1] - y_values[0]) / (x_values[1] - x_values[0])) / (x_values[2] - x_values[0])

quadratic_polynomial = f"y = {a0} + {a1}(x - {x_values[0]}) + {a2}(x - {x_values[0]})(x - {x_values[1]})"
print(f"Polynomial: {quadratic_polynomial}")

x_range = np.linspace(min(x_values) - 10, max(x_values) + 10, 400)
y_range = a0 + a1 * (x_range - x_values[0]) + a2 * (x_range - x_values[0]) * (x_range - x_values[1])

plt.scatter(x_values, y_values, color='RED')
plt.plot(x_range, y_range, label='Polynomial')
plt.xlabel('x')
plt.ylabel('y')
plt.show()

# Divided difference

import math

def proterm(i, value, x):
    
    pro = 1
    for j in range(i):
        pro = pro * (value - x[j])
    return pro
 
def dividedDiffTable(x, y, n):
    
    for i in range(1, n):
        for j in range(n - i):
            y[j][i] = ((y[j][i - 1] - y[j + 1][i - 1])/(x[j] - x[i + j]))
    return y

def applyFormula(value, x, y, n):
    
    sum = y[0][0]
    for i in range(1, n):
        sum = sum + (proterm(i, value, x) * y[0][i])     
    return sum
 
def printDiffTable(y, n):
    
    for i in range(n):
        for j in range(n - i):
            print(round(y[i][j], 4), "\t", end = " ")
        print("")

i = 0
n = 0
x = [] 

n = int(input("\nEnter the number of Variables : "))

print("Enter the Values : ")
for i in range(0, n):
    
    dummy = float(input())
    x.append(dummy)

y = [[0 for i in range(n+1)] for j in range(n)]
    
print("Enter the f(x) values respetively : ")
for i in range(0, n):

    dummy = float(input())
    y[i][0] = dummy

#choice = input("Is the value to find for in exponential ? If yes Press Y/y else Press N/n ")
value = float(input("Enter the value to find for : "))

# to check for exp input
"""
if choice == 'Y' or choice =='y':
    value = float(input("Enter the power of the exponent : "))
    value = math.exp(value)
elif choice == 'N' or choice =='n':
    value = float(input("Enter the value to find for : "))
"""
y = dividedDiffTable(x, y, n)

print("\nDivided Difference Table : \n")
printDiffTable(y, n)
 
print("\nValue at", value, "is", round(applyFormula(value, x, y, n), 4))


# ---------------------------------------------------------------------------------------------------------

# Lagrange

import sympy as sp

# -----------------------------------------------------------------------------

a = []
y = []

# -----------------------------------------------------------------------------

n = int(input("Enter the order n : "))

# -----------------------------------------------------------------------------

for i in range(n):
    temp = float(input("Enter the values of x" + str(i) + " : " ))
    a.append(temp)
    temp = float(input("Enter the values of f(x" + str(i) + ") : "))
    y.append(temp)

# -----------------------------------------------------------------------------

yp = 0
pt = float(input("Enter the value of x to substitute : "))
x  = sp.Symbol('x')

for i in range(n):
    p = 1
    for j in range(n):
        if i != j:
            p = p * (x - a[j])/(a[i] - a[j])
    
    yp = yp + p * y[i] 

print("the interpolated expression is :", yp)
print("The interpolated value is :", round(yp.subs(x, pt), 4))

# ------------------------------------------------------------------------------------------------------------
-----------------------------------------

# Inverse Lagrange

a = []
y = []

n = int(input("Enter the order n : "))

# -----------------------------------------------------------------------------

for i in range(n):
    
    temp = float(input("Enter the values of x" + str(i) + " : " ))
    a.append(temp)
    temp = float(input("Enter the values of f(x" + str(i) + ") : "))
    y.append(temp)

# -----------------------------------------------------------------------------

yp = 0
pt = float(input("Enter the value of x to substitute:"))

# -----------------------------------------------------------------------------

for i in range(n):
    p = 1
    for j in range(n):
        if i != j:
            p = p * (pt - y[j]) / (y[i] - y[j])
    yp = yp + p * a[i] 
    
# -----------------------------------------------------------------------------

print("\nThe value of x is :", round(yp, 4))

# -----------------------------------------------------------------------------

# Newtons Forward Interpolation

import math, sys

# -----------------------------------------------------------------------------

x = []
n = int(input("Enter the order n : "))
y = [[0 for i in range(n)] for j in range(n)]

# -----------------------------------------------------------------------------

for i in range(n):
    x.append(float(input("Enter the values of x" + str(i) + " : " )))

for i in range(n):
    y[i][0] = float(input("Enter the values of y" + str(i) + " : "))

# -----------------------------------------------------------------------------

pt = float(input("Enter the point of interpolation : "))

# checkin for equally spaced
space = abs(x[1] - x[0])

for i in range(0, len(x) - 1):
    
    if(abs(x[i] - x[i+1]) != space):
        print("Not Equally Spaced ! Will not be able to calculate value !")
        sys.exit()

# -----------------------------------------------------------------------------
# Calculating the forward difference
# -----------------------------------------------------------------------------

for i in range(1, n):
    for j in range(n - i):
        y[j][i] = y[j+1][i-1] - y[j][i-1]

for i in range(n):
    for j in range(n-i):
        print(y[i][j], end = "\t")
    print()

# -----------------------------------------------------------------------------

res = y[0][0]
u = (pt - x[0]) / (x[1] - x[0])

# -----------------------------------------------------------------------------

for i in range(1,n):
    temp = u
    for j in range(1,i):
        temp *= (u - j)
    res = res + (temp * y[0][i]) / math.factorial(i)
    
# -----------------------------------------------------------------------------

print("The interpolated value at %0.0f is %0.2f" %(pt,res))

# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------------------------------

# Newtons Backward Interpolation

# import math, sys

# -----------------------------------------------------------------------------

x = []
n = int(input("Enter the order n :"))
y = [[0.0 for i in range(n)] for j in range(n)]

# -----------------------------------------------------------------------------

for i in range(n):
    x.append(float(input("Enter the values of x" + str(i) + " : " )))

for i in range(n):
    y[i][0] = float(input("Enter the values of y" + str(i) + " : "))

# -----------------------------------------------------------------------------

pt = float(input("Enter the point of interpolation : "))

# checkin for equally spaced
space = abs(x[1] - x[0])

for i in range(0, len(x) - 1):
    
    if(abs(x[i] - x[i+1]) != space):
        print("Not Equally Spaced ! Will not be able to calculate value !")
        sys.exit()

# -----------------------------------------------------------------------------
# Calculating the forward difference
# -----------------------------------------------------------------------------

for i in range(1,n):
    for j in range(n-1, i-1, -1):
        y[j][i] = y[j][i-1] - y[j-1][i-1]

for i in range(n):
    for j in range(i + 1):
        print(y[i][j], end = "\t")
    print()

# -----------------------------------------------------------------------------

res = y[n - 1][0]
u = (pt - x[n - 1]) / (x[1] - x[0])

# -----------------------------------------------------------------------------

for i in range(1, n):
    temp = u
    for j in range(i):
        temp = temp * (u + j)
    res = res + (temp * y[n-1][i]) / math.factorial(i)

# -----------------------------------------------------------------------------

print("The interpolated value at %0.0f is %0.2f" %(pt, res))

# ----------------------------------------------------------------------------------------------------

# Cubic Spline

import numpy as np
import sympy as sp

# -----------------------------------------------------------------------------

def ch(n):
  if n < 0:
    return '-'
  else:
    return '+'

# -----------------------------------------------------------------------------

n = int(input("Enter number of values : "))
x  = []
fx = []

for i in range(n):
  x.append(float(input("Enter x : ")))
  fx.append(float(input("Enter f(x) : ")))

'''n=4
x=[1,2,3,4]
fx=[1,5,11,8]'''

# -----------------------------------------------------------------------------

h = x[1] - x[0]
d = []

# -----------------------------------------------------------------------------

for i in range(1, n - 1):
  d.append(6 * (fx[i-1] - 2 * fx[i] + fx[i+1]) / (h*h) )
  
# -----------------------------------------------------------------------------

a = np.zeros((n-2, n-2))

for i in range(n - 2):
  for j in range(n - 2):
    
    if i == j:
      a[i][j] = 4
      
    elif abs(i-j) == 1:
      a[i][j] = 1

m = np.linalg.solve( a, np.array(d) )
m = m.tolist()
m.insert(0, 0)
m.insert(n, 0)

px = []

# -----------------------------------------------------------------------------

print("\nThe Cubic Equations : \n")
for i in range(n - 1):
  t = "%f*(%f - x)*3/6%d %c %f*(x %c %f)*3/6%d + (%f - x)(%f %c %f%d/6)/%d + (x %c %f)(%f %c %f%d/6)/%d"%(m[i],x[i+1],h,ch(m[i+1]),abs(m[i+1]),ch(-x[i]),abs(x[i]),h,x[i+1],fx[i],ch(-m[i]),abs(m[i]),h*h,h,ch(-x[i]),abs(x[i]),fx[i+1],ch(-m[i+1]),abs(m[i+1]),h*h,h)
  px.append(str(sp.expand(t)))
  print("p%d(x) = %s"%((i+1),px[i]))

X = float(input("\nEnter the x to find f(x) :"))

for i in range(len(x)):
  if x[i] >= X:
    break

print("\nf(%.2f) = %.2f"%(X,eval((px[i-1]).replace('x',str(X)))))

# ---------------------------------------------------------------------------------------------------