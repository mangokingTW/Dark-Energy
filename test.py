from sympy import *
from constants import c
from utils import *

# Define n as a symbolic variable
n, x = symbols('n, x')

# Calculate the sum
result = solve([Eq(x, Product(exp(log(1/n)), (n, 1, 10000)).doit())], (x), dict=True)[0]

# Print the result (remains symbolic)
printE(result[x])


n, x = symbols('n, x')
result2 = solve([Eq(x, Sum(1/n, (n, 1, 10000)).doit())], (x), dict=True)[0]
printE(result2[x])