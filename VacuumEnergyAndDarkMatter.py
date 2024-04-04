from sympy import *
from constants import *
from utils import *

# printE(sqrt(hbar*c**5/4/G)*J2kg*G*2/c**2)
# printE(2.443e18*J2eV/(4/3*pi*4e-26**3))
# printE(2.443e18*J2eV/(4/3*pi*(138e9*ly2m**3)))

L = 8.8e26
printE(L)
alpha = (Mp*kg2J*J2eV)/h/c*L
printE("alpha", alpha)
x = symbols('x', positive=True)
gap = zeta(4)-Sum(1/x**4, (x, 1, 8.665485e+30)).doit()
gap = gap * 1e-90
printE("gap", gap)
eee = Mp*kg2J
printE("Planck Mass in Joule", eee)
printE("Vacuum Energy Density", (3*eee/2)*alpha*gap/(4/3*pi*lp**3))
