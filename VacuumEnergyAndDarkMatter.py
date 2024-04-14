from sympy import *
from constants import *
from utils import *

# printE(sqrt(hbar*c**5/4/G)*J2kg*G*2/c**2)
# printE(2.443e18*J2eV/(4/3*pi*4e-26**3))
# printE(2.443e18*J2eV/(4/3*pi*(138e9*ly2m**3)))

L = 8.8e26
printE("Max wave length (m)", L)
N_Planck_Mass = (Mp * kg2J * J2eV) / h / c * L * 2
printE("N of Nth energy equal to Planck mass (unitless)", N_Planck_Mass)
x = symbols('x', positive=True)
Ratio_Planck_Mass = zeta(4) - Sum(1 / x ** 4, (x, 1, 1.733097e+30)).doit()
# Alternative way to calculate value, because Sum function cannot handle upper bound larger than 1e35
Ratio_Planck_Mass = Ratio_Planck_Mass * 1e-93 / zeta(4)
printE("Ratio of energy above Plank mass (unitless)", Ratio_Planck_Mass)
Planck_Mass_J = Mp * kg2J
printE("Planck Mass (J)", Planck_Mass_J)
# E_n = n h c / wave_length / 2
# k_n = n pi / wave_length
# p_n = hbar * k_n
# x_n * p_n = hbar / 2
# x_n = h c / (pi E_n)
Length = hJ * c / (pi * Planck_Mass_J)
printE("Length (m)", Length)
VED = Planck_Mass_J * N_Planck_Mass * Ratio_Planck_Mass / ((Length) ** 3) * J2kg / 1.6
printE("Vacuum Energy Density (J/m^3)",
       VED
       )

A = (Length) ** 3 / Ratio_Planck_Mass
B = (4 / 3 * pi) * (L / 2) ** 3
printE(A)
printE(B)
printE(A/B)
