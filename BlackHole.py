from utils import *
from DarkEnergy import calculate_Casimir_Effect, calculate_ZPE
from constants import *

from sympy import *


# Schwarzschild radius (SWC_R)
def calculate_SWC_Radius(energy):
    matter = energy / c ** 2
    SWC_R = symbols('SWC', real=True, positive=True)
    eq = [
        Eq(SWC_R, 2 * G * matter)
    ]
    res = solve(eq, SWC_R, dict=True)[0]
    return res[SWC_R]


# Schwarzschild density (SWC_D)
def calculate_SWC_Density(energy):
    matter = energy / c ** 2
    SWC_D = symbols('SWC_D', real=True, positive=True)
    eq = [
        Eq(SWC_D, c ** 6 / (6 * pi * G ** 3 * matter ** 2))
    ]
    res = solve(eq, SWC_D, dict=True)[0]
    return res[SWC_D]


# Wave length won't create black hole
def calculate_SWC_Length():
    SWC_L = symbols('SWC_L', real=True, positive=True)
    eq = [
        Eq(SWC_L ** 6, ((9 * h ** 3 * G ** 3) / (2 * pi ** 3 * c ** 7)))
    ]
    res = solve(eq, SWC_L, dict=True)[0]
    return res[SWC_L]


# Black hole evaporation time
def calculate_BHET(matter):
    printE(matter)
    time = symbols('time', real=True)
    eq = [
        Eq(time, 5120 * pi * G ** 2 * matter ** 3 / (hbar * c ** 4))
    ]
    res = solve(eq, time, dict=True)[0]
    return res[time]


def doit():
    energy = mpf(h * c / lp)
    printE("SWC_R", calculate_SWC_Radius(energy))
    printE("SWC_D", calculate_SWC_Density(energy))
    SWC_L = calculate_SWC_Length()
    printE("SWC_L", SWC_L)
    ZPE = calculate_ZPE(c)
    printE("ZPE", ZPE)
    printE("Casimir Effect", calculate_Casimir_Effect(SWC_L))
    WL = 3.280000e-27
    printE(calculate_BHET((h * c / WL) / c ** 2))

    printE("SWC_D", calculate_SWC_Density(mpf(h * c / SWC_L)))
    # m = h * c / 3e-7
    # printE("ZPE", 2 * m ** 4 * c ** 5 / hbar ** 3 * eV2J)
    printE("Casimir Effect", calculate_Casimir_Effect(3e-7))


if __name__ == '__main__':
    doit()
