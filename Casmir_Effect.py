from constants import *
from sympy import abc
from sympy import *


def calculate_Casimir_Effect(gap):
    P10nm, Pa10nm = symbols(
        'P10nm, Pa10nm',
        real=True
    )

    eq = [
        Eq(P10nm,
           Sum(1 / (abc.k * Lp) ** 4, (abc.k, 1, gap / Lp)).doit() * (3 / 4 / pi) * c * h / 2 * Lp ** 4 / (
                   2 * pi ** 2 * gap ** 4)),
        Eq(Pa10nm, (- P10nm) * eV2J)
    ]
    print(solve(eq, (P10nm, Pa10nm), dict=True)[0][Pa10nm])
    print(- h * c * float(pi) / 480 / gap ** 4 * eV2J)


calculate_Casimir_Effect(mpf(1e-8))
calculate_Casimir_Effect(c)


def calculate_ZPE(d):
    ZPE = symbols('ZPE', real=True)
    eq = [
        Eq(ZPE, Sum(1 / abc.k, (abc.k, 1, d / Lp)).doit() / Lp * c * h / 2 / (4 / 3 * pi * d ** 3) * eV2J)
    ]
    print(solve(eq, ZPE, dict=True)[0])


calculate_ZPE(c)
calculate_ZPE(1e-8)


def calculate_T(d):
    T = symbols('T', real=True)
    eq = [
        Eq(T,
           Sum(1 / abc.k, (abc.k, 1, d / Lp)).doit() / Lp)
    ]
    print(solve(eq, T, dict=True)[0])


calculate_T(1e-8)


def calculate_P_vac(luminosity, radius, redshift):
    print(locals())
    print(f"Redshift/Radius: {redshift / radius}")
    P_vac, P_r = symbols('P_vac, P_r', real=True)
    abc.k = 1 / 2000  # Probability of photon interacting with ZPE
    flux = luminosity / (4 / 3 * pi * radius ** 3)
    eq = [
        Eq((luminosity - P_r * luminosity) / luminosity, 1 / (1 + redshift)),
        Eq((P_r * c / (radius * abc.k)), P_vac)
    ]
    res = solve(eq, (P_r, P_vac))
    print(f"P_vac: {res[P_vac]}")
    return res[P_vac]


# 通過測量卡西米爾效應的強度，物理學家可以估計真空中的能量密度。根據目前的測量結果，真空中的能量密度約為 10^-9 J/m^3。
# SN 1987A
# calculate_P_vac(1e36, float(1.68e5 * ly2m), 0.001067)
# SN 2014J
calculate_P_vac(1, float(3.5e6 * pc2m), 0.000677)
# SN 2008D
calculate_P_vac(1, float(27e6 * pc2m), 0.007)
# SN 2006gy
calculate_P_vac(1, float(73e6 * pc2m), 0.0192)
