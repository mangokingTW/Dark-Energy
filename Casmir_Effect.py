from constants import *
from sympy import abc
from sympy import *

VU, LU, P_vac, P10nm, Pa10nm = symbols(
    'VU LU, P_vac, P10nm, Pa10nm',
    domain=S.Reals
)
r = mpf(1e-8)
X20 = 1e0
eq = [
    Eq(LU, mpf(8.8e26)),
    Eq(VU, 1 / 6 * pi * LU ** 3),
    Eq(P_vac, Sum(1 / (abc.k * Lp) ** 4, (abc.k, 1, LU / Lp)).doit() * (3 / 4 / pi) * c * h / 2 * Lp ** 2.5 / (
            2 * pi ** 2 * LU ** 2.5)),
    Eq(P10nm,
       Sum(1 / (abc.k * Lp) ** 4, (abc.k, 1, r / Lp)).doit() * (3 / 4 / pi) * c * h / 2 * Lp ** 4 / (
                   2 * pi ** 2 * r ** 4)),
    Eq(Pa10nm, (P_vac - P10nm) * eV2J)
]
print(solve(eq, (VU, LU, P_vac, P10nm, Pa10nm)))

print(h * c * float(pi) / 480 / r ** 4 * eV2J)

if True:
    # SN 1987A
    R = mpf(1.586038e21)
    Watt = mpf(1e36)
    SN = Watt / (4 / 3 * pi * R ** 3)

    eq2 = [
        Eq((SN - P_vac * eV2J / R ** 3) / SN, 1 / 1.0031)
    ]

    print(solve(eq2, P_vac))

if True:
    # SN 2003fg
    R = mpf(4e9*ly2m)
    Watt = mpf(1e36)
    SN = Watt / (4 / 3 * pi * R ** 3)

    eq2 = [
        Eq((SN - P_vac * eV2J / R ** 3) / SN, 1 / 1.0244)
    ]

    print(solve(eq2, P_vac))

if True:
    # GN-z11
    R = mpf(9.8e9*pc2m)
    Watt = mpf(3*3.828e38)
    SN = Watt / (4 / 3 * pi * R ** 3)

    eq2 = [
        Eq((SN - P_vac * eV2J / R ** 2) / SN, 1 / 1.024934)
    ]

    print(solve(eq2, P_vac))
