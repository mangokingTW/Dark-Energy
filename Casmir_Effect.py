from constants import *
from sympy import abc
from sympy import *
import pprint


def calculate_Casimir_Effect(gap):
    # P10nm, Pa10nm = symbols(
    #     'P10nm, Pa10nm',
    #     real=True
    # )
    #
    # eq = [
    #     Eq(P10nm,
    #        Sum(1 / (abc.k * Lp) ** 4, (abc.k, 1, gap / Lp)).doit() * (3 / 4 / pi) * c * h / 2 * Lp ** 4 / (
    #                2 * pi ** 2 * gap ** 4)),
    #     Eq(Pa10nm, (- P10nm) * eV2J)
    # ]
    # print(solve(eq, (P10nm, Pa10nm), dict=True)[0][Pa10nm])
    print(- h * c * float(pi) / 480 / gap ** 4 * eV2J)


def calculate_ZPE(distance):
    ZPE, k = symbols('ZPE, k', real=True, positive=True)
    eq = [
        Eq(ZPE, Sum(1 / k, (k, 1, distance / Lp)).doit() / Lp * c * h / 2 / (4 / 3 * pi * distance ** 3) * eV2J)
    ]
    a = solve(eq, ZPE, dict=True)[0]
    print(a)
    return a[ZPE]


def calculate_Redshift_Energy(d):
    print(locals())
    luminosity, distance, redshift = d['luminosity'], d['distance'], d['redshift']
    E_rs, ratio = symbols('E_rs, ratio', real=True)
    flux = luminosity / (4 / 3 * pi * distance ** 3)
    eq = [
        Eq((luminosity - ratio * luminosity) / luminosity, 1 / (1 + redshift)),
        Eq((ratio * luminosity), E_rs)
    ]
    res = solve(eq, (ratio, E_rs), dict=True)[0]
    print(res)
    d['energy_redshift'] = res[E_rs]


def calculate_Energy_Density(d):
    energy, radius = d['energy_redshift'], d['distance']
    energy_density, volume, area = symbols('energy_density, volume, area', real=True)
    eq = [
        Eq(area, 4 * pi * radius ** 2),
        Eq(volume, 4 / 3 * pi * radius ** 3),
        Eq(energy_density, energy / volume)
    ]
    res = solve(eq, energy_density, volume, area, dict=True)
    print(res)
    d['energy_density'] = res[0][energy_density]


def calculate_Ers_Density(d):
    energy, radius = d['energy_density'], d['distance']
    ers_density, volume = symbols('ers_density, volume', real=True)
    eq = [
        Eq(volume, pi * radius ** 2),
        Eq(ers_density, energy * volume)
    ]
    res = solve(eq, ers_density, volume, dict=True)
    print(res)
    d['ers_density'] = res[0][ers_density]


def distance_from_redshift(redshift):
    distance = symbols('distance', real=True)
    eq = [
        Eq(distance, redshift * c / H0)
    ]
    res = solve(eq, distance, dict=True)
    return res[0][distance]


# 通過測量卡西米爾效應的強度，物理學家可以估計真空中的能量密度。根據目前的測量結果，真空中的能量密度約為 10^-9 J/m^3。

# calculate_Casimir_Effect(mpf(1e-8))
# calculate_Casimir_Effect(c)
#
# calculate_ZPE(c)
# calculate_ZPE(8.8e26)

def doit():
    # GN-z11
    # calculate_Redshift_Energy(1e10, float(4.1e9 * pc2m), 11.09)
    # SN 2014J
    # SN_2014J = {
    #     'luminosity': 1e44,
    #     'distance': float(26.5e6 * pc2m),
    #     'redshift': 0.0008
    # }
    # calculate_Redshift_Energy(1e44, float(3.5e6 * pc2m), 0.0008)
    # SN 2008D
    SN_2008D = {
        'luminosity': 1e44,
        'distance': float(27e6 * pc2m),
        'redshift': 0.007
    }
    # calculate_Redshift_Energy(1e44, float(27e6 * pc2m), 0.007)
    # SN 2006gy
    SN_2006gy = {
        'luminosity': 1e44,
        'distance': float(73e6 * pc2m),
        'redshift': 0.0192
    }

    # calculate_Redshift_Energy(SN_2014J)
    calculate_Redshift_Energy(SN_2008D)
    calculate_Redshift_Energy(SN_2006gy)

    # print(SN_2014J['energy_redshift'] / SN_2014J['distance'])
    print(SN_2008D['energy_redshift'] / SN_2008D['distance'])
    print(SN_2006gy['energy_redshift'] / SN_2006gy['distance'])
    print(distance_from_redshift(SN_2006gy['redshift']))

    calculate_Energy_Density(SN_2006gy)
    calculate_Energy_Density(SN_2008D)

    calculate_Ers_Density(SN_2006gy)
    calculate_Ers_Density(SN_2008D)
    pprint.pprint(locals(), indent=2)
    RU = c
    calculate_ZPE(RU)

    def printE(f):
        print(format(float(f), "e"))

    # J.s-1.m-3 功率密度
    printE(SN_2006gy['ers_density'] / (4 / 3 * pi * RU ** 3))
    printE(9.47e-27 * kg2J)
    printE(distance_from_redshift(1089))
    printE(SN_2006gy['energy_redshift']/SN_2006gy['distance'] / (4 / 3 * pi * RU ** 3))


if __name__ == '__main__':
    doit()
