from utils import *
from constants import *

from sympy import *
import pprint


def calculate_Casimir_Effect(diameter):
    P10nm, k = symbols(
        'P10nm, k',
        real=True
    )

    # Volume of 4D ball = 1/2*pi^2*radius^4
    eq = [
        Eq(P10nm,
           1 / 2 * c * h * Sum(1 / k ** 4, (k, 1, diameter / lp)).doit() / (
                   1 / 2 * pi ** 2 * (diameter * 2) ** 4) * eV2J),
    ]
    printE(solve(eq, P10nm, dict=True)[0][P10nm])
    Area = pi * (lp / diameter) ** 2
    Volume = diameter ** 3
    eq = [
        Eq(P10nm,
           1 / 2 * c * h * Sum(k, (k, 0, diameter / lp)).doit() / diameter * zeta(-3) * Area / Volume * eV2J),
    ]
    printE(solve(eq, P10nm, dict=True)[0][P10nm])
    printE(1 / 2 * c * h * (lp + 0) * lp / 2 / lp ** 3 * zeta(-3) * eV2J)
    return h * c * float(pi) / 480 / diameter ** 4 * eV2J


def calculate_ZPE(diameter, min_length=lp):
    ZPE, k = symbols('ZPE, k', real=True, positive=True)
    eq = [
        Eq(ZPE, Sum(k, (k, 1, diameter / min_length)).doit() / diameter * c * h / 2 / (
                4 / 3 * pi * (diameter / 2) ** 3) * eV2J)
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


def calculate_Vacuum_Energy(m):
    energy_density = symbols('energy_density', real=True)
    eq = [
        Eq(energy_density, 2 * m ** 4 * c ** 5 / hbar ** 3)
    ]
    res = solve(eq, energy_density, dict=True)[0]
    print(res)
    return res[energy_density]


def calculate_Power_Density(wave_length):
    frequency = 1 / wave_length
    power_density, k = symbols('power_density, k', real=True)
    eq = [
        Eq(power_density, frequency ** 5 * h / pi ** 4 / c * eV2J)
    ]
    res = solve(eq, power_density, dict=True)[0]
    return res[power_density]


def calculate_Reaction_Probability(wave_length):
    frequency = 1 / wave_length
    probability, k = symbols('probability, k', real=True)
    eq = [
        Eq(probability, h*c*pi**3*Sum(1/(k*frequency)**3, (k, 1, oo)).doit())
    ]
    res = solve(eq, probability, dict=True)[0]
    return res[probability]


def calculate_Power(wave_length):
    frequency = 1 / wave_length
    power, k = symbols('power, k', real=True)
    eq = [
        Eq(power, frequency ** 2 * h / (pi * c))
    ]
    res = solve(eq, power, dict=True)[0]
    return res[power]


def distance_from_redshift(redshift):
    distance = symbols('distance', real=True)
    eq = [
        Eq(distance, redshift * c / H0)
    ]
    res = solve(eq, distance, dict=True)
    return res[0][distance]


def distance_from_redshift_and_P_vac(redshift, P_vac):
    distance = symbols('distance', real=True, positive=True)
    eq = [
        # 9.243700e+17 ~= (c/pi)**2
        Eq(P_vac * distance, redshift * 9.243700e+17),
    ]
    res = solve(eq, distance, dict=True)
    return res[0][distance]


def distance_from_redshift_and_hubble(redshift, hubble):
    distance = symbols('distance', real=True, positive=True)
    eq = [
        Eq(distance, c * redshift / hubble)
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

    # J.s-1.m-3 功率密度
    """
    Critical density of the universe
    """
    printE("Critical density", 3 * (H0) ** 2 / 8 / pi / G * kg2J, "J/m3")
    """
    E_vac/V_vac ∝ photon/L^3, L ∝ V^0.5, h*c/L/V = h*c*L^-0.5 = h*c^0.5 if L = c
    """
    printE("Redshift Energy", h * c ** 6 * pi / 2 * eV2J * SN_2008D['distance'], "J")
    printE("Redshift Energy Density", SN_2006gy['energy_redshift'] / SN_2006gy['distance'] / (4 / 3 * pi * RU ** 3),
           "J/m3")
    printE("Redshift Power", SN_2006gy['energy_redshift'] * c / SN_2006gy['distance'], "J/s")
    printE("Redshift Power Density", SN_2006gy['energy_redshift'] * c / SN_2006gy['distance'] / (4 / 3 * pi * RU ** 3), "J/m3s")
    """
    P_vac = E*Area/Volume = (1/2)*h*c/λ * (pi*r^2) / (4/3*pi*r^3) => (if λ = c^-1,r = c^-1) => (3/8)*h*c^3 * ev2J
    """
    R_ved = c ** -1
    P_vac = (1 / 2) * h * c / R_ved * (pi * R_ved ** 2) / (4 / 3 * pi * R_ved ** 3) * eV2J
    printE("Guessed Vacuum Energy Density", P_vac, "J/m3")
    """
    = printE((3 / 8) * h * c ** 3 * eV2J)
    = printE(-calculate_Casimir_Effect(2.09922509e-5))
    == 6.694982e-09
    """
    P_vac_A = P_vac * (pi * R_ved ** 2)
    printE(P_vac_A)
    """
    Assume Density of matter:Density of vacuum energy = 2.340226e-25:2.340226e-25 = 1:1
    """
    # New Hubble constant
    H0_new = sqrt(16 * pi * G / 3 * P_vac_A)
    printE("New Hubble constant", H0_new)
    printE("Critical density of new Hubble constant", 3 * H0_new ** 2 / 8 / pi / G * kg2J, "J/m3")

    # 量子漲落時間 Qt = hbar/2/Energy
    # c * Qt = radius
    # 測算真空能量密度： 5.18e44 J/m3
    # 中微子能量 < 0.8 eV 約等同波長 1.55e-6 的光
    # 能量總和 1.675462e+04
    # Wave length =  h * c / Neutrino mass
    # Neutrino < 0.8 ev => Wave length > 1.55e-6
    printE("Neutrino's wave length", h * c / 0.8, "m")
    printE("Power density of QF to create Neutrino", calculate_Power_Density(1.55e-6), "J/m3s")
    printE("Power density of QF to create X", calculate_Power_Density(1.6e-9), "J/m3s")
    printE(h * c / 1.6e-9)
    printE("Power of QF to create X", calculate_Power(1.6e-9), "J/s")
    printE("Probability of QF to create X", calculate_Reaction_Probability(c), "J")


if __name__ == '__main__':
    doit()
