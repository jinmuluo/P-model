import numpy as np


def calc_patm(elv):

    kTo = 298.15  # base temperature, K
    kL = 0.0065  # adiabiatic temperature lapse rate, K/m
    kG = 9.80665  # gravitational acceleration, m/s^2
    kR = 8.3145  # universal gas constant, J/mol/K
    kMa = 0.028963  # molecular weight of dry air, kg/mol
    patm0 = 101325

    patm = patm0 * (1.0 - kL * elv / kTo)**(kG * kMa / (kR * kL))
    return patm


def calc_ftemp_kphio(tc):

    ftemp_kphio = 0.352 + 0.022 * tc - 3.4e-4 * tc**2

    return ftemp_kphio


def calc_soilmstress(soilm, meanalpha, apar_soilm, bpar_soilm):

    # Fixed parameters
    x0 = 0.0
    x1 = 0.6
    if soilm > x1:
        outstress = 1.0
    else:
        y0 = (apar_soilm + bpar_soilm * meanalpha)
        beta = (1.0 - y0) / (x0 - x1)**2
        outstress = 1.0 - beta * (soilm - x1)**2
    outstress = max(0.0, min(1.0, outstress))
    return outstress


def calc_ftemp_arrh( tk, dha, tkref=298.15):
    # Note that the following forms are equivalent:
    # ftemp = exp( dha * (tk - 298.15) / (298.15 * kR * tk) )
    # ftemp = exp( dha * (tc - 25.0)/(298.15 * kR * (tc + 273.15)) )
    # ftemp = exp( (dha/kR) * (1/298.15 - 1/tk) )
    # -----------------------------------------------------------------------
    kR = 8.3145     # Universal gas constant, J/mol/K
    ftemp = np.exp(dha * (tk - tkref) / (tkref * kR * tk))
    return ftemp


def calc_gammastar(tc, patm):
    # -----------------------------------------------------------------------
    # Input:    float, air temperature, degrees C (tc)
    # Output:   float, gamma-star, Pa (gammastar)
    # Features: Returns the temperature-dependent photorespiratory
    #           compensation point, Gamma star (Pascals), based on constants
    #           derived from Bernacchi et al. (2001) study.
    # Ref:      Bernacchi et al. (2001), Improved temperature response
    #           functions for models of Rubisco-limited photosynthesis,
    #           Plant, Cell and Environment, 24, 253--259.
    # -----------------------------------------------------------------------
    dha = 37830
    gs25_0 = 4.332
    gammastar = gs25_0 * patm / calc_patm(0.0) * calc_ftemp_arrh((tc + 273.15), dha)
    return gammastar


def calc_kmm(tc, patm):
    dhac = 79430  # (J/mol) Activation energy, Bernacchi et al. (2001)
    dhao = 36380  # (J/mol) Activation energy, Bernacchi et al. (2001)
    kco = 2.09476e5  # (ppm) O2 partial pressure, Standard Atmosphere

    # k25 parameters are not dependent on atmospheric pressure
    kc25 = 39.97  # Pa,converted to Pa by T. Davis assuming elevation of 227.076 m.a.s.l.
    ko25 = 27480  # Pa,converted to Pa by T. Davis assuming elevation of 227.076 m.a.s.l.

    # conversion to Kelvin
    tk = tc + 273.15
    kc = kc25 * calc_ftemp_arrh(tk, dhac)
    ko = ko25 * calc_ftemp_arrh(tk, dhao)
    po = kco * 1e-6 * patm  # O2 partial pressure
    kmm = kc * (1.0 + po / ko)
    return kmm


def calc_ftemp_inst_vcmax( tcleaf, tcgrowth = 0.001, tcref =25.0):
    if tcgrowth == 0.001:
        tcgrowth = tcleaf
    # loal parameters
    Ha = 71513  # activation energy (J/mol)
    Hd = 200000  # deactivation energy (J/mol)
    Rgas = 8.3145  # universal gas constant (J/mol/K)
    a_ent = 668.39  # offset of entropy vs. temperature relationship from Kattge & Knorr (2007) (J/mol/K)
    b_ent = 1.07  # slope of entropy vs. temperature relationship from Kattge & Knorr (2007) (J/mol/K^2)

    tkref = tcref + 273.15  # to Kelvin

    # conversion of temperature to Kelvin, tcleaf is the instantaneous leaf temperature in degrees C.
    tkleaf = tcleaf + 273.15

    # calculate entropy following Kattge & Knorr (2007), negative slope and y-axis intersect is when expressed as a
    # function of temperature in degrees Celsius, not Kelvin !!!
    dent = a_ent - b_ent * tcgrowth  # 'tcgrowth' corresponds to 'tmean' in Nicks, 'tc25' is 'to' in Nick's
    fva = calc_ftemp_arrh(tkleaf, Ha, tkref=tkref)
    fvb = (1 + np.exp((tkref * dent - Hd) / (Rgas * tkref))) / (1 + np.exp((tkleaf * dent - Hd) / (Rgas * tkleaf)))
    fv = fva * fvb
    return (fv)


def calc_ftemp_inst_rd( tc ):
    # loal parameters
    apar = 0.1012
    bpar = 0.0005
    fr = np.exp(apar * (tc - 25.0) - bpar * (tc**2 - 25.0**2))
    return fr
