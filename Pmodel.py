import numpy as np
import calc
import fuc


# ----------------------------------------------------------------------------------------------------------------------
# MAIN FUNCTION DEFINITION
# ----------------------------------------------------------------------------------------------------------------------
def pypmodel(tc, vpd, co2, fapar, ppfd, patm=np.NaN, elv=np.NaN, beta=146.0, soilm=1.0,
              meanalpha=1.0, apar_soilm=0.0, bpar_soilm=0.685, c4='FALSE', method_optci='prentice14',
              method_jmaxlim="wang17", do_ftemp_kphio='TRUE', do_soilmstress='FALSE', returnvar=None,
              verbose = 'FALSE'):
    # ------------------------------------------------------------------------------------------
    # kphio is determined by the  the simulation type.
    # ------------------------------------------------------------------------------------------
    if do_ftemp_kphio.lower() == 'false':
        kphio = 0.049977
    elif do_soilmstress.lower() == 'true':
        kphio = 0.087132
    else:
        kphio = 0.081785

    # ------------------------------------------------------------------------------------------
    # Fixed parameters
    # ------------------------------------------------------------------------------------------
    if elv == np.NaN and patm == np.NaN:
        print('Warning:Both elevation and atmospheric pressure are absenting!')
        exit('Please provide one of them!')
    elif elv != np.NaN:
        print('Found elevation but not atmospheric pressure, calculating it as a function of elevation (elv), ')
        print('assuming standard atmosphere (101325 Pa at sea level)')
        patm = calc.calc_patm(elv)

    # ------------------------------------------------------------------------------------------
    # Temperature dependence of quantum yield efficiency
    # ------------------------------------------------------------------------------------------
    c_molmass = 12.0107                 # molecular mass of carbon, g
    kPo = 101325.0                      # standard atmosphere pressure, Pa
    kTo = 25.0                          # base temperature, Celsius
    rd_to_vcmax = 0.015                 # Ratio of Rdark to Vcmax25, number from Atkin et al., 2015 for C3 herbaceous

    # ------------------------------------------------------------------------------------------
    # Temperature dependence of quantum yield efficiency
    # 'do_ftemp_kphio' is not actually a stress function, but is the temperature-dependency of
    # the quantum yield efficiency after Bernacchi et al., 2003 PCE
    # ------------------------------------------------------------------------------------------
    if do_ftemp_kphio.lower() == 'true':
        ftemp_kphio = calc.calc_ftemp_kphio(tc)
        print('(*•̀ㅂ•́)و Using temperature-dependency of the quantum yield efficiency')
    else:
        ftemp_kphio = 1.0
        print('╮(๑•́ ₃•̀๑) Using constant quantum yield efficiency:', ftemp_kphio)

    # ------------------------------------------------------------------------------------------
    # Calculate soil moisture stress as a function of soil moisture and mean alpha
    # ------------------------------------------------------------------------------------------
    if do_soilmstress.lower() == 'true':
        soilmstress = calc.calc_soilmstress(soilm, meanalpha, apar_soilm, bpar_soilm)
        print('(*•̀ㅂ•́)و Calculate soil moisture stress as a function of soil moisture and mean alpha.')
    else:
        soilmstress = 1.0
        print('╮(๑•́ ₃•̀๑)╭ Using constant soil moisture stress:', soilmstress)

    # -------------------------------------------------------------------------------------------
    # Photosynthesis model parameters depending on temperature, pressure, and CO2.
    # -------------------------------------------------------------------------------------------
    # ambient CO2 partial pressure(Pa)
    ca = fuc.co2_to_ca(co2, patm)

    # photo-respiratory compensation point - Gamma-star (Pa)
    gammastar = calc.calc_gammastar(tc, patm)

    # Michaelis-Menten coef. (Pa)
    kmm = calc.calc_kmm(tc, patm)  # replace 'NA' here with 'patm'

    # viscosity correction factor = viscosity( temp, press )/viscosity( 25 degC, 1013.25 Pa)
    ns = fuc.calc_viscosity_h2o(tc, patm)  # Pa s
    ns25 = fuc.calc_viscosity_h2o(kTo, kPo)  # Pa s
    ns_star = ns / ns25  # (unitless)

    # -----------------------------------------------------------------------
    # Optimal ci
    # The heart of the P-model: calculate ci:ca ratio (chi) and additional terms
    # -----------------------------------------------------------------------
    out_optchi = {}
    if c4.lower() == 'true':
        # "dummy" ci:ca for C4 plants
        out_optchi = fuc.calc_chi_c4()

    elif method_optci.lower() == "prentice14":
        # Full formualation (Gamma-star not zero), analytical solution
        out_optchi = fuc.calc_optimal_chi(kmm, gammastar, ns_star, ca, vpd, beta)
    else:
        print("rpmodel(): argument method_optci not idetified, program end abnormally!.")
        exit(0)

    # leaf-internal CO2 partial pressure (Pa)
    ci = out_optchi['chi'] * ca

    # -----------------------------------------------------------------------
    # Corrolary preditions
    # -----------------------------------------------------------------------
    # intrinsic water use efficiency (in Pa)
    iwue = (ca - ci) / 1.6

    # -----------------------------------------------------------------------
    # Vcmax and light use efficiency
    # -----------------------------------------------------------------------
    if c4.lower() == 'true':
        out_lue_vcmax = fuc.calc_lue_vcmax_c4(kphio, ftemp_kphio, c_molmass, soilmstress)

    elif method_jmaxlim.lower() == "wang17":

        out_lue_vcmax = fuc.calc_lue_vcmax_wang17(out_optchi, kphio, ftemp_kphio, c_molmass, soilmstress)

    elif method_jmaxlim.lower()  == "smith19":

        out_lue_vcmax = fuc.calc_lue_vcmax_smith19(out_optchi, kphio, ftemp_kphio, c_molmass, soilmstress)

    elif method_jmaxlim.lower()  == "none":

        out_lue_vcmax = fuc.calc_lue_vcmax_none(out_optchi, kphio, ftemp_kphio, c_molmass, soilmstress)

    else:
        print("rpmodel(): argument method_jmaxlim not idetified.")
        exit('Please choose one method.')

    # -----------------------------------------------------------------------
    #  Corrolary preditions
    # -----------------------------------------------------------------------
    #  Vcmax25 (vcmax normalized to 25 deg C)
    ftemp25_inst_vcmax = calc.calc_ftemp_inst_vcmax(tc, tc, tcref=25.0)
    vcmax25_unitiabs = out_lue_vcmax['vcmax_unitiabs'] / ftemp25_inst_vcmax

    #  Dark respiration at growth temperature
    ftemp_inst_rd = calc.calc_ftemp_inst_rd(tc)
    rd_unitiabs = rd_to_vcmax * (ftemp_inst_rd / ftemp25_inst_vcmax) * out_lue_vcmax['vcmax_unitiabs']

    # -----------------------------------------------------------------------
    #  Quantities that scale linearly with absorbed light
    # -----------------------------------------------------------------------
    if isinstance(out_lue_vcmax['lue'], float):
        leng = 1
    else:
        leng = len(out_optchi[[1]])
    iabs = np.full(leng, fapar * ppfd)
    # Gross primary productivity, in g C m-2 s-1
    if np.any(np.isnan(iabs)) != 0:
        gpp = iabs * out_lue_vcmax['lue']
    else:
        gpp = np.full(leng, np.NaN)

    #  Vcmax per unit ground area is the product of the intrinsic quantum
    #  efficiency, the absorbed PAR, and 'n'
    if np.any(np.isnan(iabs)) != 0:
        vcmax = iabs * out_lue_vcmax['vcmax_unitiabs']
    else:
        vcmax = np.full(leng, np.NaN)

    # (vcmax normalized to 25 deg C)
    if np.any(np.isnan(iabs)) != 0:
        vcmax25 = iabs * out_lue_vcmax['vcmax25_unitiabs']
    else:
        vcmax25 = np.full(leng, np.NaN)

    # Dark respiration
    if np.any(np.isnan(iabs)) != 0:
        rd = iabs * rd_unitiabs
    else:
        rd = np.full(leng, np.NaN)

    # construct list for output
    out = {
        'ca':np.full(leng, ca),
        'gammastar':np.full(leng, gammastar),
        'kmm':np.full(leng, kmm),
        'ns_star':np.full(leng, ns_star),
        'chi':out_optchi['chi'],
        'mj':out_optchi['mj'],
        'mc': out_optchi['mc'],
        'ci': ci,
        'lue':out_lue_vcmax['lue'],
        'gpp': gpp,
        'iwue': iwue,
        'gs': (gpp / c_molmass) / (ca - ci),
        'vcmax': vcmax,
        'vcmax25': vcmax25,
        'rd': rd
    }

    if returnvar:
        for i in range(len(returnvar)):
            out[i] = out[returnvar[i]]
    return out



