import numpy as np


def calc_optimal_chi(kmm, gammastar, ns_star, ca, vpd, beta):
    # -----------------------------------------------------------------------
    # Input:    - float, 'kmm' : Pa, Michaelis-Menten coeff.
    #           - float, 'ns_star'  : (unitless) viscosity correction factor for water
    #           - float, 'vpd' : Pa, vapor pressure deficit
    # Output:   float, ratio of ci/ca (chi)
    # Features: Returns an estimate of leaf internal to ambient CO2
    #           partial pressure following the "simple formulation".
    # Depends:  - kc
    #           - ns
    #           - vpd
    # -----------------------------------------------------------------------

    # leaf-internal-to-ambient CO2 partial pressure (ci/ca) ratio
    xi = np.sqrt((beta * (kmm + gammastar)) / (1.6 * ns_star))
    chi = gammastar / ca + (1.0 - gammastar / ca) * xi / (xi + np.sqrt(vpd))

    # Define variable substitutes:
    vdcg = ca - gammastar
    vacg = ca + 2.0 * gammastar
    vbkg = beta * (kmm + gammastar)

    # wrap if condition in a function to allow vectorization

    def calc_mj(ns_star, vpd, vbkg):
        vsr = np.sqrt(1.6 * ns_star * vpd / vbkg)
        # Based on the mc' formulation (see Regressing_LUE.pdf)
        mj = vdcg / (vacg + 3.0 * gammastar * vsr)
        return mj

    # Check for negatives, vectorized
    if ns_star > 0 and vpd > 0 and vbkg > 0:
        mj = calc_mj(ns_star, vpd, vbkg)
    elif type(vpd) == 'float':
        mj = np.NaN
    else:
        mj = np.zeros(len(vpd))

    # alternative variables
    gamma = gammastar / ca
    kappa = kmm / ca

    # mc
    mc = (chi - gamma) / (chi + kappa)

    # mj:mv
    mjoc = (chi + kappa) / (chi + 2 * gamma)

    out = {'chi': chi, 'mc': mc, 'mj': mj, 'mjoc': mjoc}
    return out



def calc_lue_vcmax_wang17(out_optchi, kphio, ftemp_kphio, c_molmass, soilmstress):

    # Include effect of Jmax limitation
    if isinstance(out_optchi['chi'], float):
        leng = 1
    else:
        leng = len(out_optchi[[1]])
    mprime = calc_mprime(out_optchi['mj'])

    out ={
    # Light use efficiency (gpp per unit absorbed light)
    'lue':kphio * ftemp_kphio * mprime * c_molmass * soilmstress,

    # Vcmax normalised per unit absorbed PPFD (assuming iabs=1), with Jmax limitation
    'vcmax_unitiabs':kphio * ftemp_kphio * out_optchi['mjoc'] * mprime / out_optchi['mj'] * soilmstress,

    # complement for non-smith19
    'omega' : np.full(leng, np.NaN),
    'omega_star': np.full(leng, np.NaN)
    }
    return out



def calc_lue_vcmax_smith19(out_optchi, kphio, ftemp_kphio, c_molmass, soilmstress):


    if isinstance(out_optchi['chi'], float):
        leng = 1
    else:
        leng = len(out_optchi[[1]])

    # Adopted from Nick Smith's code:
    # Calculate omega, see Smith et al., 2019 Ecology Letters
    def calc_omega(theta, c_cost, m):

        cm = 4 * c_cost / m  # simplification term for omega calculation
        v =1 / (cm * (1 - theta * cm)) - 4 * theta  # simplification term for omega calculation

        # account for non-linearities at low m values
        capP = (((1 / 1.4) - 0.7) **2 / (1 - theta)) + 3.4
        aquad = -1
        bquad = capP
        cquad = -(capP * theta)
        m_star = (4 * c_cost) / np.polynomial.polynomial.polyroots([aquad, bquad, cquad])

        if m < np.real(m_star[0]):
            omega = -(1 - (2 * theta)) - np.sqrt((1 - theta) * v)
        else:
            omega = -(1 - (2 * theta)) + np.sqrt((1 - theta) * v)
        return (omega)


    ## constants
    theta = 0.85  # should be calibratable?
    c_cost = 0.05336251

    ## factors derived as in Smith et al., 2019
    omega = calc_omega(theta=theta, c_cost=c_cost, m=out_optchi['mc'] )  # Eq. S4
    omega_star = 1.0 + omega - np.sqrt((1.0 + omega)**2 - (4.0 * theta * omega))  # Eq. 18

    ## Effect of Jmax limitation
    mprime = out_optchi['mj']  * omega_star / (8.0 * theta)

    ## Light use efficiency (gpp per unit absorbed light)
    lue = kphio * ftemp_kphio * mprime * c_molmass * soilmstress

    # calculate Vcmax per unit aborbed light
    vcmax_unitiabs = kphio * ftemp_kphio * out_optchi['mjoc']  * omega_star / (8.0 * theta) * soilmstress  # Eq. 19

    out ={
    'lue':  lue,
    'vcmax_unitiabs' : vcmax_unitiabs,
    'omega':  omega,
    'omega_star' : omega_star
    }
    return out



def calc_lue_vcmax_none(out_optchi, kphio, ftemp_kphio, c_molmass, soilmstress):
    ## Do not include effect of Jmax limitation
    if isinstance(out_optchi['chi'], float):
        leng = 1
    else:
        leng = len(out_optchi[[1]])

    out ={
    # Light use efficiency (gpp per unit absorbed light)
    'lue': kphio * ftemp_kphio * out_optchi['mj'] * c_molmass * soilmstress,

    # Vcmax normalised per unit absorbed PPFD (assuming iabs=1), with Jmax limitation
    'vcmax_unitiabs': kphio * ftemp_kphio * out_optchi['mjoc'] * soilmstress,

    # complement for non-smith19
    'omega': np.full(leng, np.NaN),
    'omega_star': np.full(leng, np.NaN)
     }
    return out


def calc_lue_vcmax_c4(kphio, ftemp_kphio, c_molmass, soilmstress):

    if isinstance(kphio, float):
        leng = 1
    else:
        leng = len(kphio)
    out ={
    # Light use efficiency (gpp per unit absorbed light)
    'lue': kphio * ftemp_kphio * c_molmass * soilmstress,

    ## Vcmax normalised per unit absorbed PPFD (assuming iabs=1), with Jmax limitation
    'vcmax_unitiabs':kphio * ftemp_kphio * soilmstress,

    ## complement for non-smith19
    'omega':np.full(leng, np.NaN),
    'omega_star':np.full(leng, np.NaN)
    }
    return out


def calc_chi_c4():
    # //////////////////////////////////////////////////////////////////
    # (Dummy-) ci:ca for C4 photosynthesis
    # -----------------------------------------------------------------------
    out ={'chi':9999, 'mc':1, 'mj':1, 'mjoc':1 }
    return out


def calc_mprime(mc):
    # -----------------------------------------------------------------------
    # Input:  mc   (unitless): factor determining LUE
    # Output: mpi (unitless): modiefied m accounting for the co-limitation
    #                         hypothesis after Prentice et al. (2014)
    # -----------------------------------------------------------------------
    kc = 0.41  # Jmax cost coefficient

    mpi = mc**2 - kc**(2.0 / 3.0) * (mc**(4.0 / 3.0))

    # Check for negatives:
    if mpi>0:
        mpi = np.sqrt(mpi)
    else:
        mpi=np.NaN
    return mpi


def co2_to_ca(co2, patm):

    # -----------------------------------------------------------------------
    # Input:    - float, annual atm. CO2, ppm (co2)
    #           - float, monthly atm. pressure, Pa (patm)
    # Output:   - ca in units of Pa
    # Features: Converts ca (ambient CO2) from ppm to Pa.
    # -----------------------------------------------------------------------
    ca = 1.0e-6 * co2 * patm  # Pa, atms. CO2
    return ca


def density_h2o(tc, p):

    # -----------------------------------------------------------------------
    # Input:    - float, air temperature (tc), degrees C
    #           - float, atmospheric pressure (p), Pa
    # Output:   float, density of water, kg/m^3
    # Features: Calculates density of water at a given temperature and
    #           pressure using the Tumlirz Equation
    # Ref:      F.H. Fisher and O.E Dial, Jr. (1975) Equation of state of
    #           pure water and sea water, Tech. Rept., Marine Physical
    #           Laboratory, San Diego, CA.
    # -----------------------------------------------------------------------

    # Calculate lambda, (bar cm^3)/g:
    my_lambda = 1788.316 + 21.55053 * tc + -0.4695911 * tc * tc + 3.096363e-3 * tc * tc * tc + \
                -7.341182e-6 * tc * tc * tc * tc

    # Calculate po, bar
    po = 5918.499 + 58.05267 * tc + -1.1253317 * tc * tc + 6.6123869e-3 * tc * tc * tc + \
         -1.4661625e-5 * tc * tc * tc * tc

    # Calculate vinf, cm^3/g
    vinf = 0.6980547 + -7.435626e-4 * tc + 3.704258e-5 * tc * tc + -6.315724e-7 * tc * tc * tc + \
           9.829576e-9 * tc * tc * tc * tc + -1.197269e-10 * tc * tc * tc * tc * tc + \
           1.005461e-12 * tc * tc * tc * tc * tc * tc + \
           -5.437898e-15 * tc * tc * tc * tc * tc * tc * tc + \
           1.69946e-17 * tc * tc * tc * tc * tc * tc * tc * tc + \
           -2.295063e-20 * tc * tc * tc * tc * tc * tc * tc * tc * tc

    # Convert pressure to bars (1 bar <- 100000 Pa)
    pbar = 1e-5 * p

    # Calculate the specific volume (cm^3 g^-1):
    v = vinf + my_lambda / (po + pbar)

    # Convert to density (g cm^-3) -> 1000 g/kg; 1000000 cm^3/m^3 -> kg/m^3:
    rho = (1e3 / v)

    return rho


def calc_viscosity_h2o(tc, p):

    # -----------------------------------------------------------------------
    # Input:    - float, ambient temperature (tc), degrees C
    #           - float, ambient pressure (p), Pa
    # Return:   float, viscosity of water (mu), Pa s
    # Features: Calculates viscosity of water at a given temperature and
    #           pressure.
    # Depends:  density_h2o
    # Ref:      Huber, M. L., R. A. Perkins, A. Laesecke, D. G. Friend, J. V.
    #           Sengers, M. J. Assael, ..., K. Miyagawa (2009) New
    #           international formulation for the viscosity of H2O, J. Phys.
    #           Chem. Ref. Data, Vol. 38(2), pp. 101-125.
    # -----------------------------------------------------------------------

    # Define reference temperature, density, and pressure values:
    tk_ast = 647.096  # Kelvin
    rho_ast = 322.0  # kg/m^3
    mu_ast = 1e-6  # Pa s

    # Get the density of water, kg/m^3
    rho = density_h2o(tc, p)

    # Calculate dimensionless parameters:
    tbar = (tc + 273.15) / tk_ast
    tbarx = tbar**0.5
    tbar2 = tbar**2
    tbar3 = tbar**3
    rbar = rho / rho_ast

    # Calculate mu0 (Eq. 11 & Table 2, Huber et al., 2009):
    mu0 = 1.67752 + 2.20462 / tbar + 0.6366564 / tbar2 - 0.241605 / tbar3
    mu0 = 1e2 * tbarx / mu0

    # Create Table 3, Huber et al. (2009):
    h_array = np.zeros([7, 6])
    h_array[0, :] = [0.520094, 0.0850895, -1.08374, -0.289555, 0.0, 0.0]  # hj0
    h_array[1, :] = [0.222531, 0.999115, 1.88797, 1.26613, 0.0, 0.120573]  # hj1
    h_array[2, :] = [-0.281378, -0.906851, -0.772479, -0.489837, -0.257040, 0.0]  # hj2
    h_array[3, :] = [0.161913, 0.257399, 0.0, 0.0, 0.0, 0.0]  # hj3
    h_array[4, :] = [-0.0325372, 0.0, 0.0, 0.0698452, 0.0, 0.0]  # hj4
    h_array[5, :] = [0.0, 0.0, 0.0, 0.0, 0.00872102, 0.0]  # hj5
    h_array[6, :] = [0.0, 0.0, 0.0, -0.00435673, 0.0, -0.000593264]  # hj6

    # Calculate mu1 (Eq. 12 & Table 3, Huber et al., 2009):
    mu1 = 0.0
    ctbar = (1.0 / tbar) - 1.0
    # print(paste("ctbar",ctbar))
    # for i in xrange(6):
    for i in range(0, 6):
        coef1 = ctbar**(i - 1)
    # print(paste("i, coef1", i, coef1))
        coef2 = 0.0
        for j in range(0, 7):
            coef2 = coef2 + h_array[j, i] * (rbar - 1.0) ** (j - 1)
        mu1 = mu1 + coef1 * coef2
    mu1 = np.exp(rbar * mu1)
    # print(paste("mu1",mu1))

    # Calculate mu_bar (Eq. 2, Huber et al., 2009)
    #   assumes mu2 = 1
    mu_bar = mu0 * mu1

    # Calculate mu (Eq. 1, Huber et al., 2009)
    mu = mu_bar * mu_ast  # Pa s

    return mu