
# line 12 - 212 parameter introduction
# line 213 - end references paper

# ----------------------------------------------------------------------------------------------------------------------
#  Detail information about each function used in P-MODEL
# ----------------------------------------------------------------------------------------------------------------------
#
#
#
#
# ----------------------------------------------------------------------------------------------------------------------
#  Parameters detail information in all function, the format all follow the below write method.
# ----------------------------------------------------------------------------------------------------------------------
# PARAMETER NAME      UNITS          OPTIONAL/NONE         PHYSICAL DEFINITION
# co2                 ppm            none                  Atmospheric CO2 concentration
# fapar               unitless       optional              Fraction of absorbed photosynthetically active radiation
# tc                  celsius        none                  Temperature, relevant for photosynthesis
# vpd                 pa             none                  Vapour pressure deficit
# ppfd                mol/m^2/t      optional              Photosynthetic photon flux density .
"""
   Note that the units of {ppfd} (per area and per time) determine the units of outputs {lue}, {gpp}, {vcmax}, and {rd}.
   For example, if {ppfd} is provided in units of mol m-2 month-1, then respective output variables are returned as per
   unit months.
"""

# patm                pa              none                 Atmospheric pressure
# elv                 m               none                 Elevation above sea-level (m.a.s.l.).
"""
  When {patm} provided, overrides parameter {elv}, otherwise {patm} is calculated using standard atmosphere (101325 Pa), 
  corrected for elevation {elv}, using the function \link{calc_patm}.
"""

# kphio               unitless        none                 Apparent quantum yield efficiency, Defaults to 0.0817
"""
  \code{method_jmaxlim="wang17", do_ftemp_kphio=TRUE, do_soilmstress=FALSE}, 0.0870 for
  \code{method_jmaxlim="wang17", do_ftemp_kphio=TRUE, do_soilmstress=TRUE}, and 0.0492 for
  \code{method_jmaxlim="wang17", do_ftemp_kphio=FALSE, do_soilmstress=FALSE}, corresponding to the empirically 
   fitted value as presented in Stocker et al. (2019) Geosci. Model Dev. for model setup 'BRC', 'FULL', and 'ORG'
   respectively.
"""

# beta                unitless        none                  Unit cost ratio. Defaults to 146 (see Stocker et al., 2019).
# soilm               fraction        optional              Relative soil moisture as a fraction
"""
  only if {do_soilmstress==TRUE}) 
  of field capacity (unitless). Defaults to 1.0 (no soil moisture stress). This information is used to calculate
  an empirical soil moisture stress factor (\link{calc_soilmstress}) whereby the sensitivity is determined
  by average aridity, defined by the local annual mean ratio of actual over potential evapotranspiration,
  supplied by argument \code{meanalpha}.
"""

# param               unitless        optional               Local annual mean ratio of actual over potential
#                                                            evapotranspiration, measure for average aridity.
#                                                            Defaults to 1.0
"""
 mean alpha (used only if \code{do_soilmstress==TRUE})
"""

# apar_soilm          unitless        optional               Parameter determining the sensitivity of the empirical soil
#                                                            moisture stress function
"""
  used only if \code{do_soilmstress==TRUE}, Defaults to 0.0, the empirically fitted value
  as presented in Stocker et al. (2019) Geosci. Model Dev. for model setup 'FULL' (corresponding to a setup
  with \code{method_jmaxlim="wang17", do_ftemp_kphio=TRUE, do_soilmstress=TRUE}).
"""

# bpar_soilm          unitless        optional               Parameter determining the sensitivity of the empirical soil
#                                                            moisture stress function
"""
  used only if \code{do_soilmstress==TRUE}) . Defaults to 0.685, the empirically fitted value
  as presented in Stocker et al. (2019) Geosci. Model Dev. for model setup 'FULL' (corresponding to a setup
  with \code{method_jmaxlim="wang17", do_ftemp_kphio=TRUE, do_soilmstress=TRUE}).
"""

# c4                  unitless        optional                A logical value specifying whether the C3 or C4
#                                                             photosynthetic pathway is followed.
"""
  Defaults to \code{FALSE}. If \code{TRUE}, the leaf-internal CO2 concentration is assumed to be very large
  and \eqn{m} (returned variable \code{mj}) tends to 1, and \eqn{m'} tends to 0.669 (with \code{c = 0.41}).
"""

# method_optci        unitless        optional                A character string specifying which method is to be used
#                                                             for calculating
"""
  optimal ci:ca. Defaults to \code{"prentice14"}.
  Available also \code{"prentice14_num"} for a numerical solution to the same optimization criterium as
  used for \code{"prentice14"}.
"""

# method_jmaxlim      unitless        optional                 A character string specifying which method is to be used
#                                                              for factoring in
"""
  Jmax limitation. Defaults to \code{"wang17"},
  based on Wang Han et al. 2017 Nature Plants and (Smith 1937). Available is also \code{"smith19"}, following
  the method by Smith et al., 2019 Ecology Letters,
  and \code{"none"} for ignoring effects of Jmax limitation.
"""

# do_ftemp_kphio       unitless       optional                  A logical specifying whether temperature-dependence of
#                                                               quantum yield
"""
  efficiency after Bernacchi et al., 2003 is to be accounted for. Defaults to \code{TRUE}.
"""

# do_soilmstress       unitless       optional                  A logical specifying whether an empirical soil moisture
#                                                               stress factor
"""
  is to be applied to down-scale light use efficiency (and only light use efficiency). Defaults to \code{FALSE}.
"""

# returnvar            unitless       optional                  A character string of vector of character strings
#                                                               specifying which variables are to be returned
# verbose              unitless       none                      Logical, defines whether verbose messages are printed.
#                                                               Defaults to \code{FALSE}.
# return               unitless       none                      A named list of numeric values (including temperature
#                                                               and pressure dependent parameters of the photosynthesis
#                                                               model,P-model predictions, including all its corollary).
"""
this includes :
   \itemize{
           \item \code{ca}: Ambient CO2 expressed as partial pressure (Pa)
           \item \code{gammastar}: Photorespiratory compensation point \eqn{\Gamma*}, (Pa), see \link{calc_gammastar}.
           \item \code{kmm}: Michaelis-Menten coefficient \eqn{K} for photosynthesis (Pa), see \link{calc_kmm}.
           \item \code{ns_star}: Change in the viscosity of water, relative to its value at 25 deg C (unitless).
                                 \deqn{\eta* = \eta(T) / \eta(25 deg C)}
                                 This is used to scale the unit cost of transpiration. Calculated following Huber
                                 et al. (2009).
           \item \code{chi}: Optimal ratio of leaf internal to ambient CO2 (unitless). Derived following Prentice et al.
                            (2014) as:
                            \deqn{
                                  \chi = \Gamma* / ca + (1- \Gamma* / ca) \xi / (\xi + \sqrt D )
                            }
                            with
                            \deqn{
                                 \xi = \sqrt (\beta (K+ \Gamma*) / (1.6 \eta*))
                            }
                            \eqn{\beta} is given by argument \code{beta}, \eqn{K} is \code{kmm} (see \link{calc_kmm}),
                            \eqn{\Gamma*} is \code{gammastar} (see \link{calc_gammastar}). \eqn{\eta*} is \code{ns_star}.
                           \eqn{D} is the vapour pressure deficit (argument \code{vpd}), \eqn{ca} is the 
                           '                          ambient CO2 partial pressure in Pa (\code{ca}).
          \item \code{ci}: Leaf-internal CO2 partial pressure (Pa), calculated as \eqn{(\chi ca)}. 
          '         \item \code{lue}: Light use efficiency (g C / mol photons), calculated as
                          \deqn{
                                LUE = \phi(T) \phi0 m' Mc
                           }
                           where \eqn{\phi(T)} is the temperature-dependent quantum yield efficiency modifier
                           (\link{calc_ftemp_kphio}) if \code{do_ftemp_kphio==TRUE}, and 1 otherwise. \eqn{\phi 0}
                           is given by argument \code{kphio}.
                           \eqn{m'=m} if \code{method_jmaxlim=="none"}, otherwise
                           \deqn{
                                  m' = m \sqrt( 1 - (c/m)^(2/3) )
                           }
                           with \eqn{c=0.41} (Wang et al., 2017) if \code{method_jmaxlim=="wang17"}. \eqn{Mc} is
                           the molecular mass of C (12.0107 g mol-1). \eqn{m} is given returned variable \code{mj}.
                           If \code{do_soilmstress==TRUE}, \eqn{LUE} is multiplied with a soil moisture stress factor,
                           calculated with \link{calc_soilmstress}.
           \item \code{mj}: Factor in the light-limited assimilation rate function, given by
                           \deqn{
                               m = (ci - \Gamma*) / (ci + 2 \Gamma*)
                          }
                          where \eqn{\Gamma*} is given by \code{gammastar}.
           \item \code{mc}: Factor in the Rubisco-limited assimilation rate function, given by
                           \deqn{
                               mc = (ci - \Gamma*) / (ci + K)
                          }
                          where \eqn{K} is given by \code{kmm}.
           \item \code{gpp}: Gross primary production (g C m-2), calculated as
                          \deqn{
                              GPP = Iabs LUE
                          }
                          where \eqn{Iabs} is given by \code{fapar*ppfd} (arguments), and is
                          \code{NA} if \code{fapar==NA} or \code{ppfd==NA}. Note that \code{gpp} scales with
                          absorbed light. Thus, its units depend on the units in which \code{ppfd} is given.
           \item \code{iwue}: Intrinsic water use efficiency (iWUE, Pa), calculated as
                          \deqn{
                                iWUE = ca (1-\chi)/(1.6)
                          }
           \item \code{gs}: Stomatal conductance (gs, in mol C m-2 Pa-1), calculated as
                          \deqn{
                               gs = A / (ca (1-\chi))
                          }
                          where \eqn{A} is \code{gpp}\eqn{/Mc}.
           \item \code{vcmax}: Maximum carboxylation capacity \eqn{Vcmax} (mol C m-2) at growth temperature (argument
                          \code{tc}), calculated as
                         \deqn{
                              Vcmax = \phi(T) \phi0 Iabs n
                         }
                         where \eqn{n} is given by \eqn{n=m/mc}, or
                         \deqn{
                             n = (ci + K) / (ci + 2 \Gamma*)
                         }
           \item \code{vcmax25}: Maximum carboxylation capacity \eqn{Vcmax} (mol C m-2) normalised to 25 deg C
                         following a modified Arrhenius equation, calculated as \eqn{Vcmax25 = Vcmax / fv},
                        where \eqn{fv} is the instantaneous temperature response by Vcmax and is implemented
                        by function \link{calc_ftemp_inst_vcmax}.
           \item \code{rd}: Dark respiration \eqn{Rd} (mol C m-2), calculated as
                        \deqn{
                            Rd = b0 Vcmax (fr / fv)
                        }
                        where \eqn{b0} is a constant and set to 0.015 (Atkin et al., 2015), \eqn{fv} is the
                        instantaneous temperature response by Vcmax and is implemented by function
                        \link{calc_ftemp_inst_vcmax}, and \eqn{fr} is the instantaneous temperature response
                        of dark respiration following Heskel et al. (2016) and is implemented by function
                        \link{calc_ftemp_inst_rd}.
   }
"""

# Additional variables are contained in the returned list if argument \code{method_jmaxlim=="smith19"}
# \itemize{
#         \item \code{omega}: Term corresponding to \eqn{\omega}, defined by Eq. 16 in Smith et al. (2019),
#         and Eq. E19 in Stocker et al. (2019).
#         \item \code{omega_star}: Term corresponding to \eqn{\omega^\ast}, defined by Eq. 18 in Smith et al.
#         (2019), and Eq. E21 in Stocker et al. (2019).
#         }
# ----------------------------------------------------------------------------------------------------------------------
#  List of main references
# ----------------------------------------------------------------------------------------------------------------------
# Bernacchi, C. J., Pimentel, C., and Long, S. P.:  In vivo temperature response func-tions  of  parameters
#              required  to  model  RuBP-limited  photosynthesis,  Plant  Cell Environ., 26, 1419–1430, 2003
#
# Heskel,  M.,  O’Sullivan,  O.,  Reich,  P.,  Tjoelker,  M.,  Weerasinghe,  L.,  Penillard,  A.,Egerton, J.,
#          Creek, D., Bloomfield, K., Xiang, J., Sinca, F., Stangl, Z., Martinez-De La Torre, A., Griffin, K.,
#          Huntingford, C., Hurry, V., Meir, P., Turnbull, M.,and Atkin, O.:  Convergence in the temperature response
#         of leaf respiration across biomes and plant functional types, Proceedings of the National Academy of Sciences,
#         113,  3832–3837,  doi:10.1073/pnas.1520282113,2016.
#
# Huber,  M.  L.,  Perkins,  R.  A.,  Laesecke,  A.,  Friend,  D.  G.,  Sengers,  J.  V.,  Assael,M. J.,
#          Metaxa, I. N., Vogel, E., Mares, R., and Miyagawa, K.:  New international formulation for the viscosity
#          of H2O, Journal of Physical and Chemical ReferenceData, 38, 101–125, 2009
#
# Prentice,  I. C.,  Dong,  N.,  Gleason,  S. M.,  Maire,  V.,  and Wright,  I. J.:  Balancing the costs
#            of carbon gain and water transport:  testing a new theoretical frameworkfor  plant  functional  ecology,
#            Ecology  Letters,  17,  82–91,  10.1111/ele.12211,http://dx.doi.org/10.1111/ele.12211, 2014.
#
# Wang, H., Prentice, I. C., Keenan, T. F., Davis, T. W., Wright, I. J., Cornwell, W. K.,Evans, B. J.,
#          and Peng, C.:  Towards a universal model for carbon dioxide uptake by plants, Nat Plants, 3, 734–741, 2017.
#           Atkin, O. K., et al.:  Global variability in leaf respiration in relation to climate, plant func-tional
#           types and leaf traits, New Phytologist, 206, 614–636, doi:10.1111/nph.13253,
#           https://nph.onlinelibrary.wiley.com/doi/abs/10.1111/nph.13253.
#
# Smith, N. G., Keenan, T. F., Colin Prentice, I. , Wang, H. , Wright, I. J., Niinemets, U. , Crous, K. Y.,
#           Domingues, T. F., Guerrieri, R. , Yoko Ishida, F. , Kattge, J. , Kruger, E. L., Maire, V. , Rogers, A. ,
#           Serbin, S. P., Tarvainen, L. , Togashi, H. F., Townsend, P. A., Wang, M. , Weerasinghe, L. K. and Zhou, S.
#           (2019), Global photosynthetic capacity is optimized to the environment. Ecol Lett, 22: 506-517.
#           doi:10.1111/ele.13210
#
# Stocker, B. et al. Geoscientific Model Development Discussions (in prep.)
#
# @export
#
# @examples rpmodel( tc = 20, vpd = 1000, co2 = 400, fapar = 1, ppfd = 300, elv = 0)