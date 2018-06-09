c
c ------
c comveg 
c ------
c
      real 
     >  alaiml,             ! lower canopy leaf & stem maximum area (2 sided) for normalization of drag coefficient (m2 m-2)
     >  alaimu,             ! upper canopy leaf & stem area (2 sided) for normalization of drag coefficient (m2 m-2)
     >  cgrass,             ! empirical constant in lower canopy-air aerodynamic transfer coefficient (m s-0.5) (A39a Pollard & Thompson 95)
     >  chl,                ! heat capacity of lower canopy leaves & stems per unit leaf/stem area (J kg-1 m-2)
     >  chs,                ! heat capacity of upper canopy stems per unit stem area (J kg-1 m-2)
     >  chu,                ! heat capacity of upper canopy leaves per unit leaf area (J kg-1 m-2)
     >  cleaf,              ! empirical constant in upper canopy leaf-air aerodynamic transfer coefficient (m s-0.5) (A39a Pollard & Thompson 95)
     >  cstem,              ! empirical constant in upper canopy stem-air aerodynamic transfer coefficient (m s-0.5) (A39a Pollard & Thompson 95)
     >  tblowl,             ! decay time for blowoff of snow intercepted by lower canopy leaves & stems (sec)
     >  tblows,             ! decay time for blowoff of snow intercepted by upper canopy stems (sec)
     >  tblowu,             ! decay time for blowoff of snow intercepted by upper canopy leaves (sec)
     >  tdripl,             ! decay time for dripoff of liquid intercepted by lower canopy leaves & stem (sec)
     >  tdrips,             ! decay time for dripoff of liquid intercepted by upper canopy stems (sec) 
     >  tdripu,             ! decay time for dripoff of liquid intercepted by upper canopy leaves (sec)
     >  wliqmin,            ! minimum intercepted water on unit vegetated area (kg m-2)
     >  wliqlmax,           ! maximum intercepted water on a unit lower canopy stem & leaf area (kg m-2)
     >  wliqsmax,           ! maximum intercepted water on a unit upper canopy stem area (kg m-2)
     >  wliqumax,           ! maximum intercepted water on a unit upper canopy leaf area (kg m-2)
     >  wsnomin,            ! minimum intercepted snow on unit vegetated area (kg m-2)
     >  wsnolmax,           ! intercepted snow capacity for lower canopy leaves & stems (kg m-2)
     >  wsnosmax,           ! intercepted snow capacity for upper canopy stems (kg m-2)
     >  wsnoumax,           ! intercepted snow capacity for upper canopy leaves (kg m-2)
     >  woodnorm	    ! value of woody biomass for upper canopy closure (ie when wood = woodnorm fu = 1.0) (kg_C m-2)
c
      common /comveg1/ alaiml, alaimu, cgrass, chl, chs, chu, cleaf, cstem, tblowl,
     >     tblows, tblowu, tdripl, tdrips, tdripu, wliqmin, wliqlmax,
     >     wliqsmax, wliqumax, wsnomin, wsnolmax, wsnosmax, wsnoumax, 
     >     woodnorm
c
      real
     >  q12(npoi),          ! specific humidity of air at z12
     >  q34(npoi),          ! specific humidity of air at z34
     >  sl(npoi),           ! air-vegetation transfer coefficients (*rhoa) for lower canopy leaves & stems (m s-1*kg m-3) (A39a Pollard & Thompson 1995)
     >  ss(npoi),           ! air-vegetation transfer coefficients (*rhoa) for upper canopy stems (m s-1 * kg m-3) (A39a Pollard & Thompson 1995)
     >  su(npoi),           ! air-vegetation transfer coefficients (*rhoa) for upper canopy leaves (m s-1 * kg m-3) (A39a Pollard & Thompson 1995)
     >  topparl(npoi),      ! total photosynthetically active raditaion absorbed by top leaves of lower canopy (W m-2)
     >  topparu(npoi),      ! total photosynthetically active raditaion absorbed by top leaves of upper canopy (W m-2)
     >  tl(npoi),           ! temperature of lower canopy leaves & stems(K)
     >  ts(npoi),           ! temperature of upper canopy stems (K)
     >  tu(npoi),           ! temperature of upper canopy leaves (K)
     >  tlsub(npoi),        ! temperature of lower canopy vegetation buried by snow (K)
     >  t12(npoi),          ! air temperature at z12 (K)
     >  t34(npoi),          ! air temperature at z34 (K)
     >  wliql(npoi),        ! intercepted liquid h2o on lower canopy leaf and stem area (kg m-2)
     >  wliqs(npoi),        ! intercepted liquid h2o on upper canopy stem area (kg m-2)
     >  wliqu(npoi),        ! intercepted liquid h2o on upper canopy leaf area (kg m-2)
     >  wsnol(npoi),        ! intercepted frozen h2o (snow) on lower canopy leaf & stem area (kg m-2)
     >  wsnos(npoi),        ! intercepted frozen h2o (snow) on upper canopy stem area (kg m-2)
     >  wsnou(npoi)         ! intercepted frozen h2o (snow) on upper canopy leaf area (kg m-2)
c
      common /comveg2/ q12, q34, sl, ss, su, topparl, topparu, tl, ts, tu, tlsub, 
     >     t12, t34, wliql, wliqs, wliqu, wsnol, wsnos, wsnou
c
c all photosynthesis rates are per unit leaf area
c
      real
     >  agcub(npoi),        ! canopy average gross photosynthesis rate - broadleaf  (mol_co2 m-2 s-1)
     >  agcuc(npoi),        ! canopy average gross photosynthesis rate - conifer    (mol_co2 m-2 s-1)
     >  agcls(npoi),        ! canopy average gross photosynthesis rate - shrubs     (mol_co2 m-2 s-1)
     >  agcl3(npoi),        ! canopy average gross photosynthesis rate - c3 grasses (mol_co2 m-2 s-1)
     >  agcl4(npoi),        ! canopy average gross photosynthesis rate - c4 grasses (mol_co2 m-2 s-1)
     >  ancub(npoi),        ! canopy average net photosynthesis rate - broadleaf    (mol_co2 m-2 s-1)
     >  ancuc(npoi),        ! canopy average net photosynthesis rate - conifer      (mol_co2 m-2 s-1)
     >  ancls(npoi),        ! canopy average net photosynthesis rate - shrubs       (mol_co2 m-2 s-1)
     >  ancl3(npoi),        ! canopy average net photosynthesis rate - c3 grasses   (mol_co2 m-2 s-1)
     >  ancl4(npoi),        ! canopy average net photosynthesis rate - c4 grasses   (mol_co2 m-2 s-1)
     >  totcondub(npoi),     ! 
     >  totconduc(npoi),     !
     >  totcondls(npoi),     ! 
     >  totcondl3(npoi),     !
     >  totcondl4(npoi)      !
c
      common /comveg3/ agcub, agcuc, agcls, agcl3, agcl4, 
     >     ancub, ancuc, ancls, ancl3, ancl4, 
     >     totcondub, totconduc, totcondls, totcondl3, totcondl4
c
      real
     >  ciub(npoi),         ! intercellular co2 concentration - broadleaf (mol_co2/mol_air)
     >  ciuc(npoi),         ! intercellular co2 concentration - conifer   (mol_co2/mol_air)
     >  cils(npoi),         ! intercellular co2 concentration - shrubs    (mol_co2/mol_air)
     >  cil3(npoi),         ! intercellular co2 concentration - c3 plants (mol_co2/mol_air)
     >  cil4(npoi),         ! intercellular co2 concentration - c4 plants (mol_co2/mol_air)
     >  csub(npoi),         ! leaf boundary layer co2 concentration - broadleaf (mol_co2/mol_air)
     >  csuc(npoi),         ! leaf boundary layer co2 concentration - conifer   (mol_co2/mol_air)
     >  csls(npoi),         ! leaf boundary layer co2 concentration - shrubs    (mol_co2/mol_air)
     >  csl3(npoi),         ! leaf boundary layer co2 concentration - c3 plants (mol_co2/mol_air)
     >  csl4(npoi),         ! leaf boundary layer co2 concentration - c4 plants (mol_co2/mol_air)
     >  gsub(npoi),         ! upper canopy stomatal conductance - broadleaf  (mol_co2 m-2 s-1)
     >  gsuc(npoi),         ! upper canopy stomatal conductance - conifer    (mol_co2 m-2 s-1)
     >  gsls(npoi),         ! lower canopy stomatal conductance - shrubs     (mol_co2 m-2 s-1)
     >  gsl3(npoi),         ! lower canopy stomatal conductance - c3 grasses (mol_co2 m-2 s-1)
     >  gsl4(npoi)          ! lower canopy stomatal conductance - c4 grasses (mol_co2 m-2 s-1)
c
      common /comveg4/ ciub, ciuc, cils, cil3, cil4, 
     >     csub, csuc, csls, csl3, csl4,
     >     gsub, gsuc, gsls, gsl3, gsl4
c
      real 
     >  agddl(npoi),        ! annual accumulated growing degree days for bud burst, lower canopy (day-degrees)
     >  agddu(npoi),        ! annual accumulated growing degree days for bud burst, upper canopy (day-degrees)
     >  fl(npoi),           ! fraction of snow-free area covered by lower  canopy
     >  fu(npoi),           ! fraction of overall area covered by upper canopy
     >  gdd0(npoi),         ! growing degree days > 0C 
     >  gdd0this(npoi),     ! annual total growing degree days for current year
     >  gdd5(npoi),         ! growing degree days > 5C
     >  gdd5this(npoi),     ! annual total growing degree days for current year
     >  sapfrac(npoi),      ! fraction of woody biomass that is in sapwood
     >  tc(npoi),           ! coldest monthly temperature (C)
     >  tcthis(npoi),       ! coldest monthly temperature of current year (C)
     >  tcmin(npoi),        ! coldest daily temperature of current year (C)
     >  totlail(npoi),      ! total leaf area index for the lower canopy
     >  totlaiu(npoi),      ! total leaf area index for the upper canopy
     >  totbiol(npoi),      ! total biomass in the lower canopy (kg_C m-2)
     >  totbiou(npoi),      ! total biomass in the upper canopy (kg_C m-2)
     >  tw(npoi),           ! warmest monthly temperature (C)
     >  twthis(npoi),       ! warmest monthly temperature of current year (C)
     >  disturbf(npoi),     ! annual fire disturbance regime (m2/m2/yr)
     >  disturbo(npoi),     ! fraction of biomass pool lost every year to disturbances other than fire
     >  firefac(npoi),      ! factor that respresents the annual average fuel dryness of a grid cell, and hence characterizes the readiness to burn
     >  tco2mic(npoi),      ! instantaneous microbial co2 flux from soil (mol-CO2 / m-2 / second)
     >  tco2root(npoi),     ! instantaneous fine co2 flux from soil (mol-CO2 / m-2 / second)
     >  tneetot(npoi),      ! instantaneous net ecosystem exchange of co2 per timestep (kg_C m-2/timestep)
     >  tnmin(npoi),        ! instantaneous nitrogen mineralization (kg_N m-2/timestep)
     >  tnpptot(npoi),      ! instantaneous npp (mol-CO2 / m-2 / second)
     >  tgpptot(npoi),      ! instantaneous gpp (mol-CO2 / m-2 / second)
     >  totalit(npoi),	    ! total standing aboveground litter (kg_C m-2)
     >  totanlit(npoi),	    ! total standing aboveground nitrogen in litter (kg_N m-2)
     >  totcmic(npoi),      ! total carbon residing in microbial pools (kg_C m-2)
     >  totcsoi(npoi),      ! total carbon in all soil pools (kg_C m-2)
     >  totfall(npoi),	    ! total litterfall and root turnover (kg_C m-2/year)
     >  totlit(npoi),       ! total carbon in all litter pools (kg_C m-2)
     >  totnlit(npoi),      ! total nitrogen in all litter pools (kg_N m-2)
     >  totnmic(npoi),      ! total nitrogen residing in microbial pool (kg_N m-2)
     >  totnsoi(npoi),      ! total nitrogen in soil (kg_N m-2)
     >  totrlit(npoi),      ! total root litter carbon belowground (kg_C m-2)
     >  totrnlit(npoi),     ! total root litter nitrogen belowground (kg_N m-2)
     >  tempu(npoi),        ! cold-phenology trigger for trees (non-dimensional)
     >  templ(npoi),        ! cold-phenology trigger for grasses/shrubs (non-dimensional)
     >  dropu(npoi),        ! drought-phenology trigger for trees (non-dimensional)
     >  dropls(npoi),       ! drought-phenology trigger for shrubs (non-dimensional)
     >  dropl4(npoi),       ! drought-phenology trigger for c4 grasses (non-dimensional)
     >  dropl3(npoi),       ! drought-phenology trigger for c3 grasses (non-dimensional)
     >  vegtype0(npoi)      ! annual vegetation type - ibis classification
c
      common /comveg5/ agddl, agddu, fl, fu, 
     >     gdd0, gdd0this, gdd5, gdd5this, 
     >     sapfrac, tc, tcthis, tcmin, 
     >     totlail, totlaiu, totbiol, totbiou, 
     >     tw, twthis, 
     >     disturbf, disturbo, firefac, 
     >     tco2mic, tco2root, tneetot, tnmin, tnpptot, tgpptot, 
     >     totalit, totanlit, totcmic, totcsoi, totfall, totlit, 
     >     totnlit, totnmic, totnsoi, totrlit, totrnlit, 
     >     tempu, templ, dropu, dropls, dropl4, dropl3, vegtype0
c
      real 
     >  cbiol(npoi,npft),   ! carbon in leaf biomass pool (kg_C m-2)
     >  cbior(npoi,npft),   ! carbon in fine root biomass pool (kg_C m-2)
     >  cbiow(npoi,npft)    ! carbon in woody biomass pool (kg_C m-2)
c
      common /comveg6/ cbiol, cbior, cbiow
c
      real 
     >  clitll(npoi),       ! carbon in leaf litter pool - lignin          (kg_C m-2)
     >  clitlm(npoi),       ! carbon in leaf litter pool - metabolic       (kg_C m-2)
     >  clitls(npoi),       ! carbon in leaf litter pool - structural      (kg_C m-2)
     >  clitrl(npoi),       ! carbon in fine root litter pool - lignin     (kg_C m-2)
     >  clitrm(npoi),       ! carbon in fine root litter pool - metabolic  (kg_C m-2)
     >  clitrs(npoi),       ! carbon in fine root litter pool - structural (kg_C m-2)
     >  clitwl(npoi),       ! carbon in woody litter pool - lignin         (kg_C m-2)
     >  clitwm(npoi),       ! carbon in woody litter pool - metabolic      (kg_C m-2)
     >  clitws(npoi),       ! carbon in woody litter pool - structural     (kg_C m-2)
     >  csoipas(npoi),      ! carbon in soil - passive humus               (kg_C m-2)
     >  csoislo(npoi),      ! carbon in soil - slow humus                  (kg_C m-2)
     >  csoislon(npoi),     ! carbon in soil - slow nonprotected humus     (kg_C m-2)
     >  csoislop(npoi),     ! carbon in soil - slow protected humus        (kg_C m-2)
     >  decompl(npoi),      ! litter decomposition factor                  (dimensionless)
     >  decomps(npoi),      ! soil organic matter decomposition factor     (dimensionless)
     >  falll(npoi),        ! annual leaf litter fall                      (kg_C m-2/year)
     >  fallr(npoi),        ! annual root litter input                     (kg_C m-2/year)
     >  fallw(npoi),        ! annual wood litter fall                      (kg_C m-2/year)
     >  cdisturb(npoi)      ! annual amount of vegetation carbon lost 
                            ! to atmosphere due to fire  (biomass burning) (kg_C m-2/year)
c
      common /comveg7/ clitll, clitlm, clitls, clitrl, clitrm, clitrs, clitwl,
     >     clitwm, clitws, csoipas, csoislo, csoislon, csoislop, decompl
     >     , decomps, falll, fallr, fallw, cdisturb
c
      real 
     >  biomass(npoi,npft), ! total biomass of each plant functional type  (kg_C m-2)
     >  frac(npoi,npft),    ! fraction of canopy occupied by each plant functional type
     >  plai(npoi,npft),    ! total leaf area index of each plant functional type
     >  tnpp(npoi,npft),    ! instantaneous NPP for each pft (mol-CO2 / m-2 / second)
     >  nppdummy(npoi,npft), ! canopy NPP before accounting for stem and root respiration
     >  tgpp(npoi,npft)     ! instantaneous GPP for each pft (mol-CO2 / m-2 / second)
c
      common /comveg8/ biomass, frac, plai, tnpp, tgpp
c
c 1: lower canopy
c 2: upper canopy
c
      real
     >  lai(npoi,2),        ! canopy single-sided leaf area index (area leaf/area veg)
     >  sai(npoi,2),        ! current single-sided stem area index
     >  zbot(npoi,2),       ! height of lowest branches above ground (m)
     >  ztop(npoi,2),       ! height of plant top above ground (m)
     >  ztopmx(npoi,2)      ! maximum annual height of plant top above ground (m) 
c
      common /comveg9/ lai, sai, zbot, ztop, ztopmx
c
      real
     >  dleaf(2),           ! typical linear leaf dimension in aerodynamic transfer coefficient (m)
     >  dstem(2),           ! typical linear stem dimension in aerodynamic transfer coefficient (m)
     >  orieh(2),           ! fraction of leaf/stems with horizontal orientation
     >  oriev(2)            ! fraction of leaf/stems with vertical
c
      common /comveg10/  dleaf, dstem, orieh, oriev
c
      real 
     >  rhoveg(nband,2),    ! reflectance of an average leaf/stem
     >  tauveg(nband,2)     ! transmittance of an average leaf/stem
c
      common /comveg11/ rhoveg, tauveg
c
      real 
     >  froot(nsoilay,2)    ! fraction of root in soil layer 
c
      common /comveg12/ froot
c
      real 
     >  exist(npoi,npft)    ! probability of existence of each plant functional type in a gridcell
c
      common /comveg13/ exist
c
      real 
     >  specla(npft),        ! specific leaf area (m**2/kg) 
     >  aleaf(npft),         ! carbon allocation fraction to leaves
     >  aroot(npft),         ! carbon allocation fraction to fine roots
     >  awood(npft)          ! carbon allocation fraction to wood 
c
      common /comveg14/ specla, aleaf, aroot, awood
c

