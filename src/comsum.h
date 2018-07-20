c
c ------
c comsum
c ------
c
      integer 
     >  ndtimes,	    ! counter for daily average calculations
     >  nmtimes,            ! counter for monthly average calculations
     >  nytimes             ! counter for yearly average calculations
c
c     common /comsum1/ ndtimes, nmtimes, nytimes
c
c daily average fields
c
      real 
     >  adrain(npoi),       ! daily average rainfall rate (mm/day)
     >  adsnow(npoi),       ! daily average snowfall rate (mm/day)
     >  adaet(npoi),        ! daily average aet (mm/day)
     >  adtrunoff(npoi),    ! daily average total runoff (mm/day)
     >  adsrunoff(npoi),    ! daily average surface runoff (mm/day)
     >  addrainage(npoi),   ! daily average drainage (mm/day)
     >  adrh(npoi),         ! daily average rh (percent)
     >  adsnod(npoi),       ! daily average snow depth (m)
     >  adsnof(npoi),       ! daily average snow fraction (fraction)
     >  adwsoi(npoi),       ! daily average soil moisture (fraction)
     >  adwisoi(npoi),      ! daily average soil ice (fraction)
     >  adtsoi(npoi),       ! daily average soil temperature (c)
     >  adwsoic(npoi),      ! daily average soil moisture using root profile weighting (fraction)
     >  adtsoic(npoi),      ! daily average soil temperature (c) using profile weighting
     >  adco2mic(npoi),     ! daily accumulated co2 respiration from microbes (kg_C m-2 /day)
     >  adco2root(npoi),    ! daily accumulated co2 respiration from roots (kg_C m-2 /day)
     >  adco2soi(npoi),     ! daily accumulated co2 respiration from soil(total) (kg_C m-2 /day)
     >  adco2ratio(npoi),   ! ratio of root to total co2 respiration
     >  adnmintot(npoi),    ! daily accumulated net nitrogen mineralization (kg_N m-2 /day)
     >  adtlaysoi(npoi),    ! daily average soil temperature (c) of top layer
     >  adwlaysoi(npoi),    ! daily average soil moisture of top layer(fraction)
     >     adneetot(npoi)       ! daily accumulated net ecosystem exchange of co2 in ecosystem (kg-C/m**2/day)
c
c     common /comsum3/ adrain, adsnow, adaet, adtrunoff, adsrunoff, addrainage,
     >     adrh, adsnod, adsnof, adwsoi, adwisoi, adtsoi, adwsoic, adtsoic,
     >     adco2mic, adco2root, adco2soi, adco2ratio, adnmintot,
     >     adtlaysoi, adwlaysoi, adneetot
c
c monthly average fields
c
      real
     >  amtemp(npoi),       ! monthly average air temperature (C)
     >  amrain(npoi),       ! monthly average rainfall rate (mm/day)
     >  amsnow(npoi),       ! monthly average snowfall rate (mm/day)
     >  amcloud(npoi),      ! monthly average cloudiness (percent)
     >  amrh(npoi),         ! monthly average rh (percent)
     >  amqa(npoi),         ! monthly average specific humidity (kg-h2o/kg-air)
     >  amaet(npoi),        ! monthly average aet (mm/day)
     >  amtrunoff(npoi),    ! monthly average total runoff (mm/day)
     >  amsrunoff(npoi),    ! monthly average surface runoff (mm/day)
     >  amdrainage(npoi),   ! monthly average drainage (mm/day)
     >  amwsoi(npoi),       ! monthly average 1m soil moisture (fraction)
     >  amwisoi(npoi),      ! monthly average 1m soil ice (fraction)
     >  amvwc(npoi),        ! monthly average 1m volumetric water content (fraction)
     >  amawc(npoi),        ! monthly average 1m plant-available water content (fraction)
     >  amtsoi(npoi),       ! monthly average 1m soil temperature (C)
     >  amsnod(npoi),       ! monthly average snow depth (m)
     >  amsnof(npoi),       ! monthly average snow fraction (fraction)
     >  amlaiu(npoi),       ! monthly average lai for upper canopy (m**2/m**2)
     >  amlail(npoi),       ! monthly average lai for lower canopy (m**2/m**2)
     >  amsolar(npoi),      ! monthly average incident solar radiation (W/m**2)
     >  amalbedo(npoi),     ! monthly average solar albedo (fraction)
     >  amirdown(npoi),     ! monthly average downward ir radiation (W/m**2)
     >  amirup(npoi),       ! monthly average upward ir radiation (W/m**2)
     >  amsens(npoi),       ! monthly average sensible heat flux (W/m**2)
     >  amlatent(npoi),     ! monthly average latent heat flux (W/m**2)
     >  amnpptot(npoi),     ! monthly total npp for ecosystem (kg-C/m**2/month)
     >  amneetot(npoi),	    ! monthly total net ecosystem exchange of CO2 (kg-C/m**2/month)
     >  amco2mic(npoi),     ! monthly total CO2 flux from microbial respiration (kg-C/m**2/month)
     >  amco2root(npoi),    ! monthly total CO2 flux from soil due to root respiration (kg-C/m**2/month)
     >  amco2soi(npoi),     ! monthly total soil CO2 flux from microbial and root respiration (kg-C/m**2/month)
     >  amco2ratio(npoi),   ! monthly ratio of root to total co2 flux
     >  amnmintot(npoi)     ! monthly total N mineralization from microbes (kg-N/m**2/month)
c
c     common /comsum4/ amtemp, amrain, amsnow, amcloud, amrh, amqa, amaet,
     >     amtrunoff, amsrunoff, 
     >     amdrainage, amwsoi, amwisoi, amvwc, amawc, amtsoi, amsnod,
     >     amsnof, amlaiu, amlail, amsolar, amalbedo, amirdown, amirup, 
     >     amsens, amlatent, amnpptot, amneetot, amco2mic, amco2root, 
     >     amco2soi, amco2ratio, amnmintot
c
      real
     >  amnpp(npoi,npft)    ! monthly total npp for each plant type (kg-C/m**2/month)
c
c     common /comsum5/ amnpp
c
c annual average fields
c
      real
     >  ayprcp(npoi),       ! annual average precipitation (mm/yr)
     >  ayaet(npoi),        ! annual average aet (mm/yr)
     >  aytrans(npoi),      ! annual average transpiration (mm/yr)
     >  aytrunoff(npoi),    ! annual average total runoff (mm/yr)
     >  aysrunoff(npoi),    ! annual average surface runoff (mm/yr)
     >  aydrainage(npoi),   ! annual average drainage (mm/yr)
     >  aydwtot(npoi),      ! annual average soil+vegetation+snow water recharge (mm/yr or kg_h2o/m**2/yr)
     >  aywsoi(npoi),       ! annual average 1m soil moisture (fraction)
     >  aywisoi(npoi),      ! annual average 1m soil ice (fraction)
     >  ayvwc(npoi),        ! annual average 1m volumetric water content (fraction)
     >  ayawc(npoi),        ! annual average 1m plant-available water content (fraction)
     >  aytsoi(npoi),       ! annual average 1m soil temperature (C)
     >  ayrratio(npoi),     ! annual average runoff ratio (fraction)
     >  aytratio(npoi),     ! annual average transpiration ratio (fraction)
     >  aysolar(npoi),      ! annual average incident solar radiation (w/m**2)
     >  ayalbedo(npoi),     ! annual average solar albedo (fraction)
     >  ayirdown(npoi),     ! annual average downward ir radiation (w/m**2)
     >  ayirup(npoi),       ! annual average upward ir radiation (w/m**2)
     >  aysens(npoi),       ! annual average sensible heat flux (w/m**2)
     >  aylatent(npoi),     ! annual average latent heat flux (w/m**2)
     >  aystresstu(npoi),   ! annual average soil moisture stress parameter for upper canopy (dimensionless)
     >  aystresstl(npoi),   ! annual average soil moisture stress parameter for lower canopy (dimensionless)
     >  ayanpptot(npoi),    ! annual above-ground npp for ecosystem (kg-c/m**2/yr)
     >  aynpptot(npoi),     ! annual total npp for ecosystem (kg-c/m**2/yr)
     >  aygpptot(npoi),     ! annual total gpp for ecosystem (kg-c/m**2/yr)
     >  ayalit(npoi),       ! aboveground litter (kg-c/m**2)
     >  ayblit(npoi),       ! belowground litter (kg-c/m**2)
     >  aycsoi(npoi),       ! total soil carbon (kg-c/m**2)
     >  aycmic(npoi),       ! total soil carbon in microbial biomass (kg-c/m**2)
     >  ayanlit(npoi),      ! aboveground litter nitrogen (kg-N/m**2)
     >  aybnlit(npoi),      ! belowground litter nitrogen (kg-N/m**2)
     >  aynsoi(npoi),       ! total soil nitrogen (kg-N/m**2)
     >  ynleach(npoi),
     >  ayneetot(npoi),     ! annual total NEE for ecosystem (kg-C/m**2/yr)
     >  ayco2mic(npoi),     ! annual total CO2 flux from microbial respiration (kg-C/m**2/yr)
     >  ayco2root(npoi),    ! annual total CO2 flux from soil due to root respiration (kg-C/m**2/yr)
     >  ayco2soi(npoi),     ! annual total soil CO2 flux from microbial and root respiration (kg-C/m**2/yr)
     >  aynmintot(npoi),    ! annual total nitrogen mineralization (kg-N/m**2/yr)
     >  ayrootbio(npoi)     ! annual average live root biomass (kg-C / m**2)
c
c     common /comsum6/ ayprcp, ayaet, aytrans, aytrunoff, aysrunoff, aydrainage, aydwtot,
     >     aywsoi, aywisoi, ayvwc, ayawc, aytsoi, ayrratio, aytratio,
     >     aysolar, ayalbedo, ayirdown, ayirup, aysens, aylatent,
     >     aystresstu, aystresstl, ayanpptot, aynpptot, aygpptot, 
     >     ayalit, ayblit, aycsoi, aycmic, ayanlit, aybnlit, aynsoi,ynleach,
     >     ayneetot, ayco2mic, ayco2root, ayco2soi, aynmintot, ayrootbio
c
      real 
     >  ayanpp(npoi,npft),  ! annual above-ground npp for each plant type(kg-c/m**2/yr)
     >  aynpp(npoi,npft),   ! annual total npp for each plant type(kg-c/m**2/yr)
     >  aygpp(npoi,npft)    ! annual gross npp for each plant type(kg-c/m**2/yr)
c
c     common /comsum7/ ayanpp, aynpp, aygpp
c
c other time average fields
c
      real 
     >  a10td(npoi),        ! 10-day average daily air temperature (K)
     >  a10ancub(npoi),     ! 10-day average canopy photosynthesis rate - broadleaf (mol_co2 m-2 s-1)
     >  a10ancuc(npoi),     ! 10-day average canopy photosynthesis rate - conifer (mol_co2 m-2 s-1)
     >  a10ancls(npoi),     ! 10-day average canopy photosynthesis rate - shrubs (mol_co2 m-2 s-1)
     >  a10ancl3(npoi),     ! 10-day average canopy photosynthesis rate - c3 grasses (mol_co2 m-2 s-1)
     >  a10ancl4(npoi),     ! 10-day average canopy photosynthesis rate - c4 grasses (mol_co2 m-2 s-1)
     >  a10scalparamu(npoi),    ! 10-day average day-time scaling parameter - upper canopy (dimensionless)
     >  a10scalparaml(npoi),    ! 10-day average day-time scaling parameter - lower canopy (dimensionless)
     >  a10daylightu(npoi),    ! 10-day average day-time PAR - upper canopy (micro-Ein m-2 s-1)
     >  a10daylightl(npoi)     ! 10-day average day-time PAR - lower canopy (micro-Ein m-2 s-1)
c
c     common /comsum8/ a10td, a10ancub, a10ancuc, a10ancls, a10ancl3, a10ancl4, 
     >     a10scalparamu, a10scalparaml, a10daylightu, a10daylightl
c
c biogeochem summations
c
      real 
     >  storedn(npoi),      ! total storage of N in soil profile (kg_N m-2) 
     >  yrleach(npoi)       ! annual total amount C leached from soil profile (kg_C m-2/yr)
c
c     common /comsum9/ storedn, yrleach
c
