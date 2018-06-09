c
c  ####    #####    ##     #####   ####
c #          #     #  #      #    #
c  ####      #    #    #     #     ####
c      #     #    ######     #         #
c #    #     #    #    #     #    #    #
c  ####      #    #    #     #     ####
c
c
c ---------------------------------------------------------------------
       subroutine sumnow
c ---------------------------------------------------------------------
c
c common blocks
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'comatm.h'
      include 'comhyd.h'
      include 'comsno.h'
      include 'comsoi.h'
      include 'comveg.h'
      include 'comsum.h'
c
c local variables
c
      integer i,k         ! loop indices
c
      real rwood,    ! maintenance respiration coefficient for wood (/s)
     >     rroot,    ! maintenance respiration coefficient for root (/s)
     >     rgrowth,  ! growth respiration coefficient (fraction)
     >     stemtemp, ! stem temperature
     >     roottemp, ! average root temperature for all roots
     >     funca,    ! temperature function for aboveground biomass (stems)
     >     funcb,    ! temperature function for belowground biomass (roots)
     >     zweight,  ! 10-day time averaging factor
     >     smask     ! 1 - fi
c
c
c ---------------------------------------------------------------------
c * * * define working variables * * *
c ---------------------------------------------------------------------
c
c maintenance respiration coefficients (per second)
c
c initially, we pick values for respiration coefficients that
c defined in units of  / year
c
c   rwood ~ 0.0125 
c   rroot ~ 1.2500
c
c however, we convert the unitsconvert to have resulting respiration
c fluxes in units of mol-C / m**2 / second
c
c this requires we convert the time unit to seconds and add an additional
c factor to convert biomass units from kilograms to moles
c
      rwood   = 0.0125 / (ndaypy * 86400.0) * (1000.0 / 12.0)
      rroot   = 1.2500 / (ndaypy * 86400.0) * (1000.0 / 12.0)
c
c growth respiration coefficient (fraction)
c
      rgrowth = 0.30
c
c 10-day time averaging factor
c
      zweight = exp(-1. / (10.0 * 86400.0 / dtime)) 
c
c begin global grid
c
      do 100 i = 1, npoi
c
c calculate instantaneous carbon flux parameters, including
c npp (net primary production) and nee (net ecosystem exchange)
c
c in this routine, all of the fluxes are calculated in the units
c of mol-C / m**2 / sec
c
c ---------------------------------------------------------------------
c * * * calculate instantaneous GPP * * *
c ---------------------------------------------------------------------
c
c snow masking for lower canopy vegetation
c
        smask = 1.0 - fi(i)
c
c note that the following plants types follow different physiological paths
c
c   - broadleaf trees   :  types 1, 2, 3, 5, 7, 8 
c   - conifer   trees   :  types 4, 6
c   - shrubs            :  types 9, 10
c   - c4 grasses        :  type 11
c   - c3 grasses        :  type 12
c
c note that plant type 8 is actually a deciduous conifer (e.g., Larix), but
c we are assuming that it's physiological behavior is like a broadleaf tree
c
c nppdummy is canopy npp before accounting for stem & root respirtation
c Navin Sept 02
c
        nppdummy(i,1)  = frac(i,1)  * ancub(i) * lai(i,2) * fu(i)
        nppdummy(i,2)  = frac(i,2)  * ancub(i) * lai(i,2) * fu(i)
        nppdummy(i,3)  = frac(i,3)  * ancub(i) * lai(i,2) * fu(i)
        nppdummy(i,4)  = frac(i,4)  * ancuc(i) * lai(i,2) * fu(i)
        nppdummy(i,5)  = frac(i,5)  * ancub(i) * lai(i,2) * fu(i)
        nppdummy(i,6)  = frac(i,6)  * ancuc(i) * lai(i,2) * fu(i)
        nppdummy(i,7)  = frac(i,7)  * ancub(i) * lai(i,2) * fu(i)
        nppdummy(i,8)  = frac(i,8)  * ancub(i) * lai(i,2) * fu(i)
        nppdummy(i,9)  = frac(i,9)  * ancls(i) * lai(i,1) * fl(i) * smask 
        nppdummy(i,10) = frac(i,10) * ancls(i) * lai(i,1) * fl(i) * smask
        nppdummy(i,11) = frac(i,11) * ancl4(i) * lai(i,1) * fl(i) * smask
        nppdummy(i,12) = frac(i,12) * ancl3(i) * lai(i,1) * fl(i) * smask
c
c Navin's correction to compute npp using tgpp via agXXX
c agXXX should be used 
c
        tgpp(i,1)  = frac(i,1)  * agcub(i) * lai(i,2) * fu(i)
        tgpp(i,2)  = frac(i,2)  * agcub(i) * lai(i,2) * fu(i)
        tgpp(i,3)  = frac(i,3)  * agcub(i) * lai(i,2) * fu(i)
        tgpp(i,4)  = frac(i,4)  * agcuc(i) * lai(i,2) * fu(i)
        tgpp(i,5)  = frac(i,5)  * agcub(i) * lai(i,2) * fu(i)
        tgpp(i,6)  = frac(i,6)  * agcuc(i) * lai(i,2) * fu(i)
        tgpp(i,7)  = frac(i,7)  * agcub(i) * lai(i,2) * fu(i)
        tgpp(i,8)  = frac(i,8)  * agcub(i) * lai(i,2) * fu(i)
        tgpp(i,9)  = frac(i,9)  * agcls(i) * lai(i,1) * fl(i) * smask 
        tgpp(i,10) = frac(i,10) * agcls(i) * lai(i,1) * fl(i) * smask
        tgpp(i,11) = frac(i,11) * agcl4(i) * lai(i,1) * fl(i) * smask
        tgpp(i,12) = frac(i,12) * agcl3(i) * lai(i,1) * fl(i) * smask
c
c calculate total gridcell gpp
c
        tgpptot(i) = 0.0
c
        do 110 k = 1, npft
          tgpptot(i) = tgpptot(i) + tgpp(i,k)
 110    continue
c
c ---------------------------------------------------------------------
c * * * calculate temperature functions for respiration * * *
c ---------------------------------------------------------------------
c
c calculate the stem temperature
c
        stemtemp = ts(i)
c
c calculate average root temperature (average of all roots)
c
        roottemp = 0.0
c
        do 120 k = 1, nsoilay
          roottemp = roottemp + tsoi(i,k) * 0.5 *
     >               (froot(k,1) + froot(k,2))
 120    continue
c
c calculate respiration terms on a 15 degree base
c following respiration parameterization of Lloyd and Taylor
c
        funca = exp(3500.0 * (1. / 288.16 - 1. / stemtemp))
        funcb = exp(3500.0 * (1. / 288.16 - 1. / roottemp))
c
c ---------------------------------------------------------------------
c * * * calculate instantaneous NPP * * *
c ---------------------------------------------------------------------
c
c the basic equation for npp is
c
c   npp = (1 - growth respiration term) * (gpp - maintenance respiration terms)
c
c here the respiration terms are simulated as
c
c   growth respiration = rgrowth * (gpp - maintenance respiration terms)
c
c where
c
c   rgrowth is the construction cost of new tissues
c
c and
c
c   root respiration = rroot * cbior(i,k) * funcb
c   wood respiration = rwood * cbiow(i,k) * funca * sapwood fraction
c
c where
c 
c   funca = temperature function for aboveground biomass (stems)
c   funcb = temperature function for belowground biomass (roots)
c
c note that we assume the sapwood fraction for shrubs is 1.0
c
c also note that we apply growth respiration, (1 - rgrowth), 
c throughout the year; this may cause problems when comparing
c these npp values with flux tower measurements
c
c also note that we need to convert the mass units of wood and
c root biomass from kilograms of carbon to moles of carbon
c to maintain consistent units (done in rwood, rroot)
c
c finally, note that growth respiration is only applied to 
c positive carbon gains (i.e., when gpp-rmaint is positive)
c
c Navin fix Sept 02 using nppdummy
        tnpp(i,1)  = nppdummy(i,1)                               -
     >               rwood * cbiow(i,1) * sapfrac(i) * funca -
     >               rroot * cbior(i,1)              * funcb
c
        tnpp(i,2)  = nppdummy(i,2)                               -
     >               rwood * cbiow(i,2) * sapfrac(i) * funca -
     >               rroot * cbior(i,2)              * funcb
c
        tnpp(i,3)  = nppdummy(i,3)                               -
     >               rwood * cbiow(i,3) * sapfrac(i) * funca -
     >               rroot * cbior(i,3)              * funcb
c
        tnpp(i,4)  = nppdummy(i,4)                               -
     >               rwood * cbiow(i,4) * sapfrac(i) * funca -
     >               rroot * cbior(i,4)              * funcb
c
        tnpp(i,5)  = nppdummy(i,5)                               -
     >               rwood * cbiow(i,5) * sapfrac(i) * funca -
     >               rroot * cbior(i,5)              * funcb
c
        tnpp(i,6)  = nppdummy(i,6)                               -
     >               rwood * cbiow(i,6) * sapfrac(i) * funca -
     >               rroot * cbior(i,6)              * funcb
c
        tnpp(i,7)  = nppdummy(i,7)                               -
     >               rwood * cbiow(i,7) * sapfrac(i) * funca -
     >               rroot * cbior(i,7)              * funcb
c
        tnpp(i,8)  = nppdummy(i,8)                               -
     >               rwood * cbiow(i,8) * sapfrac(i) * funca -
     >               rroot * cbior(i,8)              * funcb
c
        tnpp(i,9)  = nppdummy(i,9)                               -
     >               rwood * cbiow(i,9)              * funca -
     >               rroot * cbior(i,9)              * funcb
c
        tnpp(i,10) = nppdummy(i,10)                              -
     >               rwood * cbiow(i,10)             * funca -
     >               rroot * cbior(i,10)             * funcb
c
        tnpp(i,11) = nppdummy(i,11)                              -
     >               rroot * cbior(i,11)            * funcb
c
        tnpp(i,12) = nppdummy(i,12)                              -
     >               rroot * cbior(i,12)            * funcb
c
c apply growth respiration and calculate total gridcell npp
c
        tnpptot(i) = 0.0
c
        do 130 k = 1, npft
          if (tnpp(i,k).gt.0.0) tnpp(i,k) = tnpp(i,k)  * (1.0 - rgrowth)
          tnpptot(i) = tnpptot(i) + tnpp(i,k)
 130    continue
c
c ---------------------------------------------------------------------
c * * * calculate total fine root respiration * * *
c ---------------------------------------------------------------------
c
        tco2root(i) = 0.0
c
        do 140 k = 1, npft
          tco2root(i) = tco2root(i) + rroot * cbior(i,k) * funcb
 140    continue
c
c ---------------------------------------------------------------------
c * * * calculate instantaneous NEE * * *
c ---------------------------------------------------------------------
c
c microbial respiration is calculated in biogeochem.f
c
        tneetot(i) = tnpptot(i) - tco2mic(i)
c
c ---------------------------------------------------------------------
c * * * update 10-day running-mean parameters * * *
c ---------------------------------------------------------------------
c
c 10-day daily air temperature
c
        a10td(i)    = zweight * a10td(i)    + (1. - zweight) * td(i)
c
c 10-day canopy photosynthesis rates
c
        a10ancub(i) = zweight * a10ancub(i) + (1. - zweight) * ancub(i)
        a10ancuc(i) = zweight * a10ancuc(i) + (1. - zweight) * ancuc(i)
        a10ancls(i) = zweight * a10ancls(i) + (1. - zweight) * ancls(i)
        a10ancl3(i) = zweight * a10ancl3(i) + (1. - zweight) * ancl3(i)
        a10ancl4(i) = zweight * a10ancl4(i) + (1. - zweight) * ancl4(i)
c
 100  continue
c
c return to main program
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine sumday (istep, plens)
c ---------------------------------------------------------------------
c
c common blocks
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'comatm.h'
      include 'com1d.h'
      include 'comhyd.h'
      include 'comsno.h'
      include 'comsoi.h'
      include 'comsum.h'
      include 'comveg.h'
c
c Arguments
c
      integer
     >  istep   ! daily timestep number (passed in)
c
c local variables
c
      integer i,k      ! loop indices
c
      real 
     >  rwork,      !working time variable
     >  rwork2,     ! "
     >  rwork3,     ! " 
     >  rwork4,     ! "
     >  tconst,     ! constant for Lloyd and Taylor (1994) function
     >  bconst,     ! base temperature used for carbon decomposition
     >  btemp,      ! maximum value of decomposition factor
     >  rdepth,     ! total depth of the 4 1st soil layers
     >  rdepth2,    ! total depth of the 2 1st soil layers
     >  snodpth,    ! total snow depth
     >  soiltemp,   ! average soil temp for 2 1st layers
     >  soilmois,   ! average soil moisture (fraction of porosity) for 2 1st layers
     >  soilice,    ! average soil ice for 2 1st layers
     >  soitempc,   ! average soil temp over 6 layers
     >  soimoisc,   ! average soil moisture over 6 layers
     >  factor,     ! temperature decomposition factor for ltter/soil carbon
     >  wfps,       ! water filled pore space
     >  moist,      ! moisture effect on decomposition
     >  precipfac
c
      real plens
      integer niter
c
c ---------------------------------------------------------------------
c * * * update counters and working variables * * *
c ---------------------------------------------------------------------
c
c reset sumday if the first timestep of the day 
c
      if (istep .eq. 1) ndtimes = 0
c
c accumulate daily output (at this point for soil decomposition)
c
      ndtimes = ndtimes + 1
c
c working variables
c
      rwork  = 1. / float(ndtimes)
      rwork2 = 86400.
      rwork3 = 86400. * 12.e-3
      rwork4 = 86400. * 14.e-3
c
c constants used in temperature function for c decomposition
c (arrhenius function constant) 
c
      tconst  = 344.00  ! constant for Lloyd and Taylor (1994) function
      btemp   = 288.16  ! base temperature used for carbon decomposition
c
      bconst  = 10.0    ! maximum value of decomposition factor
c
c soil weighting factors
c
      rdepth  = 1. / (hsoi(1) + hsoi(2) + hsoi(3) + hsoi(4))
      rdepth2 = 1. / (hsoi(1) + hsoi(2))
c
c begin global grid
c
      do 100 i = 1, npoi
c
c ---------------------------------------------------------------------
c * * * daily water budget terms * * *
c ---------------------------------------------------------------------
c
        adrain(i)    = ((ndtimes-1) * adrain(i) + raina(i) * 86400.) * rwork
        adsnow(i)    = ((ndtimes-1) * adsnow(i) + snowa(i) * 86400.) * rwork
c
c Verification of the weather generator algorithm
c
        niter = int (86400./dtime)
c
        if (istep .eq. niter) then
c
          precipfac = precip(i) - adrain(i) - adsnow(i)
c
          if ((precipfac .lt. -0.1) .OR. (precipfac .gt. 0.1)) print *,
     >         'ERROR in sumday:', i, adrain(i) + adsnow(i), precip(i)
        endif
c
c End of verification
c
        adaet(i)     = ((ndtimes-1) * adaet(i)  - fvapa(i) * 86400.)
     >                 * rwork
c
        adtrunoff(i)  = ((ndtimes-1) * adtrunoff(i)  +
     >                    (grunof(i) + gdrain(i)) * 86400.) * rwork
        adsrunoff(i)  = ((ndtimes-1) * adsrunoff(i)  +
     >                     grunof(i)              * 86400.) * rwork
        addrainage(i) = ((ndtimes-1) * addrainage(i) +
     >                     gdrain(i)              * 86400.) * rwork
c
c ---------------------------------------------------------------------
c * * * daily atmospheric terms * * *
c ---------------------------------------------------------------------
c
        adrh(i) = ((ndtimes-1) * adrh(i) + rh(i)) * rwork
c
c ---------------------------------------------------------------------
c * * * daily snow parameters * * *
c ---------------------------------------------------------------------
c
        snodpth = hsno(i,1) + hsno(i,2) + hsno(i,3)
c
        adsnod(i) = ((ndtimes-1) * adsnod(i) + snodpth) * rwork
        adsnof(i) = ((ndtimes-1) * adsnof(i) + fi(i))   * rwork
c
c ---------------------------------------------------------------------
c * * * soil parameters * * *
c ---------------------------------------------------------------------
c
c initialize average soil parameters
c
        soiltemp = 0.0
        soilmois = 0.0
        soilice  = 0.0
c
	soitempc = 0.0
	soimoisc = 0.0
c
c averages for first 2 layers of soil
c
        do 110 k = 1, 2
          soiltemp =  soiltemp + tsoi(i,k)  * hsoi(k)
          soilmois =  soilmois + wsoi(i,k)  * hsoi(k)
          soilice  =  soilice  + wisoi(i,k) * hsoi(k)
 110    continue
c
c weighting on just thickness of each layer
c
        soilmois = soilmois * rdepth2
        soilice  = soilice  * rdepth2
        soiltemp = soiltemp * rdepth2
c
c calculate average root temperature, soil temperature and moisture and 
c ice content based on rooting profiles (weighted) from jackson et al
c 1996
c
c these soil moisture and temperatures are used in biogeochem.f 
c we assume that the rooting profiles approximate
c where carbon resides in the soil
c
        do 120 k = 1, nsoilay
          soitempc = soitempc + tsoi(i,k)  * 0.5 *
     >               (froot(k,1) + froot(k,2))
          soimoisc = soimoisc + wsoi(i,k)  * 0.5 *
     >               (froot(k,1) + froot(k,2))
 120    continue
c
c calculate daily average soil moisture and soil ice
c using thickness of each layer as weighting function
c
        adwsoi(i)  = ((ndtimes-1) * adwsoi(i)  + soilmois) * rwork
        adtsoi(i)  = ((ndtimes-1) * adtsoi(i)  + soiltemp) * rwork
        adwisoi(i) = ((ndtimes-1) * adwisoi(i) + soilice)  * rwork
c
c calculate daily average for soil temp/moisture of top layer
c
        adtlaysoi(i) = ((ndtimes-1) * adtlaysoi(i) + tsoi(i,1)) * rwork
        adwlaysoi(i) = ((ndtimes-1) * adwlaysoi(i) + wsoi(i,1)) * rwork
c
c calculate separate variables to keep track of weighting using 
c rooting profile information
c
c note that these variables are only used for diagnostic purposes
c and that they are not needed in the biogeochemistry code
c
        adwsoic(i)  = ((ndtimes-1) * adwsoic(i) + soimoisc) * rwork
        adtsoic(i)  = ((ndtimes-1) * adtsoic(i) + soitempc) * rwork
c
c ---------------------------------------------------------------------
c * * * calculate daily soil co2 fluxes * * *
c ---------------------------------------------------------------------
c
c increment daily total co2 respiration from microbes
c tco2mic is instantaneous value of co2 flux calculated in biogeochem.f
c
        adco2mic(i) = ((ndtimes-1) * adco2mic(i) +
     >                  tco2mic(i) * rwork3) * rwork
c
c increment daily total co2 respiration from fine roots
c tco2root is instantaneous value of co2 flux calculated in stats.f
c
        adco2root(i) = ((ndtimes-1) * adco2root(i) + 
     >                   tco2root(i) * rwork3) * rwork
c 
c calculate daily total co2 respiration from soil
c
        adco2soi(i)  = adco2root(i) + adco2mic(i)
c
c calculate daily ratio of total root to total co2 respiration
c
        if (adco2soi(i).gt.0.0) then
          adco2ratio(i) = adco2root(i) / adco2soi(i)
        else
          adco2ratio(i) = -999.99
        endif
c
c ---------------------------------------------------------------------
c * * * calculate daily litter decomposition parameters * * *
c ---------------------------------------------------------------------
c
c calculate litter carbon decomposition factors
c using soil temp, moisture and ice for top soil layer
c
c calculation of soil biogeochemistry decomposition factors 
c based on moisture and temperature affects on microbial
c biomass dynamics
c
c moisture function based on water-filled pore space (wfps)  
c williams et al., 1992 and friend et al., 1997 used in the
c hybrid 4.0 model; this is based on linn and doran, 1984
c
c temperature functions are derived from arrhenius function
c found in lloyd and taylor, 1994 with a 15 c base 
c
c calculate temperature decomposition factor
c CD impose lower limit to avoid division by zero at tsoi=227.13
c
        if (tsoi(i,1) .gt. 237.13) then
           factor = min (exp(tconst * ((1. / (btemp - 227.13)) - (1. /
     >          (tsoi(i,1)-227.13)))), bconst)
        else
           factor = exp(tconst * ((1. / (btemp - 227.13)) - (1. /
     >          (237.13-227.13))))
        end if
c
c calculate water-filled pore space (in percent)
c
c wsoi is relative to pore space not occupied by ice and water
c thus must include the ice fraction in the calculation
c	
        wfps = (1.0 - wisoi(i,1)) * wsoi(i,1) * 100.0	
c
c calculate moisture decomposition factor
c
        if (wfps .ge. 60.0) then
          moist = 0.000371 * (wfps**2) - (0.0748 * wfps) + 4.13
        else
          moist = exp((wfps - 60.0)**2 / (-800.0))	
        endif
c
c calculate combined temperature / moisture decomposition factor
c
        factor = max (0.001, min (bconst, factor * moist))
c
c calculate daily average litter decomposition factor
c
        decompl(i) = ((ndtimes-1) * decompl(i) + factor) * rwork
c
c ---------------------------------------------------------------------
c * * * calculate daily soil carbon decomposition parameters * * *
c ---------------------------------------------------------------------
c
c calculate soil carbon decomposition factors
c using soil temp, moisture and ice weighted by rooting profile scheme 
c
c calculation of soil biogeochemistry decomposition factors 
c based on moisture and temperature affects on microbial
c biomass dynamics
c
c moisture function based on water-filled pore space (wfps)  
c williams et al., 1992 and friend et al., 1997 used in the
c hybrid 4.0 model; this is based on linn and doran, 1984
c
c temperature functions are derived from arrhenius function
c found in lloyd and taylor, 1994 with a 15 c base 
c
c calculate temperature decomposition factor
c CD: impose lower limit to avoid division by zero at tsoi=227.13
c
        if (soiltemp .gt. 237.13) then
           factor = min (exp(tconst * ((1. / (btemp - 227.13)) - (1. /
     >          (soiltemp - 227.13)))), bconst)
        else
           factor = exp(tconst * ((1. / (btemp - 227.13)) - (1. /
     >          (237.13-227.13))))
        end if
c
c calculate water-filled pore space (in percent)
c
c wsoi is relative to pore space not occupied by ice and water
c thus must include the ice fraction in the calculation
c	
        wfps = (1. - soilice) * soilmois * 100.0	
c
c calculate moisture decomposition factor
c
        if (wfps .ge. 60.0) then
          moist = 0.000371 * (wfps**2) - (0.0748 * wfps) + 4.13
        else
          moist = exp((wfps - 60.0)**2 / (-800.0))
	endif
c
c calculate combined temperature / moisture decomposition factor
c
        factor = max (0.001, min (bconst, factor * moist))
c
c calculate daily average soil decomposition factor
c
        decomps(i) = ((ndtimes-1) * decomps(i) + factor) * rwork
c
c ---------------------------------------------------------------------
c * * * calculate other daily biogeochemical parameters * * *
c ---------------------------------------------------------------------
c
c increment daily total of net nitrogen mineralization
c value for tnmin is calculated in biogeochem.f
c
        adnmintot(i) = ((ndtimes-1) * adnmintot(i) +
     >                   tnmin(i) * rwork4) * rwork
c
 100  continue
c
c return to main program
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine summonth (istep, iday, imonth)
c ---------------------------------------------------------------------
c
c first convert to units that make sense for output
c
c   - convert all temperatures to deg c
c   - convert all liquid or vapor fluxes to mm/day
c   - redefine upwd directed heat fluxes as positive
c
c common blocks
c 
      include 'implicit.h'
c
      include 'compar.h'
      include 'comatm.h'
      include 'com1d.h'
      include 'comhyd.h'
      include 'comsno.h'
      include 'comsoi.h'
      include 'comveg.h'
      include 'comsum.h'
c
c Arguments (input)
c
      integer
     >  istep,     ! daily timestep number (passed in)
     >  iday,      ! day number  (passed in)
     >  imonth     ! month number (passed in)
c
c local variables
c
      integer
     >  i, k       ! loop indices
c
      real 
     >  rwork,     ! time work variable
     >  rwork2,    !
     >  rwork3,    !
     >  rwork4,    !
     >  rdepth,    ! 1/total soil depth over 4 1st layers
     >  solartot,  ! total incoming radiation (direct + diffuse, visible + nearIR)
     >  soiltemp,  ! average soil temp for 4 1st layers
     >  soilmois,  ! average soil moisture for 4 1st layers 
     >  soilice,   ! average soil ice for 4 1st layers 
     >  vwc,       ! total liquid + ice content of 4 1st layers
     >  awc,       ! total available water (+ ice) content of 4 1st layer
     >  snodpth    ! total snow depth
c
c
c ---------------------------------------------------------------------
c * * * update counters and working variables * * *
c ---------------------------------------------------------------------
c 
c if the first timestep of the month then reset averages
c
      if ((istep.eq.1).and.(iday.eq.1)) nmtimes = 0
c
c accumulate terms
c
      nmtimes = nmtimes + 1
c
c working variables
c
c rwork4 for conversion of nitrogen mineralization (moles)
c
      rwork  = 1. / float(nmtimes)
      rwork2 = float(ndaypm(imonth)) * 86400.
      rwork3 = float(ndaypm(imonth)) * 86400. * 12.e-3
      rwork4 = float(ndaypm(imonth)) * 86400. * 14.e-3
c
      rdepth = 1. / (hsoi(1) + hsoi(2) + hsoi(3) + hsoi(4))
c
c begin global grid
c
      do 100 i = 1, npoi
c
c ---------------------------------------------------------------------
c * * * monthly water budget terms * * *
c ---------------------------------------------------------------------
c
        amrain(i)    = ((nmtimes-1) * amrain(i) +
     >                   raina(i) * 86400.) * rwork
        amsnow(i)    = ((nmtimes-1) * amsnow(i) +
     >                   snowa(i) * 86400.) * rwork
        amaet(i)     = ((nmtimes-1) * amaet(i)  -
     >                   fvapa(i) * 86400.) * rwork
c
        amtrunoff(i)  = ((nmtimes-1) * amtrunoff(i)  +
     >                   (grunof(i) + gdrain(i)) * 86400.) * rwork
        amsrunoff(i)  = ((nmtimes-1) * amsrunoff(i)  +
     >                    grunof(i)              * 86400.) * rwork
        amdrainage(i) = ((nmtimes-1) * amdrainage(i) +
     .                    gdrain(i)              * 86400.) * rwork
c
c ---------------------------------------------------------------------
c * * * monthly atmospheric terms * * *
c ---------------------------------------------------------------------
c
        amtemp(i)  = ((nmtimes-1) * amtemp(i)  + ta(i) - 273.16)  * rwork
        amcloud(i) = ((nmtimes-1) * amcloud(i) + cloud(i) * 100.) * rwork
        amqa(i)    = ((nmtimes-1) * amqa(i)    + qa(i))           * rwork
        amrh(i)    = ((nmtimes-1) * amrh(i)    + rh(i))           * rwork
c
c ---------------------------------------------------------------------
c * * * energy budget terms * * *
c ---------------------------------------------------------------------
c
        solartot = solad(i,1) + solad(i,2) + solai(i,1) + solai(i,2)
c
        amsolar(i)  = ((nmtimes-1) * amsolar(i)  +
     >                  solartot)         * rwork
        amirup(i)   = ((nmtimes-1) * amirup(i)   +
     >                  firb(i))         * rwork
        amirdown(i) = ((nmtimes-1) * amirdown(i) +
     >                  fira(i))         * rwork
        amsens(i)   = ((nmtimes-1) * amsens(i)   -
     >                  fsena(i))        * rwork
        amlatent(i) = ((nmtimes-1) * amlatent(i) -
     >                  fvapa(i) * hvap) * rwork
c
c ---------------------------------------------------------------------
c * * * monthly vegetation parameters * * *
c ---------------------------------------------------------------------
c
        amlaiu(i) = ((nmtimes-1) * amlaiu(i) + fu(i) * lai(i,2)) * rwork
        amlail(i) = ((nmtimes-1) * amlail(i) + fl(i) * lai(i,1)) * rwork
c
c ---------------------------------------------------------------------
c * * * monthly soil parameters * * *
c ---------------------------------------------------------------------
c
        soiltemp = 0.0
        soilmois = 0.0
        soilice  = 0.0
c
        vwc = 0.0
        awc = 0.0
c
c averages for first 4 layers of soil (assumed to add to 1 meter depth)
c
        do 110 k = 1, 4
c
          soiltemp =  soiltemp + tsoi(i,k)  * hsoi(k)
          soilmois =  soilmois + wsoi(i,k)  * hsoi(k)
          soilice  =  soilice  + wisoi(i,k) * hsoi(k)
c
          vwc = vwc + (wisoi(i,k) + (1. - wisoi(i,k)) * wsoi(i,k)) *
     >                hsoi(k) * poros(i,k)
c
          awc = awc + max (0.0, (wisoi(i,k) +
     >                (1. - wisoi(i,k)) * wsoi(i,k)) - swilt(i,k)) *
     >                hsoi(k) * poros(i,k) * 100.0
c
 110    continue
c
        soiltemp = soiltemp * rdepth - 273.16
        soilmois = soilmois * rdepth
        soilice  = soilice  * rdepth
c
        vwc = vwc * rdepth
        awc = awc * rdepth
c
        amtsoi(i)  = ((nmtimes-1) * amtsoi(i)  + soiltemp) * rwork
        amwsoi(i)  = ((nmtimes-1) * amwsoi(i)  + soilmois) * rwork
        amwisoi(i) = ((nmtimes-1) * amwisoi(i) + soilice)  * rwork
        amvwc(i)   = ((nmtimes-1) * amvwc(i)   + vwc)      * rwork
        amawc(i)   = ((nmtimes-1) * amawc(i)   + awc)      * rwork
c
c ---------------------------------------------------------------------
c * * * snow parameters * * *
c ---------------------------------------------------------------------
c
        snodpth = hsno(i,1) + hsno(i,2) + hsno(i,3)
c
        amsnod(i) = ((nmtimes-1) * amsnod(i) + snodpth) * rwork
        amsnof(i) = ((nmtimes-1) * amsnof(i) + fi(i))   * rwork
c
c ---------------------------------------------------------------------
c * * * determine monthly npp * * *
c ---------------------------------------------------------------------
c
        amnpp(i,1)  = ((nmtimes-1) * amnpp(i,1)  +
     >                  tnpp(i,1)  * rwork3) * rwork
        amnpp(i,2)  = ((nmtimes-1) * amnpp(i,2)  +
     >                  tnpp(i,2)  * rwork3) * rwork
        amnpp(i,3)  = ((nmtimes-1) * amnpp(i,3)  +
     >                  tnpp(i,3)  * rwork3) * rwork
        amnpp(i,4)  = ((nmtimes-1) * amnpp(i,4)  +
     >                  tnpp(i,4)  * rwork3) * rwork
        amnpp(i,5)  = ((nmtimes-1) * amnpp(i,5)  +
     >                  tnpp(i,5)  * rwork3) * rwork
        amnpp(i,6)  = ((nmtimes-1) * amnpp(i,6)  +
     >                  tnpp(i,6)  * rwork3) * rwork
        amnpp(i,7)  = ((nmtimes-1) * amnpp(i,7)  +
     >                  tnpp(i,7)  * rwork3) * rwork
        amnpp(i,8)  = ((nmtimes-1) * amnpp(i,8)  +
     >                  tnpp(i,8)  * rwork3) * rwork
        amnpp(i,9)  = ((nmtimes-1) * amnpp(i,9)  +
     >                  tnpp(i,9)  * rwork3) * rwork
        amnpp(i,10) = ((nmtimes-1) * amnpp(i,10) +
     >                  tnpp(i,10) * rwork3) * rwork
        amnpp(i,11) = ((nmtimes-1) * amnpp(i,11) +
     >                  tnpp(i,11) * rwork3) * rwork
        amnpp(i,12) = ((nmtimes-1) * amnpp(i,12) +
     >                  tnpp(i,12) * rwork3) * rwork
c
        amnpptot(i) = amnpp(i,1)  + amnpp(i,2)  + amnpp(i,3)  + 
     >                amnpp(i,4)  + amnpp(i,5)  + amnpp(i,6)  + 
     >                amnpp(i,7)  + amnpp(i,8)  + amnpp(i,9)  +
     >                amnpp(i,10) + amnpp(i,11) + amnpp(i,12)
c
c ---------------------------------------------------------------------
c * * * monthly biogeochemistry parameters * * *
c ---------------------------------------------------------------------
c
c increment monthly total co2 respiration from microbes
c tco2mic is instantaneous value of co2 flux calculated in biogeochem.f
c
        amco2mic(i) = ((nmtimes-1) * amco2mic(i) +
     >                  tco2mic(i) * rwork3) * rwork
c
c increment monthly total co2 respiration from roots
c tco2root is instantaneous value of co2 flux calculated in stats.f
c
        amco2root(i) = ((nmtimes-1) * amco2root(i) +
     >                   tco2root(i) * rwork3) * rwork
c
c calculate average total co2 respiration from soil
c
        amco2soi(i)  = amco2root(i) + amco2mic(i)
c  
c  calculate ratio of root to total co2 respiration
c
        if (amco2soi(i).gt.0.0) then
          amco2ratio(i) = amco2root(i) / amco2soi(i)
        else
          amco2ratio(i) = -999.99
        endif
c 
c  monthly net ecosystem co2 flux -- npp total minus microbial respiration 
c  the npp total includes losses from root respiration
c
        amneetot(i)  = amnpptot(i) - amco2mic(i) 
c
c increment monthly total of net nitrogen mineralization
c value for tnmin is calculated in biogeochem.f
c
        amnmintot(i) = ((nmtimes-1) * amnmintot(i) + tnmin(i) *
     >                   rwork4) * rwork
c
 100  continue
c
c return to main program
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine sumyear (istep, iday, imonth)
c ---------------------------------------------------------------------
c
c common blocks
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'comatm.h'
      include 'com1d.h'
      include 'comhyd.h'
      include 'comsno.h'
      include 'comsoi.h'
      include 'comsum.h'
      include 'comveg.h'
c
c Arguments (input)
c
      integer
     >  istep,     ! daily timestep number  (passed in)
     >  iday,      ! day number  (passed in)
     >  imonth     ! month number (passed in)
c
c local variables
c
      integer
     >  i,k        ! loop indices
c
      real 
     >  rwork,     !
     >  rwork2,    !
     >  rwork3,    !
     >  rwork4,    !
     >  rdepth,    ! 1/total soil depth over 4 1st layers
     >  solartot,  ! total incoming radiation (direct + diffuse, visible + nearIR)
     >  soiltemp,  ! average soil temp for 4 1st layers
     >  soilmois,  ! average soil moisture for 4 1st layers 
     >  soilice,   ! average soil ice for 4 1st layers 
     >  vwc,       ! total liquid + ice content of 4 1st layers
     >  awc,       ! total available water (+ ice) content of 4 1st layer
     >  water,     ! fire factor: total water content of 1st layer (liquid+ice)
     >  waterfrac, ! fire factor: available water content of 1st layer
     >  fueldry,   ! fire factor
     >  allroots,  ! annual average root biomass
     >  wtotp      ! total water stored in soil+vegetation+snow at previous timestep
c
c ---------------------------------------------------------------------
c * * * update counters and working variables * * *
c ---------------------------------------------------------------------
c
c reset sumyear if the first timestep of the year
c
      if ((istep.eq.1).and.(iday.eq.1).and.(imonth.eq.1)) nytimes = 0
c
c accumulate yearly output
c
      nytimes = nytimes + 1
c
c working variables
c
c rwork4 is for nitrogen mineralization conversion
c
      rwork  = 1. / float(nytimes)
      rwork2 = float(ndaypy) * 86400.
      rwork3 = float(ndaypy) * 86400. * 12.e-3
      rwork4 = float(ndaypy) * 86400. * 14.e-3
c
      rdepth = 1. / (hsoi(1) + hsoi(2) + hsoi(3) + hsoi(4))
c
c begin global grid
c
      do 100 i = 1, npoi
c
c ---------------------------------------------------------------------
c * * * annual energy budget terms * * *
c ---------------------------------------------------------------------
c
        solartot = solad(i,1) + solad(i,2) + solai(i,1) + solai(i,2)
c
        aysolar(i)  = ((nytimes-1) * aysolar(i)  + solartot) *     
     >                  rwork
        ayirup(i)   = ((nytimes-1) * ayirup(i)   + firb(i))  *
     >                  rwork
        ayirdown(i) = ((nytimes-1) * ayirdown(i) + fira(i))  *
     >                  rwork
        aysens(i)   = ((nytimes-1) * aysens(i)   - fsena(i)) *
     >                  rwork
        aylatent(i) = ((nytimes-1) * aylatent(i) - fvapa(i)  * hvap) *
     >                  rwork
c
c ---------------------------------------------------------------------
c * * * annual water budget terms * * *
c ---------------------------------------------------------------------
c
        ayprcp(i)     = ((nytimes-1) * ayprcp(i)  +
     >                   (raina(i) + snowa(i)) * rwork2) * rwork
        ayaet(i)      = ((nytimes-1) * ayaet(i)   -
     >                    fvapa(i)             * rwork2) * rwork
        aytrans(i)    = ((nytimes-1) * aytrans(i) +
     >                    gtrans(i)            * rwork2) * rwork
c
        aytrunoff(i)  = ((nytimes-1) * aytrunoff(i)  +
     >                   (grunof(i) + gdrain(i)) * rwork2) * rwork
        aysrunoff(i)  = ((nytimes-1) * aysrunoff(i)  +
     >                    grunof(i)              * rwork2) * rwork
        aydrainage(i) = ((nytimes-1) * aydrainage(i) +
     >                    gdrain(i)  * rwork2)   * rwork
c
c---------------------------------------------------------------------
c CD
c estimate the change in soil-vegetation water content. Used to check 
c mass conservation
c---------------------------------------------------------------------
c
        wtotp = wtot(i)
c
        wtot(i) = (wliqu(i)+wsnou(i)) * fu(i) * 2.0 * lai(i,2) +
     >            (wliqs(i)+wsnos(i)) * fu(i) * 2.0 * sai(i,2) +
     >            (wliql(i)+wsnol(i)) * fl(i) * 2.0 *
     >            (lai(i,1) + sai(i,1)) * (1. - fi(i))
c
        wtot(i) = wtot(i) + wpud(i) + wipud(i)
c
        do 10 k = 1, nsoilay
          wtot(i) = wtot(i) +
     >              poros(i,k)*wsoi(i,k)*(1.-wisoi(i,k))*hsoi(k)*rhow +
     >              poros(i,k)*wisoi(i,k)*hsoi(k)*rhow
 10     continue
c
        do 20 k = 1, nsnolay
          wtot(i) = wtot(i) + fi(i)*rhos*hsno(i,k)
 20     continue
c
        aydwtot(i) = ((nytimes-1) * aydwtot(i) +
     >                  wtot(i) - wtotp) * rwork
c
c ---------------------------------------------------------------------
c * * * annual soil parameters * * *
c ---------------------------------------------------------------------
c
        soiltemp = 0.0
        soilmois = 0.0
        soilice  = 0.0
c
        vwc = 0.0
        awc = 0.0
c
c averages for first 4 layers of soil
c
        do 110 k = 1, 4
c
          soiltemp =  soiltemp + tsoi(i,k)  * hsoi(k)
          soilmois =  soilmois + wsoi(i,k)  * hsoi(k)
          soilice  =  soilice  + wisoi(i,k) * hsoi(k)
c
          vwc = vwc + (wisoi(i,k) + (1. - wisoi(i,k)) * wsoi(i,k)) *
     >                hsoi(k) * poros(i,k)
c
          awc = awc + max (0.0, (wisoi(i,k) +
     >                (1. - wisoi(i,k)) * wsoi(i,k)) - swilt(i,k)) *
     >                hsoi(k) * poros(i,k) * 100.0
c
 110    continue
c
c average soil and air temperatures
c
        soiltemp = soiltemp * rdepth - 273.16
        soilmois = soilmois * rdepth
        soilice  = soilice  * rdepth
c
        vwc = vwc * rdepth
        awc = awc * rdepth
c
c annual average soil moisture and soil ice
c
        aywsoi(i)  = ((nytimes-1) * aywsoi(i)  + soilmois) * rwork
        aywisoi(i) = ((nytimes-1) * aywisoi(i) + soilice)  * rwork
        aytsoi(i)  = ((nytimes-1) * aytsoi(i)  + soiltemp) * rwork
        ayvwc(i)   = ((nytimes-1) * ayvwc(i)   + vwc)      * rwork
        ayawc(i)   = ((nytimes-1) * ayawc(i)   + awc)      * rwork
c
c soil moisture stress
c
        aystresstu(i) = rwork * ((nytimes-1) * aystresstu(i) +
     >                  stresstu(i))
c
        aystresstl(i) = rwork * ((nytimes-1) * aystresstl(i) +
     >                  stresstl(i))
c
c ---------------------------------------------------------------------
c * * * determine annual gpp * * *
c ---------------------------------------------------------------------
c
c gross primary production of each plant type
c
        aygpp(i,1)  = ((nytimes-1) * aygpp(i,1)  +
     >                  tgpp(i,1)  * rwork3) * rwork
        aygpp(i,2)  = ((nytimes-1) * aygpp(i,2)  +
     >                  tgpp(i,2)  * rwork3) * rwork
        aygpp(i,3)  = ((nytimes-1) * aygpp(i,3)  +
     >                  tgpp(i,3)  * rwork3) * rwork
        aygpp(i,4)  = ((nytimes-1) * aygpp(i,4)  +
     >                  tgpp(i,4)  * rwork3) * rwork
        aygpp(i,5)  = ((nytimes-1) * aygpp(i,5)  +
     >                  tgpp(i,5)  * rwork3) * rwork
        aygpp(i,6)  = ((nytimes-1) * aygpp(i,6)  +
     >                  tgpp(i,6)  * rwork3) * rwork
        aygpp(i,7)  = ((nytimes-1) * aygpp(i,7)  +
     >                  tgpp(i,7)  * rwork3) * rwork
        aygpp(i,8)  = ((nytimes-1) * aygpp(i,8)  +
     >                  tgpp(i,8)  * rwork3) * rwork
        aygpp(i,9)  = ((nytimes-1) * aygpp(i,9)  +
     >                  tgpp(i,9)  * rwork3) * rwork
        aygpp(i,10) = ((nytimes-1) * aygpp(i,10) +
     >                  tgpp(i,10) * rwork3) * rwork
        aygpp(i,11) = ((nytimes-1) * aygpp(i,11) +
     >                  tgpp(i,11) * rwork3) * rwork
        aygpp(i,12) = ((nytimes-1) * aygpp(i,12) +
     >                  tgpp(i,12) * rwork3) * rwork
c
c gross primary production of the entire gridcell
c
        aygpptot(i) = aygpp(i,1)  + aygpp(i,2)  + aygpp(i,3)  + 
     >                aygpp(i,4)  + aygpp(i,5)  + aygpp(i,6)  + 
     >                aygpp(i,7)  + aygpp(i,8)  + aygpp(i,9)  +
     >                aygpp(i,10) + aygpp(i,11) + aygpp(i,12)
c
c ---------------------------------------------------------------------
c * * * determine annual npp * * *
c ---------------------------------------------------------------------
c
c net primary production of each plant type
c
        aynpp(i,1)  = ((nytimes-1) * aynpp(i,1)  +
     >                  tnpp(i,1)  * rwork3) * rwork
        aynpp(i,2)  = ((nytimes-1) * aynpp(i,2)  +
     >                  tnpp(i,2)  * rwork3) * rwork
        aynpp(i,3)  = ((nytimes-1) * aynpp(i,3)  +
     >                  tnpp(i,3)  * rwork3) * rwork
        aynpp(i,4)  = ((nytimes-1) * aynpp(i,4)  +
     >                  tnpp(i,4)  * rwork3) * rwork
        aynpp(i,5)  = ((nytimes-1) * aynpp(i,5)  +
     >                  tnpp(i,5)  * rwork3) * rwork
        aynpp(i,6)  = ((nytimes-1) * aynpp(i,6)  +
     >                  tnpp(i,6)  * rwork3) * rwork
        aynpp(i,7)  = ((nytimes-1) * aynpp(i,7)  +
     >                  tnpp(i,7)  * rwork3) * rwork
        aynpp(i,8)  = ((nytimes-1) * aynpp(i,8)  +
     >                  tnpp(i,8)  * rwork3) * rwork
        aynpp(i,9)  = ((nytimes-1) * aynpp(i,9)  +
     >                  tnpp(i,9)  * rwork3) * rwork
        aynpp(i,10) = ((nytimes-1) * aynpp(i,10) +
     >                  tnpp(i,10) * rwork3) * rwork
        aynpp(i,11) = ((nytimes-1) * aynpp(i,11) +
     >                  tnpp(i,11) * rwork3) * rwork
        aynpp(i,12) = ((nytimes-1) * aynpp(i,12) +
     >                  tnpp(i,12) * rwork3) * rwork
c
c net primary production of the entire gridcell
c
        aynpptot(i) = aynpp(i,1)  + aynpp(i,2)  + aynpp(i,3)  + 
     >                aynpp(i,4)  + aynpp(i,5)  + aynpp(i,6)  + 
     >                aynpp(i,7)  + aynpp(i,8)  + aynpp(i,9)  +
     >                aynpp(i,10) + aynpp(i,11) + aynpp(i,12)
c
c ---------------------------------------------------------------------
c * * * annual carbon budget terms * * *
c ---------------------------------------------------------------------
c
c fire factor used in vegetation dynamics calculations
c
        water     = wisoi(i,1) + (1. - wisoi(i,1)) * wsoi(i,1)
        waterfrac = (water - swilt(i,1)) / (1. - swilt(i,1))
c
        fueldry = max (0.0, min (1.0, -2.0 * (waterfrac - 0.5)))
c
        firefac(i) = ((nytimes-1) * firefac(i) + fueldry) * rwork
c
c increment annual total co2 respiration from microbes
c tco2mic is instantaneous value of co2 flux calculated in biogeochem.f
c
        ayco2mic(i) = ((nytimes-1) * ayco2mic(i) +
     >                  tco2mic(i) * rwork3) * rwork
c
c increment annual total co2 respiration from roots
c
        ayco2root(i) = ((nytimes-1) * ayco2root(i) +
     >                   tco2root(i) * rwork3) * rwork
c
c calculate annual total co2 respiration from soil
c
        ayco2soi(i)  = ayco2root(i) + ayco2mic(i)
c  
c annual net ecosystem co2 flux -- npp total minus microbial respiration 
c the npp total includes losses from root respiration
c
        ayneetot(i)  = aynpptot(i) - ayco2mic(i) 
c
c annual average root biomass
c
        allroots = cbior(i,1)  + cbior(i,2)  + cbior(i,3)  +
     >             cbior(i,4)  + cbior(i,5)  + cbior(i,6)  +
     >             cbior(i,7)  + cbior(i,8)  + cbior(i,9)  +
     >             cbior(i,10) + cbior(i,11) + cbior(i,12)
c
        ayrootbio(i) = ((nytimes-1) * ayrootbio(i) + allroots) * rwork
c
c ---------------------------------------------------------------------
c * * * annual biogeochemistry terms * * *
c ---------------------------------------------------------------------
c
c increment annual total of net nitrogen mineralization
c value for tnmin is calculated in biogeochem.f
c
        aynmintot(i) = ((nytimes-1) * aynmintot(i) +
     >                   tnmin(i) * rwork4) * rwork
c
c other biogeochemistry variables
c
        ayalit(i)  = ((nytimes-1) * ayalit(i)  + totalit(i))  * rwork
        ayblit(i)  = ((nytimes-1) * ayblit(i)  + totrlit(i))  * rwork
        aycsoi(i)  = ((nytimes-1) * aycsoi(i)  + totcsoi(i))  * rwork
        aycmic(i)  = ((nytimes-1) * aycmic(i)  + totcmic(i))  * rwork
        ayanlit(i) = ((nytimes-1) * ayanlit(i) + totanlit(i)) * rwork
        aybnlit(i) = ((nytimes-1) * aybnlit(i) + totrnlit(i)) * rwork
        aynsoi(i)  = ((nytimes-1) * aynsoi(i)  + totnsoi(i))  * rwork
c
 100  continue
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine gdiag (iyear, iyear0)
c ---------------------------------------------------------------------
c
c common blocks
c 
      include 'implicit.h'
c
      include 'compar.h'
      include 'comsum.h'
      include 'comveg.h'
c
c Arguments (input)
c
      integer iyear,   ! year counter
     >        iyear0   ! first year of simulation
c
c local variables
c
      integer i        ! loop indice
c
      real gnee,       ! domain total nee (gt-c/yr)
     >     gnpp,       ! domain total npp (gt-c/yr)
     >     ggpp,       ! domain total gpp (gt-c/yr)
     >     gbiomass,   ! domain total biomass (gt-c)
     >     galitc,     ! domain total aboveground litter carbon (gt-c)
     >     gblitc,     ! domain total belowground litter carbon (gt-c)
     >     gsoic,      ! domain total soil carbon (gt-c)
     >     gco2soi,    ! domain total soil surface co2 flux (gt-c)
     >     galitn,     ! domain total aboveground litter nitrogen (gt-c)
     >     gblitn,     ! domain total belowground litter nitrogen (gt-c)
     >     gsoin,      ! domain total soil nitrogen (gt-c)
     >     gprcp,      ! domain average annual precipitation (mm/yr)
     >     gaet,       ! domain average annual evapotranspiration (mm/yr)
     >     gt,         ! domain average annual transpiration (mm/yr)
     >     gtrunoff,   ! domain average total runoff (mm/yr)
     >     gsrunoff,   ! domain average surface runoff (mm/yr)
     >     gdrainage,  ! domain average drainage (mm/yr)
     >     gdwtot,     !   "      "     water recharge (mm/yr)
     >     gtarea,     ! total land area of the domain (m**2)
     >     aratio,     ! aet / prcp ratio
     >     rratio,     ! runoff / prcp ratio
     >     sratio,     ! surface runoff / drainage ratio
     >     tratio      ! transpiration / aet ratio
c
c initialize variables
c
      gtarea    = 0.0
      gnee      = 0.0
      gnpp      = 0.0
      ggpp      = 0.0
      gbiomass  = 0.0
      galitc    = 0.0
      gblitc    = 0.0
      gsoic     = 0.0
      gco2soi   = 0.0
      galitn    = 0.0
      gblitn    = 0.0
      gsoin     = 0.0
      gprcp     = 0.0
      gaet      = 0.0
      gt        = 0.0
      gtrunoff  = 0.0
      gsrunoff  = 0.0
      gdrainage = 0.0
      gdwtot    = 0.0
c
      do 100 i = 1, npoi
c
        ayrratio(i) = min (1.0, max (0.0, aytrunoff(i)) /
     >                max (0.1, ayprcp(i)))
c
        aytratio(i) = min (1.0, max (0.0, aytrans(i))   /
     >                max (0.1, ayaet(i)))
c
        gtarea    = gtarea    + garea(i)
c
        gnee      = gnee      + garea(i) * ayneetot(i) * 1.e-12
        gnpp      = gnpp      + garea(i) * aynpptot(i) * 1.e-12
        ggpp      = ggpp      + garea(i) * aygpptot(i) * 1.e-12
        gbiomass  = gbiomass  + garea(i) * totbiou(i)  * 1.e-12
     >                        + garea(i) * totbiol(i)  * 1.e-12
        galitc    = galitc    + garea(i) * ayalit(i)   * 1.e-12
        gblitc    = gblitc    + garea(i) * ayblit(i)   * 1.e-12
        gsoic     = gsoic     + garea(i) * aycsoi(i)   * 1.e-12
        gco2soi   = gco2soi   + garea(i) * ayco2soi(i) * 1.e-12
        galitn    = galitn    + garea(i) * ayanlit(i)  * 1.e-12
        gblitn    = gblitn    + garea(i) * aybnlit(i)  * 1.e-12
        gsoin     = gsoin     + garea(i) * aynsoi(i)   * 1.e-12
c
        gprcp     = gprcp     + garea(i) * ayprcp(i)
        gaet      = gaet      + garea(i) * ayaet(i)
        gt        = gt        + garea(i) * aytrans(i)
        gtrunoff  = gtrunoff  + garea(i) * aytrunoff(i)
        gsrunoff  = gsrunoff  + garea(i) * aysrunoff(i)
        gdrainage = gdrainage + garea(i) * aydrainage(i)
        gdwtot    = gdwtot    + garea(i) * aydwtot(i)*nytimes
c
 100  continue
c
      gprcp     = gprcp     / gtarea
      gaet      = gaet      / gtarea
      gt        = gt        / gtarea
      gtrunoff  = gtrunoff  / gtarea
      gsrunoff  = gsrunoff  / gtarea
      gdrainage = gdrainage / gtarea
      gdwtot    = gdwtot    / gtarea
c
      aratio   = gaet     / gprcp
      rratio   = gtrunoff / gprcp
      sratio   = gsrunoff / gtrunoff
      tratio   = gt       / gaet
c
      write (*,*) ' '
      write (*,*) '* * * annual diagnostic fields * * *'
      write (*,*) ' '
      write (*,9001) gnee
      write (*,9000) gnpp
      write (*,9002) ggpp
      write (*,9010) gbiomass
      write (*,9020) galitc
      write (*,9021) gblitc
      write (*,9030) gsoic
      write (*,9032) gco2soi
      write (*,9034) galitn
      write (*,9036) gblitn
      write (*,9038) gsoin
      write (*,*) ' '
      write (*,9040) gprcp
      write (*,9050) gaet
      write (*,9060) gt
      write (*,9070) gtrunoff
      write (*,9080) gsrunoff
      write (*,9090) gdrainage
      write (*,9095) gdwtot
      write (*,*) ' '
      write (*,9100) aratio
      write (*,9110) rratio
      write (*,*) ' '
      write (*,9120) tratio
      write (*,9130) sratio
      write (*,*) ' '
c
c write some diagnostic output to history file
c
      if (iyear.eq.iyear0) then
c
        open (20,file='ibis.out.global',status='unknown')
c
        write (20,*) ' '
        write (20,*) '* * * annual diagnostic fields * * *'
        write (20,*) ' '
        write (20,*) 
     >    'year       nee       npp       gpp   biomass   scarbon '// 
     >               'snitrogen   alitter   blitter    co2soi    '  //
     >               'aratio    rratio    tratio'
c
      endif
c
      write (20,9500) iyear, gnee, gnpp, ggpp, gbiomass, gsoic, 
     >                gsoin, galitc,
     >                gblitc, gco2soi, 100.0 * aratio, 100.0 * rratio,
     >                100.0 * tratio
c
      call flush (20)
c
c     close (20)
c
 9000 format (1x,'total npp             of the domain (gt-c/yr) : ',
     >        f12.3)
 9001 format (1x,'total nee             of the domain (gt-c/yr) : ',
     >        f12.5)
 9002 format (1x,'total gpp             of the domain (gt-c/yr) : ',
     >        f12.3)
 9010 format (1x,'total biomass         of the domain (gt-c)    : ',
     >        f12.3)
 9020 format (1x,'aboveground litter    of the domain (gt-c)    : ',
     >        f12.3)
 9021 format (1x,'belowground litter    of the domain (gt-c)    : ',
     >        f12.3)
 9030 format (1x,'total soil carbon     of the domain (gt-c)    : ',
     >        f12.3)
 9032 format (1x,'total soil co2 flux   of the domain (gt-c)    : ',
     >        f12.3)
 9034 format (1x,'aboveground litter n  of the domain (gt-c)    : ',
     >        f12.3)
 9036 format (1x,'belowground litter n  of the domain (gt-c)    : ',
     >        f12.3)
 9038 format (1x,'total soil nitrogen   of the domain (gt-c)    : ',
     >        f12.3)
 9040 format (1x,'average precipitation of the domain (mm/yr)   : ',
     >        f12.3)
 9050 format (1x,'average aet           of the domain (mm/yr)   : ',
     >        f12.3)
 9060 format (1x,'average transpiration of the domain (mm/yr)   : ',
     >        f12.3)
 9070 format (1x,'average runoff        of the domain (mm/yr)   : ',
     >        f12.3)
 9080 format (1x,'average surf runoff   of the domain (mm/yr)   : ',
     >        f12.3)
 9090 format (1x,'average drainage      of the domain (mm/yr)   : ',
     >        f12.3)
 9095 format (1x,'average moisture recharge of the domain (mm/yr) : ',
     >        f12.3)
 9100 format (1x,'total aet      / precipitation                : ',
     >        f12.3)
 9110 format (1x,'total runoff   / precipitation                : ',
     >        f12.3)
 9120 format (1x,'transpiration  / total aet                    : ',
     >        f12.3)
 9130 format (1x,'surface runoff / total runoff                 : ',
     >        f12.3)
 9500 format (1x,i4,12f10.2)
c
c return to main program
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine vdiag (iyear, iyear0)
c ---------------------------------------------------------------------
c
c common blocks
c 
      include 'implicit.h'
c
      include 'compar.h'
      include 'comsum.h'
      include 'comveg.h'
c
c Arguments (input)
c
      integer iyear,      ! year counter
     >        iyear0      ! first year of simulation
c
c local variables
c
      integer i, k        ! loop indices
c
      real vtarea(15),    ! total area of the vegetation type (m**2)
     >     vnee(15),      ! vegetation type average nee (kg-c/m**2/yr)
     >     vnpp(15),      ! vegetation type average npp (kg-c/m**2/yr)
     >     vgpp(15),      ! vegetation type average gpp (kg-c/m**2/yr)
     >     vbiomass(15),  ! vegetation type average biomass (kg-c/m**2)
     >     vlai(15),      ! vegetation type average lai (m**2/m**2)
     >     vsoic(15),     ! vegetation type average soil carbon (kg-c/m**2)
     >     vrunoff(15)    ! vegetation type average runoff (mm/yr)
c
c
c initialize variables
c
      do 100 k = 1, 15
c
        vtarea(k)   = 0.0
        vnee(k)     = 0.0
        vnpp(k)     = 0.0
        vgpp(k)     = 0.0
        vbiomass(k) = 0.0
        vlai(k)     = 0.0
        vsoic(k)    = 0.0
        vrunoff(k)  = 0.0
c
 100  continue
c
c sum ecosystem properties over each vegetation type
c
      do 200 i = 1, npoi
c
        k = int (max (1.0, min (15.0, vegtype0(i))))
c
        vtarea(k)   = vtarea(k)   + garea(i)
c
        vnee(k)     = vnee(k)     + garea(i) * ayneetot(i)
        vnpp(k)     = vnpp(k)     + garea(i) * aynpptot(i)
        vgpp(k)     = vgpp(k)     + garea(i) * aygpptot(i)
c
        vbiomass(k) = vbiomass(k) + garea(i) * totbiou(i)
     >                            + garea(i) * totbiol(i)
c
        vlai(k)     = vlai(k)     + garea(i) * totlaiu(i) 
     >                            + garea(i) * totlail(i)
c
        vsoic(k)    = vsoic(k)    + garea(i) * aycsoi(i) 
c
        vrunoff(k)  = vrunoff(k)  + garea(i) * aytrunoff(i)
c
 200  continue
c
c calculate area averages
c
      do 300 k = 1, 15
c
        vnee(k)     = vnee(k)     / max (1.0, vtarea(k))
        vnpp(k)     = vnpp(k)     / max (1.0, vtarea(k))
        vgpp(k)     = vgpp(k)     / max (1.0, vtarea(k))
        vbiomass(k) = vbiomass(k) / max (1.0, vtarea(k))
        vlai(k)     = vlai(k)     / max (1.0, vtarea(k))
        vsoic(k)    = vsoic(k)    / max (1.0, vtarea(k))
        vrunoff(k)  = vrunoff(k)  / max (1.0, vtarea(k))
c
 300  continue
c
c write some diagnostic output to history file
c
      if (iyear.eq.iyear0) then
        open (30,file='ibis.out.vegtype',status='unknown')
      endif
c
      write (30,*) ' '
      write (30,*) '* * annual diagnostic fields by vegetation type * *'
      write (30,*) ' '
      write (30,*) 
     > 'year    veg           area       nee       npp       gpp   '// 
     >             'biomass       lai   scarbon    runoff '
c
      do 400 k = 1, 15
c
        write (30,9000) 
     >     iyear, k, vtarea(k) / 1.0e+06, vnee(k), vnpp(k), vgpp(k),
     >                  vbiomass(k), vlai(k), vsoic(k), vrunoff(k)
c
 400  continue
c
      call flush (30)
c
c format statements
c
 9000 format (1x,i4,5x,i2,5x,1e10.3,7f10.3)
c
c return to main program
c
      return
      end
