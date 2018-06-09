c
c #####   #    #   #   #   ####      #     ####   #        ####    ####    #   #
c #    #  #    #    # #   #          #    #    #  #       #    #  #    #    # #
c #    #  ######     #     ####      #    #    #  #       #    #  #          #
c #####   #    #     #         #     #    #    #  #       #    #  #  ###     #
c #       #    #     #    #    #     #    #    #  #       #    #  #    #     #
c #       #    #     #     ####      #     ####   ######   ####    ####      #
c
c ---------------------------------------------------------------------
      subroutine stomata
c ---------------------------------------------------------------------
c
c common blocks
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'comatm.h'
      include 'comsno.h'
      include 'comsoi.h'
      include 'comveg.h'
      include 'com1d.h'
      include 'comsum.h'
      include 'compft.h'
      
c
c local variables
c
      integer i
c
      real rwork,  ! 3.47e-03 - 1. / tu(i)
     >     tau,    ! 
     >     tleaf,  ! leaf temp in celcius
     >     tempvm, !
     >     zweight !

      real esat12, ! vapor pressure in upper canopy air 
     >     qsat12, ! specific humidity in upper canopy air
     >     rh12,   ! relative humidity in upper canopy air 
     >     esat34, ! vapor pressure in lower canopy air
     >     qsat34, ! specific humidity in lower canopy air 
     >     rh34,   ! relative humidity in lower canopy air 
     >     gbco2u, ! bound. lay. conductance for CO2 in upper canopy
     >     gbco2l, ! bound. lay. conductance for CO2 in lower canopy
     >     gscub,  ! 
     >     gscuc,  !
     >     gscls,  !
     >     gscl3,  !
     >     gscl4   !
      real vmax, vmaxub, vmaxuc, vmaxls, vmaxl3, vmaxl4
      real rdarkub, rdarkuc, rdarkls, rdarkl3, rdarkl4
      real agub, aguc, agls, agl3, agl4
      real anub, anuc, anls, anl3, anl4
      real duma, dumb, dumc, dume, dumq, dump
      real pxaiu, plaiu, pxail, plail
      real cscub, cscuc, cscls, cscl3, cscl4
      real extpar, scale
c
      real kc,     ! co2 kinetic parameter (mol/mol)
     >     ko,     ! o2  kinetic parameter (mol/mol)
*     >     kc15,
*     >     ko15,   ! o2  kinetic parameter (mol/mol) at 15 degrees C
     >     kco2,   ! initial c4 co2 efficiency (mol-co2/m**2/s)
     >     je,     ! 'light limited' rate of photosynthesis (mol-co2/m**2/s)
     >     jc,     ! 'rubisco limited' rate of photosynthesis (mol-co2/m**2/s)
     >     ji,     ! 'co2 limited' rate of photosynthesis (mol-co2/m**2/s)
     >     jp,     ! model-intermediate rate of photosynthesis (mol-co2/m**2/s)
     >     js,     ! sucrose synthesis limitation
     >     gamstar ! gamma*, the co2 compensation points for c3 plants
c
c model parameters
c
c intrinsic quantum efficiency for c3 and c4 plants (dimensionless)
c
*      real alpha3, alpha4
c
*      data alpha3 /0.060/
*      data alpha4 /0.050/
c
c co2/o2 specificity ratio at 15 degrees C (dimensionless)
c
*      real tau15
c
*      data tau15 /4500.0/     
c
c o2/co2 kinetic parameters (mol/mol)
c
*      real kc15, ko12 
c
*      data kc15 /1.5e-04/ 
*      data ko15 /2.5e-01/ 
c
c leaf respiration coefficients
c
*      real gammaub, gammauc, gammals, gammal3, gammal4
c
*      data gammaub /0.0150/   ! broadleaf trees
*      data gammauc /0.0150/   ! conifer trees
*      data gammals /0.0150/   ! shrubs
*      data gammal3 /0.0150/   ! c3 grasses
*      data gammal4 /0.0300/   ! c4 grasses
c
c 'm' coefficients for stomatal conductance relationship
c
*      real coefmub, coefmuc, coefmls, coefml3, coefml4
c
*      data coefmub /10.0/     ! broadleaf trees
*      data coefmuc / 6.0/     ! conifer trees
*      data coefmls / 9.0/     ! shrubs
*      data coefml3 / 9.0/     ! c3 grasses
*      data coefml4 / 4.0/     ! c4 grasses
c
c 'b' coefficients for stomatal conductance relationship 
c (minimum conductance when net photosynthesis is zero)
c
*      real coefbub, coefbuc, coefbls, coefbl3, coefbl4
c
*      data coefbub /0.010/    ! broadleaf trees
*      data coefbuc /0.010/    ! conifer trees
*      data coefbls /0.010/    ! shrubs
*      data coefbl3 /0.010/    ! c3 grasses
*      data coefbl4 /0.040/    ! c4 grasses
c
c absolute minimum stomatal conductances
c
*      real gsubmin, gsucmin, gslsmin, gsl3min, gsl4min
c
*      data gsubmin /0.00001/  ! broadleaf trees
*      data gsucmin /0.00001/  ! conifer trees
*      data gslsmin /0.00001/  ! shrubs
*      data gsl3min /0.00001/  ! c3 grasses
*      data gsl4min /0.00001/  ! c4 grasses
c
c photosynthesis coupling coefficients (dimensionless)
c
*      real theta3
c
*      data theta3 /0.970/     ! c3 photosynthesis
c
*      real theta4, beta4
c
*      data theta4 /0.970/     ! c4 photosynthesis
*      data beta4  /0.800/     ! c4 photosynthesis
c
c maximum values for ci (for model stability)
c
*      real cimax
c
*      data cimax /2000.e-06/  ! maximum values for ci
c
c include water vapor functions
c
      include 'comsat.h'
c
c ---------------------------------------------------------------------
c * * * upper canopy physiology calculations * * *
c ---------------------------------------------------------------------
c
      do 100 i = 1, npoi
c
c calculate physiological parameter values which are a function of temperature
c
        rwork = 3.47e-03 - 1. / tu(i)
c
        tau = tau15 * exp(-4500.0 * rwork)
        kc  = kc15  * exp( 6000.0 * rwork)
        ko  = ko15  * exp( 1500.0 * rwork)
c
        tleaf = tu(i) - 273.16
c
        tempvm = exp(3500.0 * rwork ) /
     >           ((1.0 + exp(0.40 * (  5.0 - tleaf))) * 
     >            (1.0 + exp(0.40 * (tleaf - 50.0))))
c
c upper canopy gamma-star values (mol/mol)
c
        gamstar = o2conc / (2. * tau)
c
c constrain ci values to acceptable bounds -- to help ensure numerical stability
c
        ciub(i) = max (1.05 * gamstar, min (cimax, ciub(i)))
        ciuc(i) = max (1.05 * gamstar, min (cimax, ciuc(i)))
c
c calculate boundary layer parameters (mol/m**2/s) = su / 0.029 * 1.35
c
        gbco2u = min (10.0, max (0.1, su(i) * 25.5))
c 
c calculate the relative humidity in the canopy air space
c with a minimum value of 0.30 to avoid errors in the 
c physiological calculations
c
        esat12 = esat (t12(i))
        qsat12 = qsat (esat12, psurf(i))
        rh12   = max (0.30, q12(i) / qsat12)
c
c ---------------------------------------------------------------------
c broadleaf (evergreen & deciduous) tree physiology 
c ---------------------------------------------------------------------
c 
c nominal values for vmax of top leaf at 15 C (mol-co2/m**2/s)
c
c tropical broadleaf trees          60.0 e-06 mol/m**2/sec
c warm-temperate broadleaf trees    40.0 e-06 mol/m**2/sec
c temperate broadleaf trees         25.0 e-06 mol/m**2/sec
c boreal broadleaf trees            25.0 e-06 mol/m**2/sec
c
*        if (exist(i,1).gt.0.5) then
*          vmaxub = 65.0e-06
*        else if (exist(i,3).gt.0.5) then
*          vmaxub = 40.0e-06
*        else 
*          vmaxub = 30.0e-06
*        endif
*
**** DTP 2001/06/06: Following code replaces above, making initialization
*                    dependent upon parameter values read in from external
*                    canopy parameter file "params.can".
*
        if (exist(i,1).gt.0.5) then
          vmaxub = vmax_pft(1) ! 65.0e-06 ! Tropical broadleaf evergreen
        else if (exist(i,3).gt.0.5) then
          vmaxub = vmax_pft(3) ! 40.0e-06 ! Warm-temperate broadleaf evergreen
        else 
          vmaxub = vmax_pft(5) ! 30.0e-06 ! Temperate or boreal broadleaf cold deciduous
        endif
c
c vmax and dark respiration for current conditions
c
        vmax  = vmaxub * tempvm * stresstu(i)
        rdarkub = gammaub * vmaxub * tempvm
c
c 'light limited' rate of photosynthesis (mol/m**2/s)
c
        je = topparu(i) * 4.59e-06 * alpha3 * (ciub(i) - gamstar) / 
     >       (ciub(i) + 2. * gamstar)
c
c 'rubisco limited' rate of photosynthesis (mol/m**2/s)
c
        jc = vmax * (ciub(i) - gamstar) / 
     >       (ciub(i) + kc * (1. + o2conc / ko))
c
c solution to quadratic equation
c
        duma = theta3
        dumb = je + jc
        dumc = je * jc
c
        dume = max (dumb**2 - 4. * duma * dumc, 0.)
        dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15
c
c calculate the intermediate photosynthesis rate (mol/m**2/s)
c
        jp = min (dumq/duma, dumc/dumq)
c
c 'sucrose synthesis limited' rate of photosynthesis (mol/m**2/s)
c
        js = vmax / 2.2
c
c solution to quadratic equation
c
        duma = beta3
        dumb = jp + js
        dumc = jp * js
c
        dume = max (dumb**2 - 4. * duma * dumc, 0.)
        dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15
c
c calculate the net photosynthesis rate (mol/m**2/s)
c
        agub = min (dumq/duma, dumc/dumq)
        anub = agub - rdarkub
c
c calculate co2 concentrations and stomatal condutance values
c using simple iterative procedure
c
c weight results with the previous iteration's values -- this
c improves convergence by avoiding flip-flop between diffusion
c into and out of the stomatal cavities
c
c calculate new value of cs using implicit scheme
c
        csub(i) = 0.5 * (csub(i) + co2conc - anub / gbco2u)
        csub(i) = max (1.05 * gamstar, csub(i))
c
c calculate new value of gs using implicit scheme
c
        gsub(i) = 0.5 * (gsub(i)  +  (coefmub * anub * rh12 / csub(i) + 
     >                                coefbub * stresstu(i)))
c
        gsub(i) = max (gsubmin, coefbub * stresstu(i), gsub(i))
c
c calculate new value of ci using implicit scheme
c
        ciub(i) = 0.5 * (ciub(i) + csub(i) - 1.6 * anub / gsub(i))
        ciub(i) = max (1.05 * gamstar, min (cimax, ciub(i)))
c
c ---------------------------------------------------------------------
c conifer tree physiology 
c ---------------------------------------------------------------------
c 
c nominal values for vmax of top leaf at 15 C (mol-co2/m**2/s)
c
c temperate conifer trees           30.0 e-06 mol/m**2/sec
c boreal conifer trees              20.0 e-06 mol/m**2/sec
c
*        if (exist(i,4).gt.0.5) then
*          vmaxuc = 30.0e-06
*        else 
*          vmaxuc = 25.0e-06
*        endif
*
**** DTP 2001/06/06: Following code replaces above, making initialization
*                    dependent upon parameter values read in from external
*                    canopy parameter file "params.can".
*
        if (exist(i,4).gt.0.5) then
          vmaxuc = vmax_pft(4) ! 30.0e-06 ! Temperate conifer
        else 
          vmaxuc = vmax_pft(6) ! 25.0e-06 ! Boreal conifer evergreen
        endif

c
c vmax and dark respiration for current conditions
c
        vmax  = vmaxuc * tempvm * stresstu(i)
        rdarkuc = gammauc * vmaxuc * tempvm
c
c 'light limited' rate of photosynthesis (mol/m**2/s)
c
        je = topparu(i) * 4.59e-06 * alpha3 * (ciuc(i) - gamstar) / 
     >       (ciuc(i) + 2. * gamstar)
c
c 'rubisco limited' rate of photosynthesis (mol/m**2/s)
c
        jc = vmax * (ciuc(i) - gamstar) / 
     >       (ciuc(i) + kc * (1. + o2conc / ko))
c
c solution to quadratic equation
c
        duma = theta3
        dumb = je + jc
        dumc = je * jc
c
        dume = max (dumb**2 - 4. * duma * dumc, 0.)
        dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15
c
c calculate the intermediate photosynthesis rate (mol/m**2/s)
c
        jp = min (dumq/duma, dumc/dumq)
c
c 'sucrose synthesis limited' rate of photosynthesis (mol/m**2/s)
c
        js = vmax / 2.2
c
c solution to quadratic equation
c
        duma = beta3
        dumb = jp + js
        dumc = jp * js
c
        dume = max (dumb**2 - 4. * duma * dumc, 0.)
        dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15
c
c calculate the net photosynthesis rate (mol/m**2/s)
c
        aguc = min (dumq/duma, dumc/dumq) 
        anuc = aguc - rdarkuc
c
c calculate co2 concentrations and stomatal condutance values
c using simple iterative procedure
c
c weight results with the previous iteration's values -- this
c improves convergence by avoiding flip-flop between diffusion
c into and out of the stomatal cavities
c
c calculate new value of cs using implicit scheme
c
        csuc(i) = 0.5 * (csuc(i) + co2conc - anuc / gbco2u)
        csuc(i) = max (1.05 * gamstar, csuc(i))
c
c calculate new value of gs using implicit scheme
c
        gsuc(i) = 0.5 * (gsuc(i)  +  (coefmuc * anuc * rh12 / csuc(i) + 
     >                                coefbuc * stresstu(i)))
c
        gsuc(i) = max (gsucmin, coefbuc * stresstu(i), gsuc(i))
c
c calculate new value of ci using implicit scheme
c
        ciuc(i) = 0.5 * (ciuc(i) + csuc(i) - 1.6 * anuc / gsuc(i))
        ciuc(i) = max (1.05 * gamstar, min (cimax, ciuc(i)))
c
c ---------------------------------------------------------------------
c upper canopy scaling
c ---------------------------------------------------------------------
c
c the canopy scaling algorithm assumes that the net photosynthesis
c is proportional to absored par (apar) during the daytime. during night,
c the respiration is scaled using a 10-day running-average daytime canopy
c scaling parameter.
c
c apar(x) = A exp(-k x) + B exp(-h x) + C exp(h x)
c an(x) is proportional to apar(x)
c
c therefore, an(x) = an(0) * apar(x) / apar(0)
c an(x) = an(0) * (A exp(-k x) + B exp(-h x) + C exp(h x)) / 
c                 (A + B + C)
c
c this equation is further simplified to
c an(x) = an(0) * exp (-extpar * x)
c
c an(0) is calculated for a sunlit leaf at the top of the canopy using
c the full-blown plant physiology model (Farquhar/Ball&Berry, Collatz).
c then the approximate par extinction coefficient (extpar) is calculated
c using parameters obtained from the two-stream radiation calculation.
c
c an,canopy avg.= integral (an(x), from 0 to xai) / lai
c               = an(0) * (1 - exp (-extpar * xai )) / (extpar * lai)
c
c the term '(1 - exp (-extpar * xai )) / lai)' scales photosynthesis from leaf
c to canopy level (canopy average) at day time. A 10-day running mean of this
c scaling parameter (weighted by light) is then used to scale the respiration
c during night time.
c
c once canopy average photosynthesis is calculated, then the canopy average
c stomatal conductance is calculated using the 'big leaf approach',i.e. 
c assuming that the canopy is a big leaf and applying the leaf-level stomatal
c conductance equations to the whole canopy.
c
c calculate the approximate par extinction coefficient:
c
c extpar = (k * A + h * B - h * C) / (A + B + C)
c
        extpar = (termu(i,6) * scalcoefu(i,1) +
     >            termu(i,7) * scalcoefu(i,2) -
     >            termu(i,7) * scalcoefu(i,3)) /
     >            max (scalcoefu(i,4), epsilon)
c
        extpar = max (1.e-1, min (1.e+1, extpar))
c
c calculate canopy average photosynthesis (per unit leaf area):
c
        pxaiu = extpar * (lai(i,2) + sai(i,2))
        plaiu = extpar *  lai(i,2)
c
c scale is the parameter that scales from leaf-level photosynthesis to
c canopy average photosynthesis
c CD : replaced 24 (hours) by 86400/dtime for use with other timestep
c
         zweight = exp(-1. / (10.0 * 86400. / dtime))
c
c for non-zero lai
c
        if (plaiu.gt.0.0) then
c
c day-time conditions, use current scaling coefficient
c
          if (topparu(i).gt.10.) then
c
            scale = (1. - exp(-pxaiu)) / plaiu
c
c update 10-day running mean of scale, weighted by light levels
c
            a10scalparamu(i) = zweight * a10scalparamu(i) + 
     >                         (1. - zweight) * scale * topparu(i)
c
            a10daylightu(i)  = zweight * a10daylightu(i) + 
     >                         (1. - zweight) * topparu(i)
c
c night-time conditions, use long-term day-time average scaling coefficient
c
          else
c
            scale = a10scalparamu(i) / a10daylightu(i)
c
          endif
c
c if no lai present
c
        else
c
          scale = 0.0
c
        endif
c
c perform scaling on all carbon fluxes from upper canopy
c
        agcub(i) = agub * scale
        agcuc(i) = aguc * scale
c
        ancub(i) = anub * scale
        ancuc(i) = anuc * scale
c
c calculate diagnostic canopy average surface co2 concentration 
c (big leaf approach)
c
        cscub = max (1.05 * gamstar, co2conc - ancub(i) / gbco2u)
        cscuc = max (1.05 * gamstar, co2conc - ancuc(i) / gbco2u)
c
c calculate diagnostic canopy average stomatal conductance (big leaf approach)
c
        gscub = coefmub * ancub(i) * rh12 / cscub +
     >          coefbub * stresstu(i)
c
        gscuc = coefmuc * ancuc(i) * rh12 / cscuc +
     >          coefbuc * stresstu(i)
c
        gscub = max (gsubmin, coefbub * stresstu(i), gscub)
        gscuc = max (gsucmin, coefbuc * stresstu(i), gscuc)
c
c calculate total canopy and boundary-layer total conductance for 
c water vapor diffusion
c
        rwork = 1. / su(i)
        dump  = 1. / 0.029
c
        totcondub(i) = 1. / (rwork + dump / gscub)
        totconduc(i) = 1. / (rwork + dump / gscuc)
c
c multiply canopy photosynthesis by wet fraction - this calculation is
c done here and not earlier to avoid using within canopy conductance
c
        rwork = 1 - fwetu(i)
c
        agcub(i) = rwork * agcub(i)
        agcuc(i) = rwork * agcuc(i)
c
        ancub(i) = rwork * ancub(i)
        ancuc(i) = rwork * ancuc(i)
c
 100  continue
c
c ---------------------------------------------------------------------
c * * * lower canopy physiology calculations * * *
c ---------------------------------------------------------------------
c
      do 200 i = 1, npoi
c
c calculate physiological parameter values which are a function of temperature
c
        rwork = 3.47e-03 - 1. / tl(i)
c
        tau = tau15 * exp(-4500.0 * rwork)
        kc  = kc15  * exp( 6000.0 * rwork)
        ko  = ko15  * exp( 1500.0 * rwork)
c
        tleaf = tl(i) - 273.16
c
        tempvm = exp(3500.0 * rwork ) /
     >           ((1.0 + exp(0.40 * (  5.0 - tleaf))) * 
     >            (1.0 + exp(0.40 * (tleaf - 50.0))))
c
c lower canopy gamma-star values (mol/mol)
c
        gamstar = o2conc / (2. * tau)
c
c constrain ci values to acceptable bounds -- to help ensure numerical stability
c
        cils(i) = max (1.05 * gamstar, min (cimax, cils(i)))
        cil3(i) = max (1.05 * gamstar, min (cimax, cil3(i)))
        cil4(i) = max (0.0           , min (cimax, cil4(i)))
c
c calculate boundary layer parameters (mol/m**2/s) = su / 0.029 * 1.35
c
        gbco2l = min (10.0, max (0.1, sl(i) * 25.5))
c 
c calculate the relative humidity in the canopy air space
c with a minimum value of 0.30 to avoid errors in the 
c physiological calculations
c
        esat34 = esat (t34(i))
        qsat34 = qsat (esat34, psurf(i))
        rh34   = max (0.30, q34(i) / qsat34)
c
c ---------------------------------------------------------------------
c shrub physiology
c ---------------------------------------------------------------------
c 
c nominal values for vmax of top leaf at 15 C (mol-co2/m**2/s)
c
*        vmaxls = 27.5e-06
*
**** DTP 2001/06/06: Following code replaces above, making initialization
*                    dependent upon parameter values read in from external
*                    canopy parameter file "params.can".
*
        vmaxls = vmax_pft(9) ! 27.5e-06 ! Shrubs (evergreen or cold deciduous) 
c 
c vmax and dark respiration for current conditions
c
        vmax  = vmaxls * tempvm * stresstl(i)
        rdarkls = gammals * vmaxls * tempvm
c
c 'light limited' rate of photosynthesis (mol/m**2/s)
c
        je = topparl(i) * 4.59e-06 * alpha3 * (cils(i) - gamstar) / 
     >       (cils(i) + 2. * gamstar)
c
c 'rubisco limited' rate of photosynthesis (mol/m**2/s)
c
        jc = vmax * (cils(i) - gamstar) / 
     >       (cils(i) + kc * (1. + o2conc / ko))
c
c solution to quadratic equation
c
        duma = theta3
        dumb = je + jc
        dumc = je * jc
c
        dume = max (dumb**2 - 4. * duma * dumc, 0.)
        dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15
c
c calculate the intermediate photosynthesis rate (mol/m**2/s)
c
        jp = min (dumq/duma, dumc/dumq)
c
c 'sucrose synthesis limited' rate of photosynthesis (mol/m**2/s)
c
        js = vmax / 2.2
c
c solution to quadratic equation
c
        duma = beta3
        dumb = jp + js
        dumc = jp * js
c
        dume = max (dumb**2 - 4. * duma * dumc, 0.)
        dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15
c
c calculate the net photosynthesis rate (mol/m**2/s)
c
        agls = min (dumq/duma, dumc/dumq)
        anls = agls - rdarkls
c
c calculate co2 concentrations and stomatal condutance values
c using simple iterative procedure
c
c weight results with the previous iteration's values -- this
c improves convergence by avoiding flip-flop between diffusion
c into and out of the stomatal cavities
c
c calculate new value of cs using implicit scheme
c
        csls(i) = 0.5 * (csls(i) + co2conc - anls / gbco2l)
        csls(i) = max (1.05 * gamstar, csls(i))
c
c calculate new value of gs using implicit scheme
c
        gsls(i) = 0.5 * (gsls(i) + coefmls * anls * rh34 / csls(i) +
     >                             coefbls * stresstl(i))
c
        gsls(i) = max (gslsmin, coefbls * stresstl(i), gsls(i))
c
c calculate new value of ci using implicit scheme
c
        cils(i) = 0.5 * (cils(i) + csls(i) - 1.6 * anls / gsls(i))
        cils(i) = max (1.05 * gamstar, min (cimax, cils(i)))
c
c ---------------------------------------------------------------------
c c3 grass physiology
c ---------------------------------------------------------------------
c 
c nominal values for vmax of top leaf at 15 C (mol-co2/m**2/s)
c
*       vmaxl3 = 25.0e-06
*
**** DTP 2001/06/06: Following code replaces above, making initialization
*                    dependent upon parameter value read in from external
*                    canopy parameter file "params.can".
*
        vmaxl3 = vmax_pft(12) ! 25.0e-06 ! C3 grasses
c 
c vmax and dark respiration for current conditions
c
        vmax  = vmaxl3 * tempvm * stresstl(i)
        rdarkl3 = gammal3 * vmaxl3 * tempvm
c
c 'light limited' rate of photosynthesis (mol/m**2/s)
c
        je = topparl(i) * 4.59e-06 * alpha3 * (cil3(i) - gamstar) / 
     >       (cil3(i) + 2. * gamstar)
c
c 'rubisco limited' rate of photosynthesis (mol/m**2/s)
c
        jc = vmax * (cil3(i) - gamstar) / 
     >       (cil3(i) + kc * (1. + o2conc / ko))
c
c solution to quadratic equation
c
        duma = theta3
        dumb = je + jc
        dumc = je * jc
c
        dume = max (dumb**2 - 4. * duma * dumc, 0.)
        dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15
c
c calculate the intermediate photosynthesis rate (mol/m**2/s)
c
        jp = min (dumq/duma, dumc/dumq)
c
c 'sucrose synthesis limited' rate of photosynthesis (mol/m**2/s)
c
        js = vmax / 2.2
c
c solution to quadratic equation
c
        duma = beta3
        dumb = jp + js
        dumc = jp * js
c
        dume = max (dumb**2 - 4. * duma * dumc, 0.)
        dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15
c
c calculate the net photosynthesis rate (mol/m**2/s)
c
        agl3 = min (dumq/duma, dumc/dumq)
        anl3 = agl3 - rdarkl3
c
c calculate co2 concentrations and stomatal condutance values
c using simple iterative procedure
c
c weight results with the previous iteration's values -- this
c improves convergence by avoiding flip-flop between diffusion
c into and out of the stomatal cavities
c
c calculate new value of cs using implicit scheme
c
        csl3(i) = 0.5 * (csl3(i) + co2conc - anl3 / gbco2l)
        csl3(i) = max (1.05 * gamstar, csl3(i))
c
c calculate new value of gs using implicit scheme
c
        gsl3(i) = 0.5 * (gsl3(i) + coefml3 * anl3 * rh34 / csl3(i) +
     >                   coefbl3 * stresstl(i))
c
        gsl3(i) = max (gsl3min, coefbl3 * stresstl(i), gsl3(i))
c
c calculate new value of ci using implicit scheme
c
        cil3(i) = 0.5 * (cil3(i) + csl3(i) - 1.6 * anl3 / gsl3(i))
        cil3(i) = max (1.05 * gamstar, min (cimax, cil3(i)))
c
c ---------------------------------------------------------------------
c c4 grass physiology
c ---------------------------------------------------------------------
c
c nominal values for vmax of top leaf at 15 C (mol-co2/m**2/s)
c
*       vmaxl4 = 15.0e-06
*
**** DTP 2001/06/06: Following code replaces above, making initialization
*                    dependent upon parameter value read in from external
*                    canopy parameter file "params.can".
*
        vmaxl4 = vmax_pft(11) ! 15.0e-06 ! C4 grasses
c
c calculate the parameter values which are a function of temperature
c
        rwork = 3.47e-03 - 1. / tl(i)
c
        tleaf = tl(i) - 273.16
c
        tempvm = exp(3500.0 * rwork ) /
     >           ((1.0 + exp(0.40 * ( 10.0 - tleaf))) * 
     >            (1.0 + exp(0.40 * (tleaf - 50.0))))
c
c vmax and dark respiration for current conditions
c
        vmax  = vmaxl4 * tempvm * stresstl(i)
        rdarkl4 = gammal4 * vmaxl4 * tempvm
c
c initial c4 co2 efficiency (mol/m**2/s)
c
        kco2 = 18.0e+03 * vmax
c
c 'light limited' rate of photosynthesis (mol/m**2/s)
c
        je = topparl(i) * 4.59e-06 * alpha4
c
c 'rubisco limited' rate of photosynthesis
c
        jc = vmax
c
c solve for intermediate photosynthesis rate
c
        duma = theta4
        dumb = je + jc
        dumc = je * jc
c
        dume = max (dumb**2 - 4. * duma * dumc, 0.)
        dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15
c
        jp = min (dumq/duma, dumc/dumq)
c
c 'carbon dioxide limited' rate of photosynthesis (mol/m**2/s)
c
        ji = kco2 * cil4(i)
c
c solution to quadratic equation
c
        duma = beta4
        dumb = jp + ji
        dumc = jp * ji
c
        dume = max (dumb**2 - 4. * duma * dumc, 0.)
        dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15
c
c calculate the net photosynthesis rate (mol/m**2/s)
c
        agl4 = min (dumq/duma, dumc/dumq)
        anl4 = agl4 - rdarkl4
c
c calculate co2 concentrations and stomatal condutance values
c using simple iterative procedure
c
c weight results with the previous iteration's values -- this
c improves convergence by avoiding flip-flop between diffusion
c into and out of the stomatal cavities
c
c calculate new value of cs using implicit scheme
c CD: For numerical stability (to avoid division by zero in gsl4), 
c csl4 is limited to 1e-8 mol_co2/mol_air.
c  
        csl4(i) = 0.5 * (csl4(i) + co2conc - anl4 / gbco2l)
        csl4(i) = max (1.e-8, csl4(i))
c
c calculate new value of gs using implicit scheme
c
        gsl4(i) = 0.5 * (gsl4(i) + coefml4 * anl4 * rh34 / csl4(i) +
     >                   coefbl4 * stresstl(i))
c
        gsl4(i) = max (gsl4min, coefbl4 * stresstl(i), gsl4(i))
c
c calculate new value of ci using implicit scheme
c
        cil4(i) = 0.5 * (cil4(i) + csl4(i) - 1.6 * anl4 / gsl4(i))
        cil4(i) = max (0.0, min (cimax, cil4(i)))
c
c ---------------------------------------------------------------------
c lower canopy scaling
c ---------------------------------------------------------------------
c
c calculate the approximate extinction coefficient
c
        extpar = (terml(i,6) * scalcoefl(i,1) + 
     >            terml(i,7) * scalcoefl(i,2) -
     >            terml(i,7) * scalcoefl(i,3)) /
     >            max (scalcoefl(i,4), epsilon)
c
        extpar = max (1.e-1, min (1.e+1, extpar))
c
c calculate canopy average photosynthesis (per unit leaf area):
c
        pxail = extpar * (lai(i,1) + sai(i,1))
        plail = extpar *  lai(i,1)
c
c scale is the parameter that scales from leaf-level photosynthesis to
c canopy average photosynthesis
c CD : replaced 24 (hours) by 86400/dtime for use with other timestep
c
        zweight = exp(-1. / (10.0 * 86400. / dtime))
c
c for non-zero lai
c
        if (plail.gt.0.0) then
c
c day-time conditions, use current scaling coefficient
c
          if (topparl(i).gt.10.) then
c
            scale = (1. - exp(-pxail)) / plail
c
c update 10-day running mean of scale, weighted by light levels
c
            a10scalparaml(i) = zweight * a10scalparaml(i) + 
     >                         (1. - zweight) * scale * topparl(i)
c
            a10daylightl(i)  = zweight * a10daylightl(i) + 
     >                         (1. - zweight) * topparl(i)
c
c night-time conditions, use long-term day-time average scaling coefficient
c
          else
c
            scale = a10scalparaml(i) / a10daylightl(i)
c
          endif
c
c if no lai present
c
        else
c
          scale = 0.0
c
        endif
c
c perform scaling on all carbon fluxes from upper canopy
c
        agcls(i) = agls * scale
        agcl4(i) = agl4 * scale
        agcl3(i) = agl3 * scale
c
        ancls(i) = anls * scale
        ancl4(i) = anl4 * scale
        ancl3(i) = anl3 * scale
c
c calculate canopy average surface co2 concentration
c CD: For numerical stability (to avoid division by zero in gscl4),
c cscl4 is limited to 1e-8 mol_co2/mol_air.
c
        cscls = max (1.05 * gamstar, co2conc - ancls(i) / gbco2l)
        cscl3 = max (1.05 * gamstar, co2conc - ancl3(i) / gbco2l)
        cscl4 = max (1.e-8         , co2conc - ancl4(i) / gbco2l)
c
c calculate canopy average stomatal conductance
c
        gscls = coefmls * ancls(i) * rh34 / cscls +
     >          coefbls * stresstl(i)
c
        gscl3 = coefml3 * ancl3(i) * rh34 / cscl3 +
     >          coefbl3 * stresstl(i)
c
        gscl4 = coefml4 * ancl4(i) * rh34 / cscl4 +
     >          coefbl4 * stresstl(i)
c
        gscls = max (gslsmin, coefbls * stresstl(i), gscls)
        gscl3 = max (gsl3min, coefbl3 * stresstl(i), gscl3)
        gscl4 = max (gsl4min, coefbl4 * stresstl(i), gscl4)
c
c calculate canopy and boundary-layer total conductance for water vapor diffusion
c
        rwork = 1. / sl(i)
        dump =  1. / 0.029
c
        totcondls(i) = 1. / (rwork + dump / gscls)
        totcondl3(i) = 1. / (rwork + dump / gscl3)
        totcondl4(i) = 1. / (rwork + dump / gscl4)
c
c multiply canopy photosynthesis by wet fraction -- this calculation is
c done here and not earlier to avoid using within canopy conductance
c
        rwork = 1. - fwetl(i)
c
        agcls(i) = rwork * agcls(i)
        agcl3(i) = rwork * agcl3(i)
        agcl4(i) = rwork * agcl4(i)
c
        ancls(i) = rwork * ancls(i)
        ancl3(i) = rwork * ancl3(i)
        ancl4(i) = rwork * ancl4(i)
c
 200  continue
c
c return to main program
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine drystress
c ---------------------------------------------------------------------
c
c common blocks
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'comsoi.h'
      include 'comveg.h'
c
c local variables
c
      integer i, k    ! loop indices
c
      real stressfac, ! to calculate moisture stress factor 
     >     awc,       ! available water content (fraction)
     >     znorm,     ! normalizing factor
     >     zwilt      ! function of awc, =1 if awc = 1 (no stress)
c
c stressfac determines the 'strength' of the soil moisture
c stress on physiological processes
c
c strictly speaking, stresst* is multiplied to the vmax
c parameters used in the photosynthesis calculations
c
c stressfac determines the shape of the soil moisture response
c
      stressfac = -5.0
c
      znorm = 1.0 - exp(stressfac)
c
      do 100 i = 1, npoi
c
c initialize stress parameter
c
        stresstl(i) = 0.0
        stresstu(i) = 0.0
c
c fraction of soil water uptake in each layer
c
        do 110 k = 1, nsoilay
c
c plant available water content (fraction)
c
          awc = min (1.0, max (0.0,
     >              (wsoi(i,k)*(1 - wisoi(i,k))   - swilt(i,k)) /
     >              (sfield(i,k) - swilt(i,k))
     >              )         )
c
          zwilt = (1. - exp(stressfac * awc)) / znorm
c
c update for each layer
c
          stressl(i,k) = froot(k,1) * max (0.0, min (1.0, zwilt))
          stressu(i,k) = froot(k,2) * max (0.0, min (1.0, zwilt))
c
c integral over rooting profile
c
          stresstl(i) = stresstl(i) + stressl(i,k)
          stresstu(i) = stresstu(i) + stressu(i,k)
c
 110    continue
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
      subroutine co2 (co2init, co2conc, iyear)
c ---------------------------------------------------------------------
c
      include 'implicit.h'
c     include 'comveg.h'
c
c Arguments 
c
      integer iyear    ! current year
c
      real co2init,    ! input atmospheric co2 concentration
     >     co2conc     ! output " for year iyear   
c
c calculate co2 concentration for this year
c
c     if (iyear.lt.1860) then
c
        co2conc = co2init
c
c     else
c
c 1992 IPCC estimates
c
c       iyr = iyear - 1860 + 1
c       co2conc = (297.12 - 0.26716 * iyr +
c    >                      0.0015368 * iyr**2 +
c    >                      3.451e-5 * iyr**3) * 1.e-6
c
c
c M. El Maayar: 1996 IPCC estimates
c
c       iyr = iyear - 1860 + 1
c       co2conc = (303.514 - 0.57881 * iyr +
c    >                      0.00622 * iyr**2 +
c    >                      1.3e-5 * iyr**3) * 1.e-6
c
c
c     end if
c
      return
      end
c
