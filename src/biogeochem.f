c
c #####      #     ####    ####   ######   ####    ####   #    #  ######  #    #
c #    #     #    #    #  #    #  #       #    #  #    #  #    #  #       ##  ##
c #####      #    #    #  #       #####   #    #  #       ######  #####   # ## #
c #    #     #    #    #  #  ###  #       #    #  #       #    #  #       #    #
c #    #     #    #    #  #    #  #       #    #  #    #  #    #  #       #    #
c #####      #     ####    ####   ######   ####    ####   #    #  ######  #    #
c
c
c --------------------------------------------------------------------------
      subroutine soilbgc (iyear,iyear0,imonth,iday,nspinsoil,
     >                    spin,spinmax)
c --------------------------------------------------------------------------
c
      use implicit
c
      use compar
      use comsoi
      use comsum
      use comveg
      use comatm
c
c Arguments (input)
c
      integer   iday,           ! day in month
     >          iyear,          ! current year
     >          iyear0,         ! initial year
     >          imonth,         ! current month
     >          nspinsoil,      ! year when soil carbon spinup stops
     >          spin,           ! # of times soilbgc has been called in the current day
     >          spinmax         ! total # of times soilbgc is called per day (spinup)
c 
c local variables
c
      integer i            ! loop indice
c
      real totts,          ! 1/ndaypy
     >     fracll,         ! lignin fraction of leaves 
     >     fracls,         ! structural fraction of leaves 
     >     fraclm,         ! metabolic fraction of leaves 
     >     fracrl,         ! lignin fraction of roots
     >     fracrs,         ! structural fraction of roots 
     >     fracrm,         ! metabolic fraction of roots
     >     fracwl,         ! lignin fraction of wood
     >     fracws,         ! structural fraction of wood
     >     fracwm          ! metabolic fraction of wood 

      real      outclm(npoi),   ! c leaving leaf metabolic pool 
     >          outcls(npoi),   ! c leaving leaf structural pool 
     >          outcll(npoi),   ! c leaving leaf lignin pool
     >          outcrm(npoi),   ! c leaving root metabolic pool 
     >          outcrs(npoi),   ! c leaving root structural pool 
     >          outcrl(npoi),   ! c leaving root lignin pool
     >          outcwm(npoi),   ! c leaving woody metabolic carbon pool
     >          outcws(npoi),   ! c leaving woody structural carbon pool
     >          outcwl(npoi),   ! c leaving woody lignin carbon pool
     >          outcsb(npoi),   ! flow of passive c to biomass
     >          outcps(npoi),   ! flow of protected om to passive pool 
     >          outcns(npoi),   ! flow of non-protected om to passive pool
     >          outcnb(npoi),   ! flow of non-protected om to biomass 
     >          outcpb(npoi),   ! flow of protected om to biomass
     >          outcbp(npoi),   ! c leaving protected biomass pool  
     >          outcbn(npoi),   ! c leaving non-protected biomass pool
     >          totc(npoi)      ! total c in soil
c
      real      dbdt(npoi),     ! change of c in biomass pools with time 
     >          dcndt(npoi),    ! change of c in non-protected om with time
     >          dcpdt(npoi),    ! change of c in protected om with time
     >          dcsdt(npoi),    ! change of c in passive om with time
     >          totmin(npoi),   ! total nitrogen mineralization
     >          totimm(npoi),   ! total nitrogen immobilization 
     >          netmin(npoi),   ! net nitrogen mineralization
     >          nbiors(npoi),
     >          nbiols(npoi),
     >          nbiows(npoi),
     >          nbiowm(npoi),
     >          nbiolm(npoi),
     >          nbiorm(npoi),
     >          nbioslon(npoi),
     >          nbioslop(npoi),
     >          nbiopas(npoi),
     >          nminrs(npoi),
     >          nminls(npoi),
     >          nminws(npoi),
     >          nminwm(npoi),
     >          nminlm(npoi),
     >          nminrm(npoi),
     >          nminslon(npoi),
     >          nminslop(npoi),
     >          nminpas(npoi),
     >          nrelps(npoi),
     >          nrelns(npoi),
     >          nrelbn(npoi),
     >          nrelbp(npoi),
     >          nrelll(npoi),
     >          nrelrl(npoi),
     >          nrelwl(npoi),
     >          totnrel(npoi),
     >          ymintot(npoi),
     >          yminmic(npoi)
c
c nitrogen in litter and soil pools
c
      real      nlitlm(npoi),
     >          nlitls(npoi),
     >          nlitll(npoi),
     >          nlitrm(npoi),
     >          nlitrs(npoi), 
     >          nlitrl(npoi),
     >          nlitwm(npoi),
     >          nlitws(npoi),
     >          nlitwl(npoi),
     >          nsoislop(npoi),
     >          nsoipas(npoi),
     >          nsoislon(npoi),
     >          budgetn(npoi)
c
c variables controlling constraints on microbial biomass 
c
      real      cmicn(npoi),
     >          cmicp(npoi),
     >          cmicmx(npoi)
c
c variables controlling leaching, calculating co2 respiration and n deposition
c
      real      cleach(npoi),
     >          totcbegin(npoi),
     >          totcend(npoi),
     >          totcin(npoi),
     >          fixsoin(npoi),
     >          deposn(npoi)
c
      real      fleach,
     >          h20
c
c decay constants for c pools
c
      real      klm,          ! leaf metabolic litter 
     >          kls,          ! leaf structural litter
     >          kll,          ! leaf lignin
     >          krm,          ! root metabolic litter
     >          krs,          ! root structural litter
     >          krl,          ! root lignin
     >          kwm,          ! woody metabolic litter
     >          kws,          ! woody structural litter
     >          kwl,          ! wood  lignin
     >          kbn,          ! microbial biomass --> nonprotected om 
     >          kbp,          ! microbial biomass --> protected om
     >          knb,          ! nonprotected om   --> biomass
     >          kns,          ! nonprotected om   --> passive c 
     >          kpb,          ! protected om      --> biomass
     >          kps,          ! protected om      --> passive c
     >          ksb           ! passive c         --> biomass
c
c efficiencies for microbial decomposition
c
      real      ylm,          ! leaf metabolic litter decomposition 
     >          yls,          ! leaf structural litter decomposition
     >          yll,          ! leaf lignin
     >          yrm,          ! root metabolic litter decomposition
     >          yrs,          ! root structural litter decomposition
     >          yrl,          ! root lignin
     >          ywm,          ! woody metabolic litter decomposition
     >          yws,          ! woody structural litter decomposition
     >          ywl,          ! wood lignin
     >          ybn,          ! microbial biomass to nonprotected om
     >          ybp,          ! microbial biomass to protected om
     >          ynb,          ! nonprotected om to biomass
     >          yns,          ! nonprotected om to passive  c
     >          ypb,          ! protected om to biomass
     >          yps,          ! protected om to passive c
     >          ysb           ! passive c to biomass
c
      real      cnr(10)       ! c:n ratios of c and litter pools
c
c constants for calculating fraction of litterall in structural
c metabolic and lignified (resistant) fractions
c
      real      cnleaf,     ! input c:n ratio of leaf litterfall 
     >          cnwood,     ! input c:n ratio of wood litter
     >          cnroot,     ! input c:n ratio of root litter turnover
     >          rconst,     ! value set to 1200.  from Verberne model 
     >          fmax        ! maximum fraction allowed in metabolic pool
c 
c variables added to do daily time series of some values
c
      integer   gridpt,
     >          kk
c
c variables dealing with soil texture and algorithms
c
      integer   msand,
     >          mclay
c
      real      fsand,
     >          fclay,
     >          cfrac,
     >          texfact,
     >          fbpom,
     >          fbsom,
     >          rdepth,
     >          effac,
     >          lig_frac
c
       gridpt = npoi		   ! total number of gridpoints used
c
c total timesteps (daily) used to divide litterfall into daily fractions 
c
       totts=1./float(ndaypy)
c
c -------------------------------------------------------------------------------------
c specific maximum decay rate or growth constants; rates are per day
c constants are taken from Parton et al., 1987 and Verberne et al., 1990
c and special issue of Geoderma (comparison of 9 organic matter models) in Dec. 1997
c
c leaching parameterization was changed to agree with field data,
c this caused a changing of the below constants.  
c
c approximate factors for Verberne et al. model where efficiencies are 100%
c for some of the transformations: one problem was that their rate constants were
c based on 25C, and our modifying functions are based on 15 C...thus the rate constants
c are somewhat smaller compared to the Verberne et al. (1990) model parameters
c rates are based on a daily decomposition timestep (per day)
c -------------------------------------------------------------------------------------
c
c leaf litter constants
c
      klm = 0.15 		!dpm leaf --> microbial biomass
      kls = 0.01 		!spm leaf --> microbial biomass
      kll = 0.01		!rpm leaf --> non or protected om
c
c root litter constants
c
      krm = 0.10		!dpm root --> microbial biomass
      krs = 0.005 		!spm root --> microbial biomass
      krl = 0.005		!rpm root --> non or protected om 
c
c woody litter constants
c
      kwm = 0.001		!dpm wood --> microbial biomass
      kws = 0.001	 	!spm wood --> microbial biomass
      kwl = 0.001		!rpm wood --> non or protected om 
c
c biomass constants
c
      kbn = 0.045		!biomass --> non protected organic matter 
      kbp = 0.005		!biomass --> protected organic matter
c
c slow and passive c pools
c
      knb = 0.001		!non protected om --> biomass
      kns = 0.000001		!non protected om --> stablized om
      kpb = 0.0001 		!protected om     --> biomass
      kps = 0.000001		!protected om     --> stablized om
      ksb = 8.0e-07		!stablized om     --> biomass
c
c ---------------------------------------------------------------------
c  yield (efficiency) with which microbes gain biomass from c source
c  the rest is driven off as co2 respiration (microbial respiration)
c  all of the respiration produced by microbes is assumed to leave
c  the soil profile over the course of a year
c  taken primarily from the models of Verberne and CENTURY
c ---------------------------------------------------------------------
c
      ylm = 0.4       ! metabolic material efficiencies
      yrm = 0.4
      ywm = 0.4
      yls = 0.3       ! structural efficiencies
      yrs = 0.3
      yws = 0.3
c
      yll = 1.0       ! resistant fraction
      yrl = 1.0 
      ywl = 1.0 
      ybn = 1.0       ! biomass       --> non-protected pool
      ybp = 1.0       ! biomass       --> protected pool
      yps = 1.0       ! protected     --> passive
      yns = 1.0       ! non-protected --> passive
c
      ysb = 0.20       ! passive pool  --> biomass
      ypb = 0.20       ! protected     --> biomass
      ynb = 0.25       ! non-protected --> biomass
c
c -------------------------------------------------------------------
c split of lignified litter material between protected/non-protected
c slow OM pools
c -------------------------------------------------------------------
c
      lig_frac = 0.50 
c
c -------------------------------------------------------------------
c protected biomass as a fraction of total soil organic carbon
c from Verberne et al., 1990
c -------------------------------------------------------------------
c
      fbsom = 0.017
c
c ---------------------------------------------------------------------
c (effac) --> efficiency of microbial biomass reincorporated
c into biomass pool.(from NCSOIL parameterizations; Molina et al., 1983)
c ---------------------------------------------------------------------
c
       effac = 0.40 
c
c ---------------------------------------------------------------------
c define C:N ratios of substrate pools and biomass
c metabolic, structural, and lignin are for Leaves and roots
c values from Parton et al., 1987 and Whitmore and Parry, 1988
c index: 1 - biomass, 2 - passive pool, 3- slow protected c,
c 4 - slow carbon, non-protected, 5 - resistant, 6 - structural plant
c leaf and root litter, 7 - metabolic plant and root litter, 
c 8- woody biomass
c ---------------------------------------------------------------------
c
       cnr(1)  = 8.0       !c:n ratio of microbial biomass
       cnr(2)  = 15.0      !c:n ratio of passive soil carbon
       cnr(3)  = 10.0      !c:n ratio of protected slow soil carbon
       cnr(4)  = 15.0      !c:n ratio of non-protected slow soil C
       cnr(5)  = 100.0     !c:n ratio of resistant litter lignin
       cnr(6)  = 150.0     !c:n ratio of structural plant litter
       cnr(7)  = 6.0	   !c:n ratio of metabolic plant litter
       cnr(8)  = 250.0     !c:n Ratio of woody components
c
c ---------------------------------------------------------------------
c calculate the fraction of wood, roots and leaves that are structural,
c decomposable, and resistant based on equations presented in Verberne
c model discussion (Geoderma, December 1997 special issue).  fmax is the
c maximum fraction allowed in resistant fraction, rconst is a constant
c defined as 1200.  The cnratio of each plant part has to be less than
c the value of structural defined above (i.e. 150) otherwise the equations
c are unstable...thus the wood litter pool value for cnr(6) is substituted
c with a value higher than that for cnwood (i.e. 250).  this is 
c insignificant for wood since 97% is structural anyways.
c
c ** NOTE ******** 
c Would like to incorporate different C:N ratios of residue/roots for
c different biome types based on literature search
c average c:n ratio would be based on litter inputs from each pft
c ****************
c ---------------------------------------------------------------------
c
c equations were changed on 1-26-99 for erratum in literature (Whitmore
c et al. 1997) which had an error in equations to split litterfall into
c the correct three fractions
c 
c
       fmax   = 0.45
       rconst = 1200.0
       cnleaf = 40.0      ! average c:n ratio for leaf litterfall
       cnroot = 60.0      ! average c:n ratio for root turnover
       cnwood = 200.0     ! average c:n ratio for woody debris
c
c leaf litter 
c
       fracll = fmax * (cnleaf**2)/(rconst + cnleaf**2)
       fracls = (1./cnleaf - fracll/cnr(5) - (1.-fracll)/cnr(7))/
     >          (1./cnr(6) - 1./cnr(7))
       fraclm = 1.0 - fracll - fracls
c
c root litter
c
       fracrl = fmax * (cnroot**2)/(rconst + cnroot**2)
       fracrs = (1./cnroot - fracrl/cnr(5) - (1.-fracrl)/cnr(7))/
     >          (1./cnr(6) - 1./cnr(7))
       fracrm = 1.0 - fracrl - fracrs
c
c wood litter
c
       fracwl = fmax * (cnwood**2)/(rconst + cnwood**2)
       fracws = (1./cnwood - fracwl/cnr(5) - (1.-fracwl)/cnr(7))/
     >          (1./cnr(8) - 1./cnr(7))
       fracwm = 1.0 - fracwl - fracws
c
c
      do 100 i = 1, npoi
c
c ---------------------------------------------------------------------
c fraction of decomposing microbial biomass into protected organic
c matter; taken from the model of Verberne et al., 1990
c this is the proportion of decomposing dead microbial biomass that
c is transferred to a protected pool vs. a non-protected pool
c related to the clay content of the soil. in sandy soils, fbpom = 0.3,
c whereas in clay soils fbpom = 0.7.  created a linear function based
c on clay fraction of soil to adjust according to amount of clay in
c the top 1 m of profile (weighted average according to depth of each
c layer)
c
c also take care of calculation of texfact, which is a leaching
c parameter based on the average sand fraction of the top 1 m of
c soil
c ---------------------------------------------------------------------
c
        rdepth   = 1./(hsoi(1) + hsoi(2) + hsoi(3) + hsoi(4))
        cfrac    = 0.0
        texfact  = 0.0 
c
        do 90 kk = 1, 4                    ! top 1 m of soil -- 4 layers
c
          msand    = nint(sand(i,kk)) 
          mclay    = nint(clay(i,kk)) 
          fclay    = 0.01 * mclay
          fsand    = 0.01 * msand 
          cfrac    = cfrac   + fclay * hsoi(kk)
          texfact  = texfact + fsand * hsoi(kk)
c
 90     continue
c
        cfrac   = cfrac   * rdepth
        texfact = texfact * rdepth
c
c if cfrac is greater than 0.4, set fbpom = 0.7, if cfrac is less
c than 0.17, set fbpom = 0.30 (sandy soil)
c
c        fbpom = min(max(0.3, cfrac/0.4 * 0.7),0.7)      
         fbpom = 0.50
c
c ------------------------------------------------------------------------
c total soil carbon initialized to 0 at beginning of model run
c used in calculation of soil co2 respiration from microbial decomposition 
c ------------------------------------------------------------------------
c
       if (iday .eq. 1 .and. imonth .eq. 1 .and. iyear .eq. iyear0) then
         totcbegin(i) = 0.0
         storedn(i)   = 0.0 
       endif
c
c ------------------------------------------------------------------------
c initialize yearly summation of net mineralization and co2 respiration
c to 0 at beginning of each year; because these quantities are usually 
c reported on a yearly basis, we wish to do the same in the model so we
c can compare easily with the data.
c ------------------------------------------------------------------------
c
       if (iday .eq. 1 .and. imonth .eq. 1) then
         yrleach(i) = 0.0
         cleach(i)  = 0.0
         ynleach(i) = 0.0
         ymintot(i) = 0.0
         yminmic(i) = 0.0
       endif
c       
c determine amount of substrate available to microbial growth
c
c calculate the total amount of litterfall entering soil(C)
c
       totcin(i) =  falll(i)*totts + fallr(i)*totts
     >            + fallw(i)*totts
c
c
c calculate the current total amount of carbon at each grid cell
c
       totc(i) = clitlm(i) + clitls(i) + clitrm(i) + clitrs(i) +
     >       clitwm(i) + clitws(i) + csoislop(i) + csoislon(i) +
     > csoipas(i) + totcmic(i) + clitll(i) + clitrl(i) + clitwl(i)
c
c beginning amount of soil C at each timestep (used for respiration
c calculation)
c
       totcbegin(i) = totc(i)
c
c ------------------------------------------------------------------------
c split current amount of total soil microbes
c maximum amount of biomass is a function of the total soil C
c from Verberne et al., 1990
c ------------------------------------------------------------------------
c
c      totcmic(i) = cmicp(i) + cmicn(i)
       cmicmx(i) = fbsom * totc(i) 
c
c calculate the amount of protected and unprotected biomass
c
       if (totcmic(i) .ge. cmicmx(i)) then
c
         cmicp(i) = cmicmx(i)
         cmicn(i) = totcmic(i) - cmicmx(i)
c
       else
c
         cmicn(i) = 0.0
         cmicp(i) = totcmic(i)
c
       endif
c
c ---------------------------------------------------------------
c litter pools 
c
c add in the amount of litterfall, and root turnover
c ---------------------------------------------------------------
c
       clitlm(i) = clitlm(i) + (fraclm * falll(i)*totts)  
       clitls(i) = clitls(i) + (fracls * falll(i)*totts)  
       clitll(i) = clitll(i) + (fracll * falll(i)*totts)  
       clitrm(i) = clitrm(i) + (fracrm * fallr(i)*totts)  
       clitrs(i) = clitrs(i) + (fracrs * fallr(i)*totts)  
       clitrl(i) = clitrl(i) + (fracrl * fallr(i)*totts)  
       clitwm(i) = clitwm(i) + (fracwm * fallw(i)*totts)  
       clitws(i) = clitws(i) + (fracws * fallw(i)*totts)  
       clitwl(i) = clitwl(i) + (fracwl * fallw(i)*totts)  
c
c ---------------------------------------------------------------
c calculate microbial growth rates based on available C sources
c to microbes (substrate : litter, C in slow, passive pools)
c the amount of biomass added cannot be larger than the amount of
c available carbon from substrates and other pools at this point.
c ---------------------------------------------------------------
c
c
       outcrs(i) = min(decomps(i) * krs * clitrs(i),clitrs(i))
       outcws(i) = min(decompl(i) * kws * clitws(i),clitws(i))
       outcls(i) = min(decompl(i) * kls * clitls(i),clitls(i))
       outclm(i) = min(decompl(i) * klm * clitlm(i),clitlm(i))
       outcrm(i) = min(decomps(i) * krm * clitrm(i),clitrm(i))
       outcwm(i) = min(decompl(i) * kwm * clitwm(i),clitwm(i))
       outcnb(i) = min(decomps(i) * knb * csoislon(i),
     >             csoislon(i))
c
       outcpb(i) = min(decomps(i) * kpb * csoislop(i),
     >             csoislop(i))
c
       outcsb(i) = min(decomps(i) * ksb * csoipas(i),
     >             csoipas(i))
c
c ---------------------------------------------------------------
c calculate turnover of microbial biomass
c two disctinct pools: one with rapid turnover, and one with slow
c turnover rate
c ---------------------------------------------------------------
c
       outcbp(i) = min(kbp * cmicp(i),cmicp(i))
       outcbn(i) = min(kbn * cmicn(i),cmicn(i))
c
c ---------------------------------------------------------------------
c recycle microbes back to respective microbial pools based on effac as
c discussed in NCSOIL model from Molina et al., 1983
c ---------------------------------------------------------------------
c
       outcbp(i) = outcbp(i) *  effac
       outcbn(i) = outcbn(i) *  effac
c
c -------------------------------------------------------------------------
c have to adjust inputs into microbial pool for the slow
c and passive carbon amounts that are leaving their respective
c pools at an increased rate during the spinup procedure.
c these values should be decreased by the respective spinup factors
c because the microbial pools will otherwise become larger without
c scientific reason due to the spinup relationships used.
c 3 main pools: outcpb, outcnb, outcsb
c -------------------------------------------------------------------------
c
       dbdt(i) =  outcrs(i) * yrs + outcws(i) * yws +
     >            outcls(i) * yls + outclm(i) * ylm +
     >            outcrm(i) * yrm + outcwm(i) * ywm +
     >            outcnb(i) * ynb + 
     >            outcpb(i) * ypb +
     >            outcsb(i) * ysb - outcbp(i) -
     >            outcbn(i)
c 
c -------------------------------------------------------------------------
c change in non-protected organic matter from growth in microbial
c biomass, lignin input, and stablized organic matter pool
c the flow out of the pool from its decomposition is always less
c the yield--which is factored into the pool it is flowing into
c -------------------------------------------------------------------------
c
       outcll(i) = min(decompl(i) * kll * clitll(i),clitll(i))
       outcrl(i) = min(decomps(i) * krl * clitrl(i),clitrl(i))
       outcwl(i) = min(decompl(i) * kwl * clitwl(i),clitwl(i))
       outcns(i) = min(decomps(i) * kns * csoislon(i),
     >                 csoislon(i))
c
c ------------------------------------------------------------ 
c the lig_frac  factor only applies to lignin content...half goes to
c protected slow OM, and half goes to non protected slow OM
c ------------------------------------------------------------
c
       dcndt(i) =  (lig_frac * (outcll(i) * yll + outcrl(i) * yrl +
     >             outcwl(i) * ywl) +
     >             (1. - fbpom) * (ybn * outcbn(i) +
     >             ybp * outcbp(i))) - outcnb(i) - outcns(i)

c
c ------------------------------------------------------------
c change in protected organic matter from growth in microbial 
c biomass, lignin input, and stablized organic matter pool
c ------------------------------------------------------------
c
       outcps(i) = min(decomps(i) * kps * csoislop(i),
     >                 csoislop(i))
c
c ------------------------------------------------------------
c the lig_frac factor only applies to lignin content...half goes to
c protected slow OM, and half goes to non protected slow OM
c ------------------------------------------------------------
c
       dcpdt(i) = (lig_frac * (outcll(i)*yll+outcrl(i)*yrl +
     >            outcwl(i) * ywl) +
     >            fbpom * (ybn * outcbn(i) +
     >            ybp * outcbp(i))) - outcpb(i) - outcps(i)
c
c ----------------------------------------------------------------------
c change in stablized organic matter (passive pool) from growth
c in microbial biomass, and changes in protected and unprotected
c SOM
c
c add a loss of C due to leaching out of the profile, based
c on approximation of CENTURY model below 1 m in depth
c based on water in the profile, and texture of soil
c tuned to known outputs or leaching that has been measured in the field
c at Arlington-WI (Courtesy K. Brye, MS) and applied to the global scale
c on average, this calibration yields about 10-50 Kg C ha-1 yr-1 leaching
c depending on C in soil...will need to be tied to an amount of water
c flowing through the profile based upon precipitation eventually
c ----------------------------------------------------------------------
c
         h20    = 0.30e-03
c
c h20 is a constant relating to the flow of water through the top 1 m of the
c profile 
c use texfact -- the % sand -- or texture factor effect on leaching (see Parton
c et al. (1991) calculated from the average sand content of top 1 m of soil
c in the model
c
        fleach = h20/18.0 * (0.01 + 0.04 * texfact)
c
c --------------------------------------------------------------------
c change in passive organic carbon pool
c ---------------------------------------------------------------------
c
       dcsdt(i) = ((yns * outcns(i)) + (yps * outcps(i))) -
     >            outcsb(i) -  (fleach * csoipas(i))
c
       cleach(i) = fleach * csoipas(i) + fleach * csoislop(i) +
     >             fleach * csoislon(i)
c
       ynleach(i) = ynleach(i) + fleach * csoipas(i)/cnr(2) +
     >              fleach * csoislop(i)/cnr(3) +
     >              fleach * csoislon(i)/cnr(4)
c
c update slow pools of carbon for leaching losses
c
       dcndt(i) = dcndt(i) - fleach * csoislon(i)		
       dcpdt(i) = dcpdt(i) - fleach * csoislop(i)		
c
      if (spin .eq. spinmax) then

         yrleach(i) =  cleach(i) + yrleach(i)
c
      endif
c
c ---------------------------------------------------------------------
c calculate the amount of net N mineralization or immobilization
c ---------------------------------------------------------------------
c
c uptake of n by growth of microbial biomass
c
c immobilized n used for requirements of microbial growth
c is based on flow of carbon and the difference of C/N ratio of
c the microbes and their efficiency versus the C/N ratio of the
c material that is being decomposed 
c
c
c ------------------------------
c structural root decomposition 
c ------------------------------
c
       if (yrs/cnr(1) .gt. 1./cnr(6)) then
         nbiors(i) = (1./cnr(6) - yrs/cnr(1))
     >                * outcrs(i)
         nminrs(i) = 0.0
c
       else
         nminrs(i) = (1./cnr(6) - yrs/cnr(1))
     >                * outcrs(i)
         nbiors(i) = 0.0
       endif
c
c ------------------------------
c structural leaf decomposition
c ------------------------------
c
       if (yls/cnr(1) .gt. 1./cnr(6)) then
         nbiols(i) = (1./cnr(6) - yls/cnr(1))
     >                * outcls(i)
         nminls(i) = 0.0
c
       else
         nminls(i) = (1./cnr(6) - yls/cnr(1))
     >                * outcls(i)
         nbiols(i) = 0.0
       endif
c
c
c ------------------------------
c structural wood decomposition
c ------------------------------
c	
       if (yws/cnr(1) .gt. 1./cnr(8)) then
         nbiows(i) = (1./cnr(8) - yws/cnr(1))
     >                * outcws(i)
         nminws(i) = 0.0
c
       else
         nminws(i) = (1./cnr(8) - yws/cnr(1))
     >                * outcws(i)
         nbiows(i) = 0.0
       endif
c
c ------------------------------
c metabolic wood decomposition
c ------------------------------
c
       if (ywm/cnr(1) .gt. 1./cnr(8)) then
         nbiowm(i) = (1./cnr(8) - ywm/cnr(1))
     >                * outcwm(i)
         nminwm(i) = 0.0
c
       else
         nminwm(i) = (1./cnr(8) - ywm/cnr(1))
     >                * outcwm(i)
         nbiowm(i) = 0.0
       endif
c
c ------------------------------
c metabolic leaf decomposition
c ------------------------------
c
       if (ylm/cnr(1) .gt. 1./cnr(7)) then
         nbiolm(i) = (1./cnr(7) - ylm/cnr(1))
     >                * outclm(i)
         nminlm(i) = 0.0
c
       else
         nminlm(i) = (1./cnr(7) - ylm/cnr(1))
     >                * outclm(i)
         nbiolm(i) = 0.0
       endif
c
c
c ------------------------------
c metabolic root decomposition
c ------------------------------
c
       if (yrm/cnr(1) .gt. 1./cnr(7)) then
         nbiorm(i) = (1./cnr(7) - yrm/cnr(1))
     >                * outcrm(i)
         nminrm(i) = 0.0
c
       else
         nminrm(i) = (1./cnr(7) - yrm/cnr(1))
     >                * outcrm(i)
         nbiorm(i) = 0.0
       endif
c
c
c ----------------------------------------------
c non-protected organic matter decomposition
c ----------------------------------------------
c
       if (ynb/cnr(1) .gt. 1./cnr(4)) then
         nbioslon(i) = (1./cnr(4) - ynb/cnr(1))
     >                  * outcnb(i)
         nminslon(i) = 0.0
c
       else
         nminslon(i) = (1./cnr(4) - ynb/cnr(1))
     >                  * outcnb(i)
         nbioslon(i) = 0.0
       endif
c
c
c ----------------------------------------------
c protected organic matter decomposition
c ----------------------------------------------
c
       if (ypb/cnr(1) .gt. 1./cnr(3)) then
         nbioslop(i) = (1./cnr(3) - ypb/cnr(1))
     >                  * outcpb(i)
         nminslop(i) = 0.0
c
       else
         nminslop(i) = (1./cnr(3) - ypb/cnr(1))
     >                  * outcpb(i)
         nbioslop(i) = 0.0
       endif
c
c
c ----------------------------------------------
c stablized organic matter decomposition
c ----------------------------------------------
c
       if (ysb/cnr(1) .gt. 1./cnr(2)) then
         nbiopas(i) = (1./cnr(2) - ysb/cnr(1))
     >                 * outcsb(i)
         nminpas(i) = 0.0
c
       else
         nminpas(i) = (1./cnr(2) - ysb/cnr(1))
     >                 * outcsb(i)
         nbiopas(i) = 0.0
       endif
c
c ----------------------------------------------
c total immobilized N used for biomass growth
c ----------------------------------------------
c
       totimm(i) = nbiors(i) + nbiols(i) + nbiows(i) + nbiowm(i)
     >           + nbiolm(i) + nbiorm(i) + nbioslon(i) + nbioslop(i)
     >           + nbiopas(i)
c
c -----------------------------------------------------------------------------
c gross amount of N mineralized by decomposition of C by microbial biomass
c assume that N is attached to the flow of C by the C/N ratio of the substrate
c also assume that the amount of N attached to CO2 that is respired is also
c mineralized (i.e. the amount of N mineralized is related to the total outflow
c of carbon, and not the efficiency or yield)..see Parton et al., 1987
c -----------------------------------------------------------------------------
c
       totmin(i) = nminrs(i) + nminls(i) + nminws(i) + nminwm(i)
     >           + nminlm(i) + nminrm(i) + nminslon(i) + nminslop(i)
     >           + nminpas(i)
c
c -----------------------------------------------------------------------------
c when carbon is transferred from one pool to another, each pool has a distinct
c C:N ratio.  In the case of pools where carbon is moving from the pool to 
c the microbial biomass (used for growth/assimilation), net mineralization
c takes place (N is released) after the requirements of building the biomass
c are met.  In the cases of other transformations of C, N is not conserved
c if it follows from one pool to another which has a different C:N ratio;
c either N is released or is needed to make the transformation and keep N
c conserved in the model. 
c
c other calculations of either N release or immobilization to keep track of
c the budget
c
        nrelps(i) = outcps(i) * (1./cnr(3) - 1./cnr(2))
        nrelns(i) = outcns(i) * (1./cnr(4) - 1./cnr(2))
        nrelbn(i) = (1.-fbpom) * outcbn(i) * (1./cnr(1) - 1./cnr(4)) +
     >              (1.-fbpom) * outcbp(i) * (1./cnr(1) - 1./cnr(4))
        nrelbp(i) = fbpom * outcbp(i) * (1./cnr(1) - 1./cnr(3)) +
     >              fbpom * outcbn(i) * (1./cnr(1) - 1./cnr(3))
        nrelll(i) = lig_frac * outcll(i) * (1./cnr(5) - 1./cnr(3)) +
     >              lig_frac * outcll(i) * (1./cnr(5) - 1./cnr(4))
        nrelrl(i) = lig_frac * outcrl(i) * (1./cnr(5) - 1./cnr(3)) +
     >              lig_frac * outcrl(i) * (1./cnr(5) - 1./cnr(4))
        nrelwl(i) = lig_frac * outcwl(i) * (1./cnr(5) - 1./cnr(3)) +
     >              lig_frac * outcwl(i) * (1./cnr(5) - 1./cnr(4))
c
        totnrel(i) = nrelps(i) + nrelns(i) + nrelbn(i) +
     >               nrelbp(i) + nrelll(i) + nrelrl(i) + nrelwl(i)
c
c -----------------------------------------------------------------------------
c calculate whether net mineralization or immobilization occurs
c on a grid cell basis -- tnmin is an instantaneous value for each time step
c it is passed along to stats to calculate, daily, monthly and annual totals
c of nitrogen mineralization
c
c this is for mineralization/immobilization that is directly related to 
c microbial processes (oxidation of carbon)
c
c the value of totnrel(i) would need to be added to complete the budget
c of N in the model. Because it can add/subtract a certain amount of N
c from the amount of net mineralization.  However, these transformations
c are not directly related to microbial decomposition, so do we add them
c into the value or not?
c -----------------------------------------------------------------------------
c
           netmin(i) = totmin(i) + totimm(i) + totnrel(i) 
           if (netmin(i) .gt. 0.0) then
c
              tnmin(i) = netmin(i)
c	
           else
c	  
              tnmin(i) = 0.0 
c
           endif
c
c convert value of tnmin of Kg-N/m2/dtime to mole-N/s
c based on N = .014 Kg/mole -- divide by the number of seconds in daily timestep
c
            tnmin(i) = tnmin(i)/(86400. * 0.014)
c
c ---------------------------------------------------
c update soil c pools for transformations of c and n
c ---------------------------------------------------
c
           totcmic(i)  = max(totcmic(i)  + dbdt(i), 0.0)
           csoislon(i) = max(csoislon(i) + dcndt(i),0.0)
           csoislop(i) = max(csoislop(i) + dcpdt(i),0.0)
           csoipas(i)  = max(csoipas(i)  + dcsdt(i),0.0)
           clitlm(i)   = max(clitlm(i)  - outclm(i),0.0)
           clitls(i)   = max(clitls(i)  - outcls(i),0.0)
           clitll(i)   = max(clitll(i)  - outcll(i),0.0)
           clitrm(i)   = max(clitrm(i)  - outcrm(i),0.0)
           clitrs(i)   = max(clitrs(i)  - outcrs(i),0.0)
           clitrl(i)   = max(clitrl(i)  - outcrl(i),0.0)
           clitwm(i)   = max(clitwm(i)  - outcwm(i),0.0)
           clitws(i)   = max(clitws(i)  - outcws(i),0.0)
           clitwl(i)   = max(clitwl(i)  - outcwl(i),0.0)
c
c -----------------------------------------------------------
c update soil n pools based on c:n ratios of each pool
c this approach is assuming that the c:n ratios are remaining
c constant through the simulation. flow of nitrogen is attached
c to carbon 
c -----------------------------------------------------------
c
           totnmic(i)  = totcmic(i) /cnr(1)
           nsoislon(i) = csoislon(i)/cnr(4)
           nsoislop(i) = csoislop(i)/cnr(3)
           nsoipas(i)  = csoipas(i) /cnr(2)
           nlitlm(i)   = clitlm(i)  /cnr(7)
           nlitls(i)   = clitls(i)  /cnr(6)
           nlitll(i)   = clitll(i)  /cnr(5)
           nlitrm(i)   = clitrm(i)  /cnr(7)
           nlitrs(i)   = clitrs(i)  /cnr(6)
           nlitrl(i)   = clitrl(i)  /cnr(5)
           nlitwm(i)   = clitwm(i)  /cnr(8)
           nlitws(i)   = clitws(i)  /cnr(8)
           nlitwl(i)   = clitwl(i)  /cnr(8)
c
c total above and belowground litter
c
           totlit(i) =  clitlm(i) + clitls(i) + clitll(i) +
     >                  clitrm(i) + clitrs(i) + clitrl(i) +
     >                  clitwm(i) + clitws(i) + clitwl(i)
c
c sum total aboveground litter (leaves and wood)
c
           totalit(i) = clitlm(i) + clitls(i) + clitwm(i) +
     >                  clitll(i) + clitws(i) + clitwl(i)
c
c sum total belowground litter (roots) 
c
           totrlit(i) = clitrm(i) + clitrs(i) + clitrl(i)
c
c determine total soil carbon amounts (densities are to 1 m depth; Kg/m-2)
c	
           totcsoi(i) = csoipas(i) + csoislop(i) +
     >                  totcmic(i) + csoislon(i)
c
c calculate total amount of litterfall occurring (total for year)
c
           totfall(i) = falll(i) + fallr(i) + fallw(i)
c
c nitrogen 
c
c total nitrogen in litter pools (above and belowground)
c
          totnlit(i) =  nlitlm(i) + nlitls(i) + nlitrm(i) + nlitrs(i) +
     >                  nlitwm(i) + nlitws(i) + nlitll(i) + nlitrl(i) +
     >                  nlitwl(i)
c
c sum total aboveground litter   (leaves and wood)
c
          totanlit(i) = nlitlm(i) + nlitls(i) + nlitwm(i) +
     >                  nlitll(i) + nlitws(i) + nlitwl(i)
c
c sum total belowground litter  (roots)
c
          totrnlit(i) = nlitrm(i) + nlitrs(i) + nlitrl(i)
c
c total soil nitrogen to 1 m depth (kg-N/m**2)
c
          totnsoi(i) = nsoislop(i) + nsoislon(i) +
     >                 nsoipas(i)  + totnmic(i) + totnlit(i)
c
c --------------------------------------------------------------------------
c calculate running sum of yearly net mineralization, and nitrogen in pool
c available to plants for uptake--during spin up period, can only count one
c of the cycles for each timestep--otherwise false additions will result
c values of yearly mineralization are in Kg/m-2
c --------------------------------------------------------------------------
c
        if (spin .eq. spinmax) then
c
          storedn(i)  = storedn(i) + tnmin(i)
c
        endif
c
c calculate total amount of carbon in soil at end of cycle
c this is used to help calculate the amount of carbon that is respired
c by decomposing microbial biomass
c
        totcend(i) = totlit(i) + totcsoi(i)
c
c --------------------------------------------------------------------------
c the amount of co2resp(i) is yearly value and is dependent on the amount
c of c input each year, the amount in each pool at beginning of the year,
c and the amount left in the pool at the end of the year
c along with the amount of root respiration contributing to the flux from
c calculations performed in stats.f
c --------------------------------------------------------------------------
c
        if (spin .eq. spinmax) then 
c
c --------------------------------------------------------------------------
c only count the last cycle in the spin-up for co2soi
c when the iyear is less than the nspinsoil value...otherwise
c an amount of CO2 respired will be about 10 times the actual
c value because this routine is called articially 10 extra times
c each time step to spin up the soil carbon
c
c add n-deposition due to rainfall once each day, and
c the amount of N fixed through N-fixers.  These equations
c are based on the annual precip input (cm) and are from
c the CENTURY model...Parton et al., 1987.
c The base equations are in units of (g) N m-2 so have to
c divide by 1000 to put in units of Kg.
c
c the values in the equation of 0.21 and -0.18 were adjusted to reflect
c average daily inputs when no precipitation was falling - the original
c constants are for the entire year 
c --------------------------------------------------------------------------
c
            deposn(i)    = (0.0005753  + 0.0028 * (precip(i)*0.1))*1.e-3
            fixsoin(i)   = (-0.0004932 + 0.14   * (precip(i)*0.1))*1.e-3
c
c --------------------------------------------------------------------------
c add to the daily total of co2 flux leaving the soil from microbial
c respiration -- instantaneous value for each timestep
c since this subroutine gets called daily...instantaneous fluxes
c the fluxes need to be put on a per second basis, which will be dependent
c on the timestep.  Furthermore, because the biogeochem subroutine does
c not get called each timestep...an approximation for a timestep average
c microbial flux and nmineralization rate will be applied
c --------------------------------------------------------------------------
c
c calculate daily co2 flux due to microbial decomposition
c
          tco2mic(i) = totcbegin(i) + totcin(i) - totcend(i) - cleach(i)
c
c convert co2 flux from kg C/day  (seconds in a daily timestep) to mol-C/s
c based on .012 Kg C/mol
c
          tco2mic(i) = tco2mic(i)/(86400. * 0.012)
c
        endif
c
 100  continue
c
c return to main
c
        return
        end
c
