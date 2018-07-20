c
c #    #  ######   ####   ######   #####    ##     #####     #     ####   #    #
c #    #  #       #    #  #          #     #  #      #       #    #    #  ##   #
c #    #  #####   #       #####      #    #    #     #       #    #    #  # #  #
c #    #  #       #  ###  #          #    ######     #       #    #    #  #  # #
c  #  #   #       #    #  #          #    #    #     #       #    #    #  #   ##
c   ##    ######   ####   ######     #    #    #     #       #     ####   #    #
c
c ---------------------------------------------------------------------
      subroutine pheno
c ---------------------------------------------------------------------
c
c common blocks
c
      use implicit
c
      use compar
      use comatm
      use comsoi
      use comsum
      use comveg
c
c local variables
c
      integer
     >  i,          !
     >  imonth,     !
     >  iday        !
c
      real
     >  ddays,        !
     >  ddfac,        !
     >  tthreshold,   ! temperature threshold for budburst and senescence
     >  gthreshold,   ! temperature threshold for budburst and senescence
     >  avglaiu,      ! average lai of upper canopy 
     >  avglail,      ! average lai of lower canopy 
**** DTP 2000/06/28 Modified this following discussion with Navin. We now
*    retain fu(i) derived from dynaveg and constrain it to a local value
*    "fu_phys" in the range 0.25 to 0.975 used in the canopy physics calcs.
     >  fu_phys       ! Local value of fu(i) constrained to range 0.25 to 0.975
                      ! to keep physics calculations stable.
c
c define 'drop days' -- number of days to affect phenology change
c
      ddays = 15.0
      ddfac = 1.0 / ddays
c
c begin global grid
c
      do 100 i = 1, npoi
c
c ---------------------------------------------------------------------
c * * * upper canopy winter phenology * * *
c ---------------------------------------------------------------------
c
c temperature threshold for budburst and senescence
c
c temperature threshold is assumed to be 0 degrees C 
c or 5 degrees warmer than the coldest monthly temperature
c
        tthreshold = max (0.0         + 273.16,
     >                    tc(i) + 5.0 + 273.16)
c
c gdd threshold temperature for leaf budburst
c with a growing degree threshold of 100 units
c
        gthreshold = 0.0 + 273.16
c 
c determine if growing degree days are initiated
c
        if (a10td(i).lt.gthreshold) then
          agddu(i)  = 0.0
        else
          agddu(i) = agddu(i) + td(i) - gthreshold
        endif
c
c determine leaf display
c
        if (a10td(i).lt.tthreshold) then
          tempu(i)  = max (0.0, tempu(i) - ddfac)
        else
          tempu(i) = min (1., max (0.0, agddu(i) - 100.0) / 50.0)
        endif
c
c ---------------------------------------------------------------------
c * * * lower canopy winter phenology * * *
c ---------------------------------------------------------------------
c
c temperature threshold for budburst and senescence
c
c temperature threshold is assumed to be 0 degrees C 
c
        tthreshold = 0.0 + 273.16
c
c gdd threshold temperature for leaf budburst
c with a growing degree threshold of 150 units
c
        gthreshold = -5.0 + 273.16
c 
c determine if growing degree days are initiated
c
        if (a10td(i).lt.gthreshold) then
          agddl(i)  = 0.0
        else
          agddl(i) = agddl(i) + td(i) - gthreshold
        endif
c
c determine leaf display
c
        if (a10td(i).lt.tthreshold) then
          templ(i)  = max (0.0, templ(i) - ddfac)
        else
          templ(i) = min (1., max (0.0, agddl(i) - 150.0) / 50.0)
        endif
c
c ---------------------------------------------------------------------
c * * * drought canopy winter phenology * * *
c ---------------------------------------------------------------------
c
        if (a10ancub(i).lt.0.0) dropu(i) = max (0.1, dropu(i) - ddfac)
        if (a10ancub(i).ge.0.0) dropu(i) = min (1.0, dropu(i) + ddfac)
c
        if (a10ancls(i).lt.0.0) dropls(i) = max (0.1, dropls(i) - ddfac)
        if (a10ancls(i).ge.0.0) dropls(i) = min (1.0, dropls(i) + ddfac)
c
        if (a10ancl4(i).lt.0.0) dropl4(i) = max (0.1, dropl4(i) - ddfac)
        if (a10ancl4(i).ge.0.0) dropl4(i) = min (1.0, dropl4(i) + ddfac)
c
        if (a10ancl3(i).lt.0.0) dropl3(i) = max (0.1, dropl3(i) - ddfac)
        if (a10ancl3(i).ge.0.0) dropl3(i) = min (1.0, dropl3(i) + ddfac)
c
c ---------------------------------------------------------------------
c * * * update lai and canopy fractions * * *
c ---------------------------------------------------------------------
c
c upper canopy single sided leaf area index (area-weighted)
c
        avglaiu = plai(i,1)             +
     >            plai(i,2) * dropu(i)  +
     >            plai(i,3)             +
     >            plai(i,4)             +
     >            plai(i,5) * tempu(i)  +
     >            plai(i,6)             +
     >            plai(i,7) * tempu(i)  +
     >            plai(i,8) * tempu(i)
c
c upper canopy fractions
c
        frac(i,1) = plai(i,1)            / max (avglaiu, epsilon)
        frac(i,2) = plai(i,2) * dropu(i) / max (avglaiu, epsilon)
        frac(i,3) = plai(i,3)            / max (avglaiu, epsilon)
        frac(i,4) = plai(i,4)            / max (avglaiu, epsilon)
        frac(i,5) = plai(i,5) * tempu(i) / max (avglaiu, epsilon)
        frac(i,6) = plai(i,6)            / max (avglaiu, epsilon)
        frac(i,7) = plai(i,7) * tempu(i) / max (avglaiu, epsilon)
        frac(i,8) = plai(i,8) * tempu(i) / max (avglaiu, epsilon)
c
c lower canopy single sided leaf area index (area-weighted)
c
        avglail = plai(i,9)                              +
     >            plai(i,10) * min (templ(i), dropls(i)) +
     >            plai(i,11) * min (templ(i), dropl4(i)) +
     >            plai(i,12) * min (templ(i), dropl3(i))
c
c lower canopy fractions
c
        frac(i,9)  = plai(i,9)                              /
     >               max (avglail, epsilon)
c
        frac(i,10) = plai(i,10) * min (templ(i), dropls(i)) /
     >               max (avglail, epsilon)
c
        frac(i,11) = plai(i,11) * min (templ(i), dropl4(i)) /
     >               max (avglail, epsilon)
c
        frac(i,12) = plai(i,12) * min (templ(i), dropl3(i)) /
     >               max (avglail, epsilon)
c
c calculate the canopy leaf area index using the fractional vegetation cover
c
        lai(i,1) = avglail / fl(i)

**** DTP 2000/06/28 Modified this following discussion with Navin. We now
*    retain fu(i) derived from dynaveg and constrain it to a local value
*    "fu_phys" in the range 0.25 to 0.975 used in the canopy physics calcs.

        fu_phys = max (0.25, min (0.975, fu(i)))
        lai(i,2) = avglaiu / fu_phys
        lai(i,2) = avglaiu / fu(i)
c
c put a fix on canopy lais to avoid problems in physics
c
        lai(i,1) = min (lai(i,1), 12.0)
        lai(i,2) = min (lai(i,2), 12.0)
c
c ---------------------------------------------------------------------
c * * * update canopy height parameters * * *
c ---------------------------------------------------------------------
c
c update lower canopy height parameters
c
c note that they are based on vegetation fraction and not
c averaged over the entire gridcell
c
        zbot(i,1)   =  0.05
        ztop(i,1)   =  max (0.25, lai(i,1) * 0.25)
c        
c constrain ztop to be at least 0.5 meter lower than 
c zbot for upper canopy
c
        ztop(i,1) = min (ztop(i,1), zbot(i,2) - 0.5)
c
c end of loop
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
      subroutine dynaveg (isimfire) ! , isim_ac, year)
c ---------------------------------------------------------------------
c
      use implicit
c
      use compar
      use comsoi
      use comsum
      use comveg
      use compft
      use comage
c
c Arguments
c
      integer isimfire  ! fire switch
!     >        isim_ac,   ! age-class dynamics switch
!     >        year      ! year of simulation

      real    pfire,     ! probability of fire -- should be determined externally.
     >        pdist      ! probability of other disturbance types....
       
      PARAMETER (pfire = 1.0) ! for now we just assume it occurs all the time
      
c
c local variables
c
      integer
     >  i, j           ! gridcell counter
c
      real
     >  sapspeed,      ! in mm/day
     >  trans,         ! (2.5 mm/day) 
     >  saparea,       ! in m**2
     >  sapvolume,     ! in m**3
     >  denswood,      ! kg/m**3
     >  wood,          ! total amount of woody biomass in gridcell
     >  taufin         !
*     >  xminlai        !
c
*      real
*     >  aleaf(npft),   ! allocation fraction to leaves
*     >  aroot(npft),   ! allocation fraction to fine roots
*     >  awood(npft),   ! allocation fraction to wood
*     >  tauleaf(npft), ! turnover time of carbon in leaves (years)
*     >  tauroot(npft), ! turnover time of carbon in fine roots (years)
*     >  tauwood(npft)   ! turnover time of carbon in wood (years)
*     >  tauwood0(npft) ! normal (unstressed) turnover time
c
c ibis uses a small number of plant functional types:
c
c  1: tropical broadleaf evergreen tree
c  2: tropical broadleaf drought-deciduous trees
c  3: warm-temperate broadleaf evergreen tree
c  4: temperate conifer evergreen tree
c  5: temperate broadleaf cold-deciduous tree
c  6: boreal conifer evergreen tree
c  7: boreal broadleaf cold-deciduous tree
c  8: boreal conifer cold-deciduous tree
c  9: evergreen shrub
c 10: deciduous shrub
c 11: warm (c4) grass
c 12: cool (c3) grass
c
c ---------------------------------------------------------------------
c * * * specify biomass turnover parameters (years) * * *
c ---------------------------------------------------------------------
c
*      data tauleaf / 1.01,   ! tropical broadleaf evergreen trees
*     >               1.00,   ! tropical broadleaf drought-deciduous trees
*     >               1.00,   ! warm-temperate broadleaf evergreen trees
*     >               2.00,   ! temperate conifer evergreen trees
*     >               1.00,   ! temperate broadleaf cold-deciduous trees
*     >               2.50,   ! boreal conifer evergreen trees
*     >               1.00,   ! boreal broadleaf cold-deciduous trees
*     >               1.00,   ! boreal conifer cold-deciduous trees
*     >               1.50,   ! evergreen shrubs
*     >               1.00,   ! deciduous shrubs
*     >               1.25,   ! warm (c4) grasses
*     >               1.50 /  ! cool (c3) grasses
c
*      data tauwood0 / 25.0,  ! tropical broadleaf evergreen trees
*     >                25.0,  ! tropical broadleaf drought-deciduous trees
*     >                25.0,  ! warm-temperate broadleaf evergreen trees
*     >                50.0,  ! temperate conifer evergreen trees
*     >                50.0,  ! temperate broadleaf cold-deciduous trees
*     >               100.0,  ! boreal conifer evergreen trees
*     >               100.0,  ! boreal broadleaf cold-deciduous trees
*     >               100.0,  ! boreal conifer cold-deciduous trees
*     >                 5.0,  ! evergreen shrubs
*     >                 5.0,  ! deciduous shrubs
*     >               999.0,  ! warm (c4) grasses
*     >               999.0 / ! cool (c3) grasses
c
c begin global grid
c
 999  do 100 i = 1, npoi
c
c ---------------------------------------------------------------------
c * * * initialize vegetation dynamics pools * * *
c ---------------------------------------------------------------------
c
c zero out litter fall fields
c
        falll(i) = 0.0
        fallr(i) = 0.0
        fallw(i) = 0.0
c
c zero out carbon lost due to disturbance
c 
        cdisturb(i) = 0.0
c
        wood = 0.001
c
c ---------------------------------------------------------------------
c * * * update npp, and pool losses  * * *
c ---------------------------------------------------------------------
c
c go through all the pfts
c
        do 110 j = 1, npft
c
c apply this year's existence arrays to npp
c
          aynpp(i,j)  = exist(i,j) * aynpp(i,j)
c
c determine above-ground npp for each plant type
c
          ayanpp(i,j) = (aleaf(j) + awood(j)) * aynpp(i,j)
c
c determine turnover rates for woody biomass:
c
c if pft can exist,    then tauwood = tauwood0 (normal turnover),
c if pft cannot exist, then tauwood = taufin years (to kill off trees)
c
c          taufin     = 5.0
           taufin     = tauwood0(j)/2.0
c
          tauwood(j) = tauwood0(j) - (tauwood0(j) - taufin) *
     >                               (1.0 - exist(i,j))
c
c assume a constant fine root turnover time
c
*          tauroot(j) = 1.0
c
c determine litter fall rates
c
          falll(i) = falll(i) + cbiol(i,j) / tauleaf(j)
          fallr(i) = fallr(i) + cbior(i,j) / tauroot(j)
          fallw(i) = fallw(i) + cbiow(i,j) / tauwood(j)
c
c ---------------------------------------------------------------------
c * * * update biomass pools  * * *
c ---------------------------------------------------------------------
c
c update carbon reservoirs using an analytical solution
c to the original carbon balance differential equation
c
          cbiol(i,j) = cbiol(i,j) * exp(-1./tauleaf(j))  +
     >                 aleaf(j) * tauleaf(j) * max (0., aynpp(i,j)) *
     >                 (1. - exp(-1./tauleaf(j)))
c
          cbiow(i,j) = cbiow(i,j) * exp(-1./tauwood(j))  +
     >                 awood(j) * tauwood(j) * max (0., aynpp(i,j)) *
     >                 (1. - exp(-1./tauwood(j)))
c
          cbior(i,j) = cbior(i,j) * exp(-1./tauroot(j))  +
     >                 aroot(j) * tauroot(j) * max (0., aynpp(i,j)) *
     >                 (1. - exp(-1./tauroot(j)))
c
          if (j.le.8) wood = wood + max (0.0, cbiow(i,j))
c
 110    continue
c
c ---------------------------------------------------------------------
c * * * apply disturbances * * *
c ---------------------------------------------------------------------
c
c set fixed disturbance regime
c
        disturbf(i) = 0.005
        disturbo(i) = 0.005
        
**** DTP 2000/08/10. One can do a decent test of ACME by setting
*    these disturbance rates to zero. With these values, the area
*    disturbed each year will be zero so the distribution of biomass
*    and PFTs across the domain should be identical to those 
*    resulting from a run of standard IBIS with zero disturbance

*        disturbf(i) = 0.0  ! Test with zero disturbance rate 
*        disturbo(i) = 0.0  ! (This should equal reference sim).

c
c call fire disturbance routine
**** DTP 2001/03/06: In general isimfire should be set to zero if isim_ac
*    is set to 1 (but what does isimfire really do?)
c
        if (isimfire.eq.1) call fire

 
          do 116 j = 1, npft 
c
c calculate biomass (vegetations) carbon lost to atmosphere   
c used to balance net ecosystem exchange  
c
c ---------------------------------------------------------------------
**** DTP 2000/04/22 QUESTION: 
c ---------------------------------------------------------------------
* Shouldn't a portion of the destroyed material be added to litter fall?

            cdisturb(i) = cdisturb(i) + 
     >                    cbiol(i,j) * (disturbf(i) + disturbo(i)) +
     >                    cbiow(i,j) * (disturbf(i) + disturbo(i)) +
     >                    cbior(i,j) * (disturbf(i) + disturbo(i))                  
c          
c adjust biomass pools due to disturbances
c
            cbiol(i,j) = cbiol(i,j) * (1. - disturbf(i) - disturbo(i))
            cbiow(i,j) = cbiow(i,j) * (1. - disturbf(i) - disturbo(i))
            cbior(i,j) = cbior(i,j) * (1. - disturbf(i) - disturbo(i))
c
c constrain biomass fields to be positive
c
            cbiol(i,j) = max (0.0, cbiol(i,j))
            cbiow(i,j) = max (0.0, cbiow(i,j))
            cbior(i,j) = max (0.0, cbior(i,j))

 116      continue



c ---------------------------------------------------------------------
c * * * check and update biomass pools following disturbance * * *
c ---------------------------------------------------------------------
c
        do 120 j = 1, npft
c
c maintain minimum value of leaf carbon in areas that plants exist
c
*          xminlai = 0.010
c
          cbiol(i,j) = max (exist(i,j) * xminlai / specla(j),
     >                      cbiol(i,j))
c
c update vegetation's physical characteristics
c
          plai(i,j)    = cbiol(i,j) * specla(j)
          biomass(i,j) = cbiol(i,j) + cbiow(i,j) + cbior(i,j)
c
 120    continue
c
c ---------------------------------------------------------------------
c * * * update annual npp, lai, and biomass * * *
c ---------------------------------------------------------------------
c
c adjust annual net ecosystem exchange (calculated in stats.f) 
c by loss of carbon to atmosphere due to biomass burning (fire)
c
        ayneetot(i) = ayneetot(i) - cdisturb(i)
c
c determine total ecosystem above-ground npp
c
        ayanpptot(i) = ayanpp(i,1)  + ayanpp(i,2) +
     >                 ayanpp(i,3)  + ayanpp(i,4) +
     >                 ayanpp(i,5)  + ayanpp(i,6) +
     >                 ayanpp(i,7)  + ayanpp(i,8) +
     >                 ayanpp(i,9)  + ayanpp(i,10) +
     >                 ayanpp(i,11) + ayanpp(i,12)
c
c update total canopy leaf area
c
        totlaiu(i) = plai(i,1)  + plai(i,2) +
     >               plai(i,3)  + plai(i,4) +
     >               plai(i,5)  + plai(i,6) +
     >               plai(i,7)  + plai(i,8)
c
        totlail(i) = plai(i,9)  + plai(i,10) +
     >               plai(i,11) + plai(i,12)
c
c update total biomass
c
        totbiou(i) = biomass(i,1) +
     >               biomass(i,2) +
     >               biomass(i,3) +
     >               biomass(i,4) +
     >               biomass(i,5) +
     >               biomass(i,6) +
     >               biomass(i,7) +
     >               biomass(i,8)
c
        totbiol(i) = biomass(i,9)  +
     >               biomass(i,10) +
     >               biomass(i,11) +
     >               biomass(i,12)
c
c ---------------------------------------------------------------------
c * * * update fractional cover and vegetation height parameters * * *
c ---------------------------------------------------------------------
c
c
**** Added these in temporarily for comparison with original code.
**** Delete these from production version....
        fu(i) = (1.0 - exp(-wood)) / (1.0 - exp(-woodnorm))
        fu(i) = fu(i) * (1. - disturbf(i) - disturbo(i))

c constrain the fractional cover (upper canopy)
c
        fu(i) = max (0.25, min (0.975, fu(i)))
c
c update fractional cover of herbaceous (lower) canopy:
c 
        fl(i) = totlail(i) / 1.0
c
c apply disturbances to fractional cover (lower canopy)
c
        fl(i) = fl(i) * (1. - disturbf(i) - disturbo(i))
c
c constrain the fractional cover (lower canopy)
c
        fl(i) = max (0.25, min (0.975, fl(i)))
c
c
c annual update upper canopy height parameters
c should be calculated based on vegetative fraction and not the
c average over the entire grid cell
c
        zbot(i,2) = 3.0
        ztop(i,2) = max(zbot(i,2) + 1.00, 2.50 *
     >                  totbiou(i) / fu(i) * 0.75)
c
c ---------------------------------------------------------------------
c * * * update stem area index and sapwood fraction * * *
c ---------------------------------------------------------------------
c
c estimate stem area index (sai) as a fraction of the lai
c
        sai(i,1) = 0.050 * totlail(i)
        sai(i,2) = 0.250 * totlaiu(i)
c
c estimate sapwood fraction of woody biomass
c
        sapspeed  = 25.0                        ! (m/day)
        trans     = 0.0025                      ! (2.5 mm/day) 
        saparea   = (trans / sapspeed)          ! m**2
c
        sapvolume = saparea * ztop(i,2) * 0.75  ! m**3
c
        denswood  = 400.0                       ! kg/m**3
c
        sapfrac(i) = min (0.50, max (0.05, sapvolume * denswood / wood))
c
 100  continue
c
c ---------------------------------------------------------------------
c * * * map out vegetation classes for this year * * *
c ---------------------------------------------------------------------
c
      call vegmap
c
c
c return to the main program
c
      return

      end  ! DYNAVEG
c
c
c ---------------------------------------------------------------------
      subroutine fire
c ---------------------------------------------------------------------
c
      use implicit
c
      use compar
      use comveg
c
c local variables
c
      integer i
c
      real burn
c
c begin global grid
c
      do 100 i = 1, npoi
c
        burn = firefac(i) * min (1.0, totlit(i) / 0.200)
c
        disturbf(i) = 1.0 - exp(-0.5 * burn)
c
        disturbf(i) = max (0.0, min (1.0, disturbf(i)))
c
 100  continue
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine vegmap
c ---------------------------------------------------------------------
c
      use implicit
c
      use compar
      use comveg
c
c local variables
c
      integer 
     >     i, j,            ! loop indice
     >     domtree          ! dominant tree
c
      real maxlai,          ! maximum lai
     >     totlai,          ! total ecosystem lai
     >     grassfrac,       ! fraction of total lai in grasses
     >     treefrac,        ! fraction of total lai in trees
     >     treelai,         ! lai of trees
     >     shrublai,        ! lai of shrubs
     >     grasslai,        ! lai of grass
     >     ratio
c
c classify vegetation cover into standard ibis vegetation classes 
c
c ---------------------------------------------------
c  1: tropical evergreen forest / woodland
c  2: tropical deciduous forest / woodland
c  3: temperate evergreen broadleaf forest / woodland
c  4: temperate evergreen conifer forest / woodland
c  5: temperate deciduous forest / woodland
c  6: boreal evergreen forest / woodland
c  7: boreal deciduous forest / woodland
c  8: mixed forest / woodland
c  9: savanna
c 10: grassland / steppe 
c 11: dense shrubland
c 12: open shrubland
c 13: tundra
c 14: desert 
c 15: polar desert / rock / ice
c ---------------------------------------------------
c
c begin global grid
c
      do 100 i = 1, npoi
c
c determine total lai and tree, shrub, and grass fractions
c
        treelai   = totlaiu(i) 
        shrublai  = plai(i,9)  + plai(i,10)
        grasslai  = plai(i,11) + plai(i,12)
c
        totlai    = max (0.01, totlail(i) + totlaiu(i))
c
c determine dominant tree type by lai dominance
c
        domtree = 0
        maxlai = 0.0
c
        do 110 j = 1, 8
          if (plai(i,j).gt.maxlai) then
            domtree = j
            maxlai = plai(i,j)
          endif
 110    continue
c
c assign initial vegetation type
c
        vegtype0(i) = -999.99
c
c dominant type:  tropical broadleaf evergreen tree
c
        if (domtree.eq.1) then
          if (treelai.gt.2.5)         vegtype0(i) =  1.0  ! tropical evergreen forest / woodland
          if (treelai.le.2.5)         vegtype0(i) =  9.0  ! savanna
          if (treelai.le.0.5) then
            if (grasslai.ge.shrublai) vegtype0(i) = 10.0  ! grassland
            if (shrublai.ge.grasslai) vegtype0(i) = 11.0  ! closed shrubland
          endif
        endif
c
c dominant type:  tropical broadleaf drought-deciduous tree
c
        if (domtree.eq.2) then
          if (treelai.gt.2.5)         vegtype0(i) =  2.0  ! tropical deciduous forest / woodland
          if (treelai.le.2.5)         vegtype0(i) =  9.0  ! savanna
          if (treelai.le.0.5) then
            if (grasslai.ge.shrublai) vegtype0(i) = 10.0  ! grassland
            if (shrublai.ge.grasslai) vegtype0(i) = 11.0  ! closed shrubland
          endif
        endif
c
c dominant type:  warm-temperate broadleaf evergreen tree
c
        if (domtree.eq.3) then
          if (treelai.gt.2.5)         vegtype0(i) =  3.0  ! temperate evergreen broadleaf forest / woodland
          if (treelai.le.2.5)         vegtype0(i) =  9.0  ! savanna
          if (treelai.le.0.5) then
            if (grasslai.ge.shrublai) vegtype0(i) = 10.0  ! grassland
            if (shrublai.ge.grasslai) vegtype0(i) = 11.0  ! closed shrubland
          endif
        endif
c
c dominant type:  temperate conifer evergreen tree
c
        if (domtree.eq.4) then
          if (treelai.gt.1.5)         vegtype0(i) =  4.0  ! temperate evergreen conifer forest / woodland
          if (treelai.le.1.5)         vegtype0(i) =  9.0  ! savanna
          if (treelai.le.0.5) then
            if (grasslai.ge.shrublai) vegtype0(i) = 10.0  ! grassland
            if (shrublai.ge.grasslai) vegtype0(i) = 11.0  ! closed shrubland
          endif
        endif
c
c dominant type:  temperate broadleaf deciduous tree
c
        if (domtree.eq.5) then
          if (treelai.gt.1.5)         vegtype0(i) =  5.0  ! temperate deciduous forest / woodland
          if (treelai.le.1.5)         vegtype0(i) =  9.0  ! savanna
          if (treelai.le.0.5) then
            if (grasslai.ge.shrublai) vegtype0(i) = 10.0  ! grassland
            if (shrublai.ge.grasslai) vegtype0(i) = 11.0  ! closed shrubland
          endif
        endif
c
c dominant type:  boreal conifer evergreen tree
c
        if (domtree.eq.6)             vegtype0(i) =  6.0  ! boreal evergreen forest / woodland
c
c       if (domtree.eq.6) then
c         if (treelai.gt.1.0)         vegtype0(i) =  6.0  ! boreal evergreen forest / woodland
c         if (treelai.le.1.0) then
c           if (grasslai.ge.shrublai) vegtype0(i) = 10.0  ! grassland
c           if (shrublai.ge.grasslai) vegtype0(i) = 11.0  ! closed shrubland
c         endif
c       endif
c
c dominant type:  boreal broadleaf cold-deciduous tree
c
        if (domtree.eq.7)             vegtype0(i) =  7.0  ! boreal deciduous forest / woodland
c
c       if (domtree.eq.7) then
c         if (treelai.gt.1.0)         vegtype0(i) =  7.0  ! boreal deciduous forest / woodland
c         if (treelai.le.1.0) then
c           if (grasslai.ge.shrublai) vegtype0(i) = 10.0  ! grassland
c           if (shrublai.ge.grasslai) vegtype0(i) = 11.0  ! closed shrubland
c         endif
c       endif
c
c dominant type:  boreal conifer cold-deciduous tree
c
        if (domtree.eq.8)             vegtype0(i) =  7.0  ! boreal deciduous forest / woodland
c
c       if (domtree.eq.8) then
c         if (treelai.gt.1.0)         vegtype0(i) =  7.0  ! boreal deciduous forest / woodland
c         if (treelai.le.1.0) then
c           if (grasslai.ge.shrublai) vegtype0(i) = 10.0  ! grassland
c           if (shrublai.ge.grasslai) vegtype0(i) = 11.0  ! closed shrubland
c         endif
c       endif
c
c temperate/boreal forest mixtures
c
        if ((domtree.ge.4).and.(domtree.le.8)) then
          ratio = (plai(i,5) + plai(i,7) + plai(i,8)) / 
     >            (plai(i,4) + plai(i,5) + plai(i,6) + 
     >             plai(i,7) + plai(i,8))
          if (treelai.gt.1.0) then
            if ((ratio.gt.0.45).and.(ratio.lt.0.55)) vegtype0(i) = 8.
          endif
          if ((domtree.le.5).and.(treelai.le.1.0)) then
            if (grasslai.ge.shrublai) vegtype0(i) = 10.0  ! grassland
            if (shrublai.ge.grasslai) vegtype0(i) = 11.0  ! closed shrubland
          endif
        endif
c
c no tree is dominant
c
        if (domtree.eq.0) then
          if (treelai.gt.1.0)         vegtype0(i) =  9.0  ! savanna
          if (treelai.le.1.0) then
            if (grasslai.ge.shrublai) vegtype0(i) = 10.0  ! grassland
            if (shrublai.ge.grasslai) vegtype0(i) = 11.0  ! closed shrubland
          endif
        endif
c
c overriding vegtation classifications
c
        if (totlai.lt.1.0)            vegtype0(i) = 12.0  ! open shrubland
        if (totlai.le.0.4)            vegtype0(i) = 14.0  ! desert
c
c overriding climatic rules
c
        if (gdd5(i).lt.350.0) then
          if (totlai.ge.0.4)          vegtype0(i) = 13.0  ! tundra
          if (totlai.lt.0.4)          vegtype0(i) = 15.0  ! polar desert
        endif
c
        if (gdd0(i).lt.100.0)         vegtype0(i) = 15.0  ! polar desert
c
 100  continue
c
c return to the main program
c
      return
      end
c
