c
c  ####    ####      #    #
c #       #    #     #    #
c  ####   #    #     #    #
c      #  #    #     #    #
c #    #  #    #     #    #
c  ####    ####      #    ######
c
c ---------------------------------------------------------------------
      subroutine setsoi
c ---------------------------------------------------------------------
c
c sets diagnostic soil quantities
c
      use implicit
c
      use compar
      use comatm
      use comsno
      use comsoi
      use comveg
c
c Local variables
c
      integer i, k,           ! loop indices
     >        msand,          ! % of sand in grid point
     >        mclay           ! % of clay in grid point
c
      real fsand,             ! fraction of sand in grid point
     >     fsilt,             ! fraction of silt in grid point
     >     fclay,             ! fraction of clay in grid point
* MEM: added forganic for organic soils.
     >     forganic,          ! fraction of organic soil in grid point
     >     powliq,            ! liquid water content in fraction of soil depth
     >     powice,            ! ice water content in fraction of soil depth
     >     zcondry,           ! dry-soil conductivity
     >     zvap,              ! latent heat of vaporisation at soil temp
     >     zsub,              ! latent heat of sublimation at soil temp
     >     zwpud,             ! fraction of soil surface covered by puddle
*     >     zwpmax,            ! assumed maximum value of zwpud
     >     zwsoi,             ! volumetric water content of top soil layer 
     >     rwork1,
     >     rwork2
c
      use comsat    
c
c set soil layer quantities
c
      do 100 k = 1, nsoilay
c
        do 110 i = 1, npoi
c
c Convert input sand and clay percents to fractions
c
          msand = nint(sand(i,k))
          mclay = nint(clay(i,k))
c
          fsand = 0.01 * msand
          fclay = 0.01 * mclay
          fsilt = 0.01 * (100 - msand - mclay)
c
c update thermal conductivity (w m-1 k-1)
c
c based on c = c1**v1 * c2**v2 * c3**v3 * c4**v4 where c1,c2..
c are conductivities of soil grains, air, liquid and ice
c respectively, and v1,v2... are their volume fractions 
c (so v1 = 1-p where p is the porosity, and v1+v2+v3+v4 = 1).
c then condry = c1**(1-p) * c2**p  is the dry-soil
c conductivity, and c = condry * (c3/c2)**v3 * (c4/c2)**v4, 
c where c2 = conductivity of air = .025 w m-1 k-1.
c however this formula agrees better with williams+smith
c table 4 for wet (unfrozen) sand and clay if c2 is decreased
c to ~.005. (for peat in next section, ok if c2 = .025).
c also see lachenbruch etal,1982,jgr,87,9301 and refs therein.
c
          powliq = poros(i,k) * wsoi(i,k) * (1. - wisoi(i,k))
          powice = poros(i,k) * wisoi(i,k)
c
          zcondry = fsand * 0.300 +
     >              fsilt * 0.265 +
     >              fclay * 0.250 ! +
* M. El Maayar added this to account for contribution of organic soils
c     >              forganic * 0.026
c
          consoi(i,k) = zcondry * ((0.56*100.)**powliq)
     >                          * ((2.24*100.)**powice)
c
 110    continue
c
 100  continue
c
c set qglif - the fraction of soil sfc evaporation from soil liquid,
c soil ice, puddle liquid, and puddle ice (relative to total sfc evap)
c
c zwpud:   fraction of surface area covered by puddle (range: 0 - zwpmax)
c zwpmax:  maximum value of zwpud (currently assumed to be 0.5)
c 1-zwpud: fraction of surface area covered by soil (range: (1-zwpmax) - 1.0)
c zwsoi:   volumetric water content of top soil layer (range: 0 - 1.0)
c
c qglif(i,1): fraction of soil evap (fvapg) from soil liquid
c qglif(i,2): fraction of soil evap (fvapg) from soil ice
c qglif(i,3): fraction of soil evap (fvapg) from puddle liquid
c qglif(i,4): fraction of soil evap (fvapg) from puddle ice
c
      do 200 i = 1, npoi
c
*        zwpmax = 0.5
        zwpud = max (0.0, min (zwpmax, zwpmax*(wpud(i)+wipud(i))/wpudmax) )
        zwsoi = (1. - wisoi(i,1)) * wsoi(i,1) + wisoi(i,1)
c
        if (zwsoi.ge.epsilon) then
c
          rwork1 = 1./zwsoi
c
          if (zwpud.ge.epsilon) then
            rwork2 = 1./(wpud(i) + wipud(i))
            qglif(i,1) = (1. - zwpud) * (1. - wisoi(i,1)) * wsoi(i,1) * rwork1
            qglif(i,2) = (1. - zwpud) * wisoi(i,1) * rwork1
            qglif(i,3) = zwpud * wpud(i) * rwork2
            qglif(i,4) = zwpud * wipud(i) * rwork2
          else
            qglif(i,1) = (1. - wisoi(i,1)) * wsoi(i,1) * rwork1
            qglif(i,2) = wisoi(i,1) * rwork1
            qglif(i,3) = 0.0
            qglif(i,4) = 0.0
          endif
c
        else
c
c for a 100% dry soil surface, assign all soil evap to the puddles.
c Note that for small puddle sizes, this could lead to negative
c puddle depths. However, for a 100% dry soil with small puddles,
c evaporation is likely to be very small or less than zero
c (condensation), so negative puddle depths are not likely to occur.
c
          if (zwpud.ge.epsilon) then
            rwork2 = 1./(wpud(i) + wipud(i))
            qglif(i,1) = 0.0
            qglif(i,2) = 0.0
            qglif(i,3) = zwpud * wpud(i) * rwork2
            qglif(i,4) = zwpud * wipud(i) * rwork2
          else
            if (tsoi(i,1).ge.tmelt) then
c
c above freezing
c
              qglif(i,1) = 0.
              qglif(i,2) = 0.
              qglif(i,3) = 1.
              qglif(i,4) = 0.
c
            else
c
c below freezing
c
              qglif(i,1) = 0.
              qglif(i,2) = 0.
              qglif(i,3) = 0.
              qglif(i,4) = 1.
            endif
          endif
c
        endif
c
c set latent heat values
c
        zvap = hvapf (tsoi(i,1), ta(i))
        zsub = hsubf (tsoi(i,1), ta(i))
c
        hvasug(i) = (qglif(i,1) + qglif(i,3)) * zvap +
     >              (qglif(i,2) + qglif(i,4)) * zsub 
c
        hvasui(i) = hsubf(tsno(i,1),ta(i))
c
 200  continue   
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine soilctl
c ---------------------------------------------------------------------
c
c steps soil/seaice model through one timestep
c
      use implicit
c
      use compar
      use comhyd
      use comatm
      use comsno
      use comsoi
      use comveg
      use com1d
c
      integer i, k          ! loop indices
c
      real zfrez,           ! factor decreasing runoff fraction for tsoi < tmelt
     >     zrunf,           ! fraction of rain that doesn't stay in puddle (runoff fraction)
     >     rwork,           ! 
     >     wipre,           ! storing variable
     >     zdpud,           ! used to compute transfer from puddle to infiltration
     >     cx,              ! average specific heat for soil, water and ice
     >     chav,            ! average specific heat for water and ice
     >     zwsoi      
c
      real
     >  owsoi(npoi,nsoilay),    ! old value of wsoi
     >  otsoi(npoi,nsoilay),    ! old value of tsoi
     >  c0pud(npoi,nsoilay),    ! layer heat capacity due to puddles (=0 except for top)
     >  c1pud(npoi,nsoilay),    ! updated av. specifilayer heat capacity due to  puddle
     >  wflo(npoi,nsoilay+1),   ! = drainage at the bottom, returned by soilh2o
     >  fwtop(npoi),            ! evaporation rate from soil (for soilh2o)
     >  fhtop(npoi),            ! heat flux through soil surface (for soilheat)
     >  fwpud(npoi),            ! portion of puddle that infiltrates in soil (rate)
     >  fsqueez(npoi),          ! excess amount of water (soilh2o) 
     >  dh(npoi),               ! correction if water at tsoi < tmelt or ice at temp > tmelt
     >  dw(npoi),               ! '
     >  zporos(npoi)            
c
      use comsat
c
      call const (c0pud, npoi*nsoilay, 0.0)
      call const (c1pud, npoi*nsoilay, 0.0)
c
c for soil, set soil infiltration rate fwtop (for 
c soilh2o) and upper heat flux fhtop (for soilheat)
c
c also step puddle model wpud, wipud
c
c procedure is:
c
c   (0) immediately transfer any excess puddle liq to runoff
c
c   (1) apportion raing btwn puddle liquid(wpud) or runoff(grunof)
c
c   (2) apportion evap/condens (fvapg) btwn infil rate(fwtop), soil
c       ice (wisoi(i,1)), puddle liq (wpud), or puddle ice (wipud)
c
c   (3) transfer some puddle liquid to fwtop
c
c   (4) compute upper heat flx fhtop: includes fwtop*ch2o*tsoi(i,1)
c       to be consistent with whflo in soilheat, and accounts for
c       changing rain temp from traing to tsoi(i,1) and runoff temp
c       from tsoi to max(tsoi(i,1),tmelt)
c
      do 100 i = 1, npoi
c
c (0) immediately transfer any excess puddle liq to runoff
c
c the following runoff formulation could give rise to very
c small amounts of negative runoff
c
        grunof(i) = min (wpud(i), max (0., wpud(i) + wipud(i) -
     >                                     wpudmax)) / dtime
c
        wpud(i) = wpud(i) - grunof(i) * dtime
c
c (1) apportion sfc-level rain between puddle liquid and runoff
c
c linear dependence of runoff fraction on wpud+wipud assumes
c uniform statistical distribution of sub-grid puddle 
c capacities between 0 and wpudmax. runoff fraction is 
c reduced linearly for tsoi < tmelt (by zfrez) in which case
c any rain will increase wpud which will be frozen to wipud
c below
c
        zfrez = max (0., min (1., (tsoi(i,1) - tmelt + .5) * 2.))
c
c always have some minimal amount of runoff (0.10) even if 
c puddles are dry or soil is cold
c
        zrunf = zfrez * max (0., min (1., (wpud(i) + wipud(i)) / 
     >                                     wpudmax))
c
        wpud(i) = wpud(i) + (1. - zrunf) * raing(i) * dtime
c
        grunof(i) = grunof(i) + zrunf * raing(i)
c
c (2) apportion evaporation or condensation between 4 h2o stores:
c
        rwork = fvapg(i) * dtime
c
        if (fvapg(i).ge.0.) then
c
c evaporation: split according to qglif
c
          fwtop(i)   =            - qglif(i,1)*fvapg(i)
          wpud(i)    = wpud(i)    - qglif(i,3)*rwork
          wipud(i)   = wipud(i)   - qglif(i,4)*rwork
c
          wipre = wisoi(i,1)
          wisoi(i,1) = max (0., wipre - qglif(i,2)*rwork /
     >                          (rhow*poros(i,1)*hsoi(1)))
c
          if (1.-wisoi(i,1).gt.epsilon)
     >      wsoi(i,1) = wsoi(i,1)*(1.-wipre)/(1.-wisoi(i,1))
c
        else
c
c condensation: give all to puddles (to avoid wsoi, wisoi > 1)
c
          fwtop(i) = 0.
          wpud(i) = wpud(i)  - (qglif(i,1)+qglif(i,3))*rwork
          wipud(i)= wipud(i) - (qglif(i,2)+qglif(i,4))*rwork
c
        endif
c
c (3) transfer some puddle liquid to infiltration; can lead
c     to small amounts of negative wpud (in soilh2o) due to
c     round-off error
c
        zdpud = rhow * dtime * max (0., 1.-wisoi(i,1))**2 *
     >          hydraul(i,1)
c
        fwpud(i) = max (0., min (wpud(i), zdpud)) / dtime
        c0pud(i,1) = ch2o*wpud(i) + cice*wipud(i)
c
c (4) compute upper soil heat flux
c
        fhtop(i) = heatg(i)
     >           + raing(i)*ch2o*(traing(i)-tsoi(i,1))
     >           - grunof(i)*ch2o*max(tmelt-tsoi(i,1), 0.)
c
c update diagnostic variables
c
        gadjust(i) = 0.0
c
 100  continue
c
c reduce soil moisture due to transpiration (upsoi[u,l], from
c turvap).need to do that before other time stepping below since 
c specific heat of this transport is neglected
c
c first set porosflo, reduced porosity due to ice content, used
c as the effective porosity for uptake here and liquid hydraulics
c later in soilh2o. to avoid divide-by-zeros, use small epsilon
c limit; this will always cancel with epsilon or 0 in numerators
c
c also increment soil temperature to balance transpired water
c differential between temps of soil and leaf. physically
c should apply this to the tree, but would be awkward in turvap.
c 
c also, save old soil moisture owsoi and temperatures otsoi so
c implicit soilh2o and soilheat can aposteriori deduce fluxes.
c
      do 120 k = 1, nsoilay
c
        do 130 i = 1, npoi
c
          porosflo(i,k) = poros(i,k) * max (epsilon, (1.-wisoi(i,k)))
c
c next line just for ice whose poros(i,k) is 0.0
c
          porosflo(i,k) = max (porosflo(i,k), epsilon)
c
          wsoi(i,k) = wsoi(i,k) - dtime * (upsoiu(i,k) + upsoil(i,k)) /
     >                            (rhow * porosflo(i,k) * hsoi(k))
c
          cx = c0pud(i,k) + 
     >         (   (1.-poros(i,k))*csoi(i,k)*rhosoi(i,k)
     >           + poros(i,k)*(1.-wisoi(i,k))*wsoi(i,k)*ch2o*rhow
     >           + poros(i,k)*wisoi(i,k)*cice*rhow
     >         ) * hsoi(k)
c
          tsoi(i,k) = tsoi(i,k) - dtime * ch2o * 
     >                (  upsoiu(i,k)*(tu(i)-tsoi(i,k))
     >                 + upsoil(i,k)*(tl(i)-tsoi(i,k)) ) / cx
c
          owsoi(i,k)  = wsoi(i,k)
          otsoi(i,k)  = tsoi(i,k)
c
 130    continue
c
 120  continue
c
c step soil moisture calculations
c
      call soilh2o (owsoi, fwtop, fwpud, fsqueez, wflo)
c
c update drainage and puddle
c
      do 200 i = 1, npoi
c
        gdrain(i)  = wflo(i,nsoilay+1)
        c1pud(i,1) = ch2o*wpud(i) + cice*wipud(i)
c
 200  continue
c
c step temperatures due to conductive heat transport
c
      call soilheat (otsoi, owsoi, c0pud, fhtop, wflo, c1pud)
c
c set wsoi, wisoi to exactly 0 or 1 if differ by negligible 
c amount (needed to avoid epsilon errors in loop 400 below)
c
      call wadjust
c
c heat-conserving adjustment for liquid/ice below/above melt
c point. uses exactly the same logic as for intercepted veg h2o
c in steph2o2. we employ the fiction here that soil liquid and
c soil ice both have density rhow, to avoid "pot-hole"
c difficulties of expansion on freezing. this is done by 
c dividing all eqns through by rhow(*hsoi).
c
c the factor (1-wsoi(old))/(1-wisoi(new)) in the wsoi increments
c results simply from conservation of h2o mass; recall wsoi is
c liquid content relative to ice-reduced pore space.
c
      do 400 k = 1, nsoilay
        do 410 i = 1, npoi
c
c next line is just to avoid divide-by-zero for ice with
c poros = 0
c
          zporos(i) = max (poros(i,k), epsilon)
          rwork = c1pud(i,k)/rhow/hsoi(k)
     >           + (1.-zporos(i))*csoi(i,k)*rhosoi(i,k)/rhow
c
          chav = rwork
     >           + zporos(i)*(1.-wisoi(i,k))*wsoi(i,k)*ch2o
     >           + zporos(i)*wisoi(i,k)*cice
c
c if liquid exists below melt point, freeze some to ice
c
c (note that if tsoi>tmelt or wsoi=0, nothing changes.)
c (also note if resulting wisoi=1, either dw=0 and prev
c wisoi=1, or prev wsoi=1, so use of epsilon is ok.)
c
          zwsoi = min (1., wsoi(i,k))
c
          dh(i) = chav * (tmelt-tsoi(i,k))
          dw(i) = min ( zporos(i)*(1.-wisoi(i,k))*zwsoi,
     >                  max (0.,dh(i)/hfus) )
c
          wisoi(i,k) = wisoi(i,k) +  dw(i)/zporos(i)
          wsoi(i,k)  = wsoi(i,k)  - (dw(i)/zporos(i))*(1.-zwsoi)
     >                              / max (epsilon,1.-wisoi(i,k))
c
          chav = rwork
     >           + zporos(i)*(1.-wisoi(i,k))*wsoi(i,k)*ch2o
     >           + zporos(i)*wisoi(i,k)*cice
c
          tsoi(i,k) = tmelt - (dh(i)-hfus*dw(i)) / chav
c
c if ice exists above melt point, melt some to liquid
c
c note that if tsoi<tmelt or wisoi=0, nothing changes
c
c also note if resulting wisoi=1, dw=0 and prev wisoi=1,
c so use of epsilon is ok
c
          dh(i) = chav * (tsoi(i,k) - tmelt)
          dw(i) = min ( zporos(i)*wisoi(i,k), max (0., dh(i)/hfus) )
c
          wisoi(i,k) = wisoi(i,k) -  dw(i)/zporos(i)
          wsoi(i,k)  = wsoi(i,k)  + (dw(i)/zporos(i))
     >                 * (1.-wsoi(i,k)) / max(epsilon,1.-wisoi(i,k))
c
          chav = rwork
     >           + zporos(i)*(1.-wisoi(i,k))*wsoi(i,k)*ch2o 
     >           + zporos(i)*wisoi(i,k)*cice
c
          tsoi(i,k) = tmelt + (dh(i)-hfus*dw(i)) / chav
c
c reset porosflo (although not used after this)
c
          porosflo(i,k) = zporos(i) * max (epsilon, 1.-wisoi(i,k))
c
  410   continue
  400 continue
c 
c set wsoi, wisoi to exactly 0 or 1 if differ by negligible 
c amount (roundoff error in loop 400 above can produce very
c small negative amounts)
c
      call wadjust
c
c repeat ice/liquid adjustment for upper-layer puddles (don't 
c divide through by rhow*hsoi). upper-layer moistures wsoi,wisoi
c are already consistent with tsoi(i,1) > or < tmelt, and will 
c remain consistent here since tsoi(i,1) will not cross tmelt
c
      k = 1
c
      do 500 i = 1, npoi
c
c if any puddle liquid below tmelt, freeze some to puddle ice
c
        rwork = ( (1.-poros(i,k))*csoi(i,k)*rhosoi(i,k)
     >           + poros(i,k)*(1.-wisoi(i,k))*wsoi(i,k)*ch2o*rhow
     >           + poros(i,k)*wisoi(i,k)*cice*rhow
     >         ) * hsoi(k)
c
        chav = ch2o*wpud(i) + cice*wipud(i) + rwork
c
        dh(i) = chav * (tmelt-tsoi(i,k))
        dw(i) = min (wpud(i), max (0., dh(i)/hfus))
        wipud(i) = wipud(i) + dw(i)
        wpud(i)  = wpud(i)  - dw(i)
        chav = ch2o*wpud(i) + cice*wipud(i) + rwork
        tsoi(i,k) = tmelt - (dh(i)-hfus*dw(i)) / chav
c
c if any puddle ice above tmelt, melt some to puddle liquid
c
        dh(i) = chav * (tsoi(i,k)-tmelt)
        dw(i) = min (wipud(i), max (0., dh(i)/hfus))
        wipud(i) = wipud(i) - dw(i)
        wpud(i)  = wpud(i)  + dw(i)
        chav = ch2o*wpud(i) + cice*wipud(i) + rwork
        tsoi(i,k) = tmelt + (dh(i)-hfus*dw(i)) / chav
c  
  500 continue
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine soilh2o (owsoi, fwtop, fwpud, fsqueez, wflo)
c ---------------------------------------------------------------------
c
c sets up call to tridia to solve implicit soil moisture eqn,
c using soil temperatures in wsoi (in comsoi)
c
c lower bc can be no h2o flow or free drainage, set by bperm below
c
      use implicit
c
      use compar
      use comsoi
      use com1d
c
c Arguments : all are supplied except wflo (returned):
c
      real owsoi(npoi,nsoilay),  ! soil moistures at start of timestep
     >     fwtop(npoi),          ! evaporation from top soil layer
     >     fwpud(npoi),          ! h2o flux into top layer (infiltration from puddle)
     >     fsqueez(npoi),        ! excess water at end of time step in soil column  
     >     wflo(npoi,nsoilay+1)  ! downward h2o flow across layer boundaries
c
c local variables
c     
      integer k, i,    ! loop indices
     >        km1,
     >        kka, 
     >        kkb 
c
      integer 
     >  m(npoi),                ! exponents 
     >  n(npoi)
c
      real dmin,                ! minimum diffusivity for dry soils (m**2 s-1) 
     >     rimp,                ! implicit fraction of the calculation (0 to 1)
*     >     bperm,               ! = 0 for impermeable (no drainage) lower bc,
     >                          ! = 1 for free drainage lower bc
     >     zbex, z, dt, zz
c 
      real hsoim(nsoilay+1)  ! vertical distances between centers of layers
c
      real
     >  wsoim(npoi,nsoilay+1),    ! interpolated moisture values at layer boundaries
     >  wsoia(npoi,nsoilay+1),    ! 
     >  wsoib(npoi,nsoilay+1),    ! '
     >  weim(npoi,nsoilay+1),     ! '
     >  weip(npoi,nsoilay+1),     ! '
     >  a(npoi),                  ! intermediate terms (const for each pt)
     >  b(npoi),                  ! 
     >  bwn(npoi),                ! 
     >  bwn1(npoi),               ! 
     >  e(npoi,nsoilay+1),        ! intermediate terms in algebraic devel 
     >  f(npoi,nsoilay+1),        ! '
     >  g(npoi,nsoilay+1),        ! '
     >  d1(npoi,nsoilay),         ! diagonals of tridiagonal systems of equations 
     >  d2(npoi,nsoilay),         !  '
     >  d3(npoi,nsoilay),         !  '
     >  rhs(npoi,nsoilay),        ! right-hand sides of systems of equations
     >  w1(npoi,nsoilay),         ! work arrays needed by tridia
     >  w2(npoi,nsoilay)          !  '
c
      save dmin, rimp
      data dmin, rimp /1.e-9, 1.0/
c
c set lower boundary condition for the soil
c (permeability of the base)
c
c     bperm = 0.00  ! e.g. fully impermeable base
c     bperm = 1.00  ! e.g. fully permeable base
c
*      bperm = 0.10
c
c set level vertical distances, interpolated moistures, and
c interpolation weights
c
c top layer
c
      k = 1
c
      do 100 i = 1, npoi
c
        hsoim(k) = 0.5 * hsoi(k)
c
        weip(i,k) = 0.0
        weim(i,k) = 1.0
c
        wsoim(i,k) = wsoi(i,k)
        wsoia(i,k) = min (wsoim(i,k), 1.0)
        wsoib(i,k) = min (wsoim(i,k), 1.0)
c
  100 continue
c
c middle layers
c
      do 110 k = 2, nsoilay
c
        do 120 i = 1, npoi
c
          hsoim(k) = 0.5 * (hsoi(k-1) + hsoi(k))
c
          weip(i,k) = 0.5 * hsoi(k) / hsoim(k) 
          weim(i,k) = 1.0 - weip(i,k)
c
          wsoim(i,k) = weim(i,k) * wsoi(i,k-1) + weip(i,k) * wsoi(i,k)
          wsoia(i,k) = min (wsoim(i,k), 1.0)
          wsoib(i,k) = min (wsoim(i,k), 1.0)
c
  120   continue
c
  110 continue
c
c bottom layer
c
      k = nsoilay + 1
c
      do 130 i = 1, npoi
c
        hsoim(k) = 0.5 * hsoi(k-1)
c
        weip(i,k) = 1.0
        weim(i,k) = 0.0
c
        wsoim(i,k) = wsoi(i,k-1)
        wsoia(i,k) = min (wsoim(i,k), 1.0)
        wsoib(i,k) = min (wsoim(i,k), 1.0)
c
  130 continue
c
c set intermediate quantities e,f,g. these are terms in the
c expressions for the fluxes at boundaries between layers,
c so are zero for k=1. use bwn1 to account for minimum 
c diffusivity dmin. bperm is used for k=nsoilay+1 to set the
c type of the lower bc.
c
c top layer
c
      k = 1
c
      call const (e(1,k), npoi, 0.0)
      call const (f(1,k), npoi, 0.0)
      call const (g(1,k), npoi, 0.0)
c
c middle layers
c
      do 200 k = 2, nsoilay
c
        do 210 i = 1, npoi
c
c now that hydraul, suction and ibex can vary with depth,
c use averages of surrounding mid-layer values
c
c (see notes 8/27/93)
c
          a(i) = weim(i,k) * hydraul(i,k-1) +
     >           weip(i,k) * hydraul(i,k  )
c
          b(i) = weim(i,k) * hydraul(i,k-1) * 
     >           suction(i,k-1) * bex(i,k-1) +
     >           weip(i,k) * hydraul(i,k  ) *
     >           suction(i,k  ) * bex(i,k  )
c
          zbex = weim(i,k) * bex(i,k-1) +
     >           weip(i,k) * bex(i,k  ) 
c
          m(i) = 2 * nint(zbex) + 3
          n(i) =     nint(zbex) + 2
c
          bwn1(i) = b(i) * (wsoib(i,k)**(n(i)-1))
          bwn(i)  = bwn1(i) * wsoib(i,k)
c
          if (bwn(i).lt.dmin) bwn1(i) = 0.0
          bwn(i) = max (bwn(i), dmin)
c
          e(i,k) =  (-1.+rimp*m(i))*a(i)*(wsoia(i,k)**m(i))
     >            + ((1.-rimp)*bwn(i) - rimp*n(i)*bwn1(i)*wsoib(i,k))
     >              * (wsoi(i,k)-wsoi(i,k-1)) / hsoim(k)
c
          f(i,k) = - rimp*m(i)*a(i)*(wsoia(i,k)**(m(i)-1))
     >             + rimp*n(i)*bwn1(i)
     >               * (wsoi(i,k)-wsoi(i,k-1)) / hsoim(k)
c
          g(i,k) = rimp*bwn(i)
c
  210     continue
c
  200 continue
c
c bottom layer
c
      k = nsoilay + 1
c
      do 220 i = 1, npoi
c
        a(i) = hydraul(i,nsoilay) 
        b(i) = hydraul(i,nsoilay)*suction(i,nsoilay)*ibex(i,nsoilay)
c
        m(i) = 2*ibex(i,nsoilay) + 3
        n(i) = ibex(i,nsoilay)   + 2
c
        e(i,k) =                -a(i)*(wsoia(i,k)**m(i))*bperm
        f(i,k) = 0.0
        g(i,k) = 0.0
c
  220 continue
c
c deduce all e,f,g in proportion to the minimum of the two 
c adjacent layers' (1-wisoi), to account for restriction of flow
c by soil ice. this will cancel in loop 300  with the factor 
c 1-wisoi in (one of) the layer's porosflo, even if wisoi=1 by 
c the use of epsilon limit. so a layer with wisoi=1 will form a 
c barrier to flow of liquid, but still have a predicted wsoi
c
      do 230 k = 1, nsoilay+1
c
        kka = max (k-1,1)
        kkb = min (k,nsoilay)
c
        do 240 i=1,npoi
c
c multiply by an additional factor of 1-wisoi for stability
c
          z = max(0.,1.-max(wisoi(i,kka),wisoi(i,kkb)))**2
c
          e(i,k) = z * e(i,k)
          f(i,k) = z * f(i,k)
          g(i,k) = z * g(i,k)
c
  240   continue
c
  230 continue
c
c set matrix diagonals and right-hand sides
c
      do 300 k = 1, nsoilay
c
        do 310 i = 1, npoi
c
          dt = dtime / (porosflo(i,k)*hsoi(k))
          d1(i,k) = dt*(   f(i,k)*0.5*hsoi(k)/hsoim(k)
     >                   - g(i,k)/hsoim(k) )
          rhs(i,k) = wsoi(i,k) + dt*( e(i,k+1) - e(i,k) )
c
  310   continue
c
        if (k.eq.1) then
c
          do 320 i=1,npoi
c
            dt = dtime / (porosflo(i,k)*hsoi(k))
            rhs(i,k) = rhs(i,k) + dt*(fwtop(i)+fwpud(i))/rhow
c
  320     continue
c
        endif
c
        if (k.lt.nsoilay) then
c
          km1 = max (k-1,1)
c
          do 330 i=1,npoi
c
            dt = dtime / (porosflo(i,k)*hsoi(k))
            d2(i,k) = 1. + dt*( - f(i,k+1)*0.5*hsoi(k+1)/hsoim(k+1)
     >                          + f(i,k)  *0.5*hsoi(km1)/hsoim(k)
     >                          + g(i,k+1)/hsoim(k+1)
     >                          + g(i,k)  /hsoim(k) )
            d3(i,k) = dt*( - f(i,k+1)*0.5*hsoi(k)/hsoim(k+1)
     >                     - g(i,k+1)            /hsoim(k+1) )
c
  330     continue
c
        else if (k.eq.nsoilay) then
c
          do 340 i=1,npoi
c
            dt = dtime / (porosflo(i,k)*hsoi(k))
            d2(i,k) = 1. + dt*( - f(i,k+1)
     >                          + f(i,k)  *0.5*hsoi(k-1)/hsoim(k)
     >                          + g(i,k)  /hsoim(k) )
            d3(i,k) = 0.0
c
  340     continue
c
        endif
c
  300 continue
c
c solve the systems of equations
c
      call tridia (npoi, npoi, nsoilay, d1,d2,d3, rhs, wsoi, w1,w2)
c
      do 400 i = 1, npoi
c
        fsqueez(i) = 0.0
        wflo(i,nsoilay+1) = - rhow * e(i,nsoilay+1)
c
  400 continue
c
      do 500 k = nsoilay, 1, -1
c
        do 510 i = 1, npoi
c
          zz = rhow * poros(i,k) * 
     >         max(epsilon, (1.-wisoi(i,k))) * hsoi(k)
c
          wsoi(i,k) = wsoi(i,k) + dtime * fsqueez(i) / zz 
          fsqueez(i) = max (wsoi(i,k)-1.,0.) * zz / dtime
          wsoi(i,k) = min (wsoi(i,k),1.)
c
          wflo(i,k) = wflo(i,k+1) + (wsoi(i,k)-owsoi(i,k)) * zz / dtime
c
  510   continue
c
  500 continue
c
c step puddle liquid due to fsqueez and fwpud
c
c also subtract net puddle-to-top-layer flux from wflo(i,1),
c since puddle and top soil layer are lumped together in soilheat
c so upper wflo should be external flux only (evap/condens)

      do 600 i = 1, npoi
c
        wpud(i)   = wpud(i)   + (fsqueez(i) - fwpud(i)) * dtime
        wflo(i,1) = wflo(i,1) + (fsqueez(i) - fwpud(i))
c
  600 continue
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine soilheat (otsoi, owsoi, c0pud, fhtop, wflo, c1pud)
c ---------------------------------------------------------------------
c
c        sets up call to tridia to solve implicit soil/ice heat 
c        conduction, using layer temperatures in tsoi (in comsoi).
c        the heat flux due to liquid flow previously calculated
c        in soilh2o is accounted for. lower bc is conductive flux = 0
c        for soil (although the flux due to liquid drainage flow can
c        be > 0)
c
      use implicit
c
      use compar
      use comsoi
c
c Arguments
c
      real
     >  otsoi(npoi,nsoilay),    ! soil/ice temperatures at start of timestep (redundant
c                                 with tsoi, but passed to be consistent with soilh2o)
     >  owsoi(npoi,nsoilay),    ! soil moistures at start of timestep (before soilh2o)
     >  c0pud(npoi,nsoilay),    ! layer heat capacity due to puddles (=0 except for top)
     >  c1pud(npoi,nsoilay),    ! updated c0pud
     >  fhtop(npoi),            ! heat flux into top layer from atmos
     >  wflo(npoi,nsoilay+1)    ! downward h2o  flow across layer boundaries
c
c local variables
c
      integer k, i, km1, kp1    ! loop indices
c
      real rimp,                ! implicit fraction of the calculation (0 to 1)
     >     rwork, rwork1,       ! work variables
     >     rwork2, rwork3,      ! work variables
     >     t
c
      real
     >  whflo(npoi,nsoilay+1),   ! downward heat fluxes across layer bdries due to h2o
     >                           ! movement calculated in soilh2o
     >  con(npoi,nsoilay+1),     ! conduction coeffs between layers
     >  c0(npoi,nsoilay),        ! specific heats at start of timestep
     >  c1(npoi,nsoilay),        ! specific heats at end of timestep
     >  d1(npoi,nsoilay),        ! diagonals of tridiagonal systems of equations
     >  d2(npoi,nsoilay),        !  ''
     >  d3(npoi,nsoilay),        !  ''
     >  rhs(npoi,nsoilay),       ! right-hand sides of systems of equations
     >  w1(npoi,nsoilay),        ! work arrays needed by tridia
     >  w2(npoi,nsoilay)
c
      data rimp /1.0/
c
c set conduction coefficient between layers, and heat fluxes
c due to liquid transport
c
c top layer
c
      k = 1
c
      do 100 i = 1, npoi
        con(i,k) = 0.0
        whflo(i,k) = wflo(i,k) * ch2o * tsoi(i,k)
  100 continue
c
c middle layers
c
      do 110 k = 2, nsoilay
c
        do 120 i = 1, npoi
c
          con(i,k) =  1. / (0.5 * (hsoi(k-1) / consoi(i,k-1) +
     >                             hsoi(k)   / consoi(i,k)))
c
          t = (hsoi(k) * tsoi(i,k-1) + hsoi(k-1) * tsoi(i,k)) /
     >        (hsoi(k-1)             + hsoi(k))
c
          whflo(i,k) = wflo(i,k) * ch2o * t
c
  120   continue
c
  110 continue
c
c bottom layer
c
      k = nsoilay + 1
c
      do 130 i = 1, npoi
        con(i,k) = 0.0
        whflo(i,k) = wflo(i,k) * ch2o * tsoi(i,k-1)
  130 continue
c
c set diagonals of matrix and right-hand side. use old and
c new heat capacities c0, c1 consistently with moisture fluxes
c whflo computed above, to conserve heat associated with 
c changing h2o amounts in each layer
c
      do 200 k = 1, nsoilay
c
        km1 = max (k-1,1)
        kp1 = min (k+1,nsoilay)
c
        do 210 i=1,npoi
c
          rwork1 = (1.-poros(i,k))*csoi(i,k)*rhosoi(i,k)
          rwork2 = poros(i,k)*(1.-wisoi(i,k))*ch2o*rhow
          rwork3 = poros(i,k)*wisoi(i,k)*cice*rhow
c
          c0(i,k) = c0pud(i,k) + 
     >              (   rwork1
     >                + rwork2 * owsoi(i,k)
     >                + rwork3
     >              ) * hsoi(k)
c
          c1(i,k) = c1pud(i,k) +  
     >              (   rwork1
     >                + rwork2 * wsoi(i,k)
     >                + rwork3
     >              ) * hsoi(k)
c
          rwork = dtime/c1(i,k)
c
          d1(i,k) =    - rwork * rimp * con(i,k)
          d2(i,k) = 1. + rwork * rimp * (con(i,k)+con(i,k+1))
          d3(i,k) =    - rwork * rimp * con(i,k+1)
c
          rhs(i,k) = (c0(i,k)/c1(i,k))*tsoi(i,k) + rwork
     >               * ( (1.-rimp)*con(i,k)  *(tsoi(i,km1)-tsoi(i,k))
     >                 + (1.-rimp)*con(i,k+1)*(tsoi(i,kp1)-tsoi(i,k))
     >                 + whflo(i,k) - whflo(i,k+1) )
c
  210   continue
c
        if (k.eq.1) then
          do 220 i=1,npoi
            rhs(i,k) = rhs(i,k) + (dtime/c1(i,k))*fhtop(i)
  220     continue
        endif
c
  200 continue
c
c solve systems of equations
c
      call tridia (npoi, npoi, nsoilay, d1,d2,d3, rhs, tsoi, w1,w2)
c
c deduce downward heat fluxes between layers
c
      call scopy (npoi, fhtop, hflo(1,1))
c
      do 300 k=1,nsoilay
        do 310 i=1,npoi
          hflo(i,k+1) = hflo(i,k) -
     >                  (c1(i,k)*tsoi(i,k) - c0(i,k)*otsoi(i,k)) / dtime
  310   continue
  300 continue
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine wadjust
c ---------------------------------------------------------------------
c
c set wsoi, wisoi to exactly 0 if differ by negligible amount, 
c to protect epsilon logic in soilctl and soilh2o
c
c ice-liquid transformations in soilctl loop 400 can produce very
c small -ve amounts due to roundoff error, and very small -ve or +ve
c amounts can cause (harmless) "underflow" fpes in soilh2o
c
      use implicit
c
      use compar
      use comhyd
      use comsoi
      use com1d
c 
c local variables
c
      integer k, i
      real ztot0, ztot1
c
      do 100 k = 1, nsoilay
        do 110 i = 1, npoi
c
c initial total soil water
c
         ztot0 = hsoi(k) * poros(i,k) * rhow *
     >           ((1. - wisoi(i,k)) * wsoi(i,k) + wisoi(i,k))
c
c set bounds on wsoi and wisoi
c
         if (wsoi(i,k).lt.epsilon)  wsoi(i,k)  = 0.0
         if (wisoi(i,k).lt.epsilon) wisoi(i,k) = 0.0
c
         wsoi(i,k)  = min (1., wsoi(i,k))
         wisoi(i,k) = min (1., wisoi(i,k))
c
         if (wisoi(i,k).ge.1-epsilon) wsoi(i,k) = 0.0
c
c for diagnosis of total adjustment
c
         ztot1 = hsoi(k) * poros(i,k) * rhow *
     >           ((1. - wisoi(i,k)) * wsoi(i,k) + wisoi(i,k))
c
         gadjust(i) = gadjust(i) + (ztot1 - ztot0) / dtime
c
  110   continue
  100 continue
c
      return
      end
c
