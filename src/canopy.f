c
c  ####     ##    #    #   ####   #####    #   #
c #    #   #  #   ##   #  #    #  #    #    # #
c #       #    #  # #  #  #    #  #    #     #
c #       ######  #  # #  #    #  #####      #
c #    #  #    #  #   ##  #    #  #          #
c  ####   #    #  #    #   ####   #          #
c
c ---------------------------------------------------------------------
      subroutine canopy
c ---------------------------------------------------------------------
c
c calculates sensible heat and moisture flux coefficients,
c and steps canopy temperatures through one timestep
c
c atmospheric conditions at za are supplied in comatm
c arrays ta, qa, psurf and scalar siga (p/ps)
c
c downward sensible heat and moisture fluxes at za
c are returned in com1d arrays fsena, fvapa
c
c sensible heat and moisture fluxes from solid objects to air
c are stored (for other models and budget) in com1d arrays
c fsen[u,s,l,g,i], fvap[u,s,l,g,i]
c
c the procedure is first to compute wind speeds and aerodynamic
c transfer coefficients in turcof, then call turvap to solve an
c implicit linear system for temperatures and specific
c humidities and the corresponding fluxes - this is iterated
c niter times for non-linearities due to stratification,
c implicit/explicit (h2o phase), dew, vpd and max soil
c moisture uptake - t12 and q12 are changed each iteration,
c and tu, ts, tl, tg, ti can be adjusted too
c
c initialize aerodynamic quantities
c
       include 'implicit.h'
c
c Local variables
c
      integer niter,       ! total number of ierations
     >        iter         ! number of iteration

      call canini
c
c estimate soil moisture stress parameters
c
      call drystress
c
c iterate the whole canopy physics solution niter times:
c
      niter = 3
c
      do 100 iter = 1, niter
c
c calculate wind speeds and aerodynamic transfer coeffs
c
        call turcof (iter)
c
c calculate canopy photosynthesis rates and conductance
c
        call stomata
c
c solve implicit system of heat and water balance equations
c
        call turvap (iter, niter)
c
  100 continue
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine canini
c ---------------------------------------------------------------------
c
c initializes aerodynamic quantities that remain constant 
c through one timestep
c
c note that some quantities actually are
c constant as long as the vegetation amounts and fractional
c coverage remain unchanged, so could re-arrange code for
c efficiency - currently all arrays initialized here are in
c com1d which can be overwritten elsewhere
c
c rwork is used throughout as a scratch variable to reduce number of
c computations
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'comatm.h'
      include 'comsno.h'
      include 'comsoi.h'
      include 'comveg.h'
      include 'com1d.h'
c
c Local variables
c
      integer i       ! loop indice
c
      real siga,      ! sigma level of atmospheric data
     >     pa,        ! pressure at level of atmospheric data
     >     x,         ! density of vegetation (without distinction between
     >                ! lai,sai)
     >     x1,        ! density of vegetation (different max)
     >     rwork,     ! difference between top and bottom of canopy 
     >     cvegl,     !
     >     dvegl,     ! diffusion coefficient for lower canopy
     >     bvegl,     ! e-folding depth in canopy for lower canopy
     >     cvegu,     !
     >     dvegu,     ! diffusion coefficient for upper canopy
     >     bvegu      ! e-folding depth in canopy for upper canopy
c
c define sigma level of atmospheric data
c
c currently, the value of siga is set to 0.999. This is roughly 10 meters
c above ground, which is the typical height for the CRU05 input wind speed data
c
      siga = 0.999
c
      tfac = 1.0 / (siga**cappa)
c
c atmospheric conditions at za
c za is variable, although siga = p/ps is constant
c
      do 100 i = 1, npoi
c
        pa = psurf(i) * siga
c
        rhoa(i) = pa / ( rair * ta(i) * 
     >            (1.0 + (rvap / rair - 1.0) * qa(i)) )
c
        cp(i) = cair * (1.0 + (cvap / cair - 1.0) * qa(i))
c
        za(i) = (psurf(i) - pa) / (rhoa(i) * grav)
c
c make sure that atmospheric level is higher than canopy top
c
        za(i) = max (za(i), ztop(i,2) + 1.0)
c
 100  continue 
c
c aerodynamic coefficients for the lower story
c
c cvegl (drag coeff for momentum) is proportional, and dvegl
c (diffusion coeff for momentum) inversely proportional,
c to x = density of vegetation (without distinction between
c lai,sai and fl*(1-fi)) - x is not allowed to be exactly
c zero to avoid divide-by-zeros, and for x>1 dvegl is 
c proportional to 1/x**2 so that roughness length tends to
c zero as x tends to infinity
c
c also the top, bottom and displacement heights z3(i),z4(i),
c displ(i) tend to particular values as the density tends to
c zero, to give same results as equations for no veg at all.
c
      do 200 i = 1, npoi
c
        x = fl(i) * (1.0 - fi(i)) * 2.0 * (lai(i,1) + sai(i,1)) / alaiml
c
        x  = min (x, 3.0)
        x1 = min (x, 1.0)
c
        rwork = max(ztop(i,1)-zbot(i,1),0.01)
        cvegl = (0.4 / rwork) *
     >           max(1.e-5, x)
c
        dvegl = (0.1 * rwork) / 
     >           max(1.e-5, x, x**2)
c
c e-folding depth in canopy
c
        bvegl = sqrt (2.0 * cvegl / dvegl )
c
c [(tau/rho)/u**2] for inf canopy
c
        bdl(i) = 0.5 * bvegl * dvegl
c
c 1 / diffusion coefficient
c
        dil(i) = 1. / dvegl
c
        rwork = (1.0 - x1) * (max (z0soi(i),z0sno) + 0.01) 
c
        z3(i) = x1 * ztop(i,1) + rwork
c
        z4(i) = x1 * zbot(i,1) + rwork
c
        z34(i) = 0.5 * (z3(i) + z4(i))
c
        exphl(i) = exp (0.5 * bvegl * (z3(i)-z4(i)))
        expl(i)  = exphl(i)**2
c
        displ(i) = x1 * 0.7 * z3(i)
c
 200  continue 
c
c aerodynamic coefficients for the upper story
c same comments as for lower story
c
      do 300 i = 1, npoi
c
        x = fu(i) * 2.0 * (lai(i,2)+sai(i,2)) / alaimu
c
        x  = min (x, 3.0)
        x1 = min (x, 1.0)
c
        rwork = max(ztop(i,2)-zbot(i,2),.01)
        cvegu = (0.4 / rwork) * 
     >           max(1.e-5,x)
c
        dvegu = (0.1 * rwork) / 
     >           max(1.e-5,x,x**2)
c
        rwork = 1. / dvegu
        bvegu  = sqrt (2.0 * cvegu * rwork)
        bdu(i) = 0.5 * bvegu * dvegu
        diu(i) = rwork
c
        rwork = (1.0 - x1) * (z3(i) + 0.01)
        z1(i) = x1 * ztop(i,2) + rwork
        z2(i) = x1 * zbot(i,2) + rwork
c
        z12(i) = 0.5 * (z1(i) + z2(i))
c
        exphu(i) = exp (0.5 * bvegu * (z1(i) - z2(i)))
        expu(i)  = exphu(i)**2
c
        dispu(i) = x1 * 0.7 * z1(i) + (1.0 - x1) * displ(i)
c
 300  continue 
c
c mixing-length logarithms
c
      do 400 i = 1, npoi
c
        alogg(i)  = alog (z0soi(i))
        alogi(i)  = alog (z0sno)
        alogav(i) = (1.0 - fi(i)) * alogg(i) + fi(i) * alogi(i)
c
c alog4 must be > z0soi, z0sno to avoid possible problems later 
c
        alog4(i) = alog ( max (z4(i), 1.1*z0soi(i), 1.1*z0sno) )
        alog3(i) = alog (z3(i)-displ(i))
        alog2(i) = alog (z2(i)-displ(i))
        alog1(i) = alog (z1(i)-dispu(i))
        aloga(i) = alog (za(i)-dispu(i))
c
c initialize u2, alogu, alogl for first iteration's fstrat
c
        u2(i)    = ua(i)/exphu(i)
        alogu(i) = alog (max(.01, .1*(z1(i)-z2(i))))
        alogl(i) = alog (max(.01, .1*(z3(i)-z4(i))))
c
  400 continue
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine turcof (iter)
c ---------------------------------------------------------------------
c
c solves for wind speeds at various levels
c
c also computes upper and lower-region air-air transfer coefficients
c and saves them in com1d arrays cu and cl for use by turvap,
c and similarly for the solid-air transfer coefficients
c su, ss, sl, sg and si
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'comatm.h'
      include 'comsoi.h'
      include 'comsno.h'
      include 'comveg.h'
      include 'com1d.h'
c
c Arguments (input)
c
      integer iter          !current iteration number
c
c
c Local variables
c
      integer i             ! loop indice
c
      real xfac,            !
     >     x,               !
     >     rwork,           ! working variable
     >     cdmax,           ! max value for cd
     >     tauu,            !
     >     a,b,c,d,         !
     >     taul,            !
     >     ca,              ! to compute inverse air-air transfer coeffs
     >     cai, cbi, cci,   !
     >     cdi, cei, cfi,   !
     >     sg0,             ! to compute air-solid transfer coeff for soil
     >     si0              ! to compute air-solid transfer coeff for ice
c
      real yu(npoi), yl(npoi)
c
c set stratification factors for lower and upper regions
c using values from the previous iteration
c
      xfac = 1.0
c
      call fstrat (t34, t12, xfac, q34, q12, z3, z2, 
     >             alogl, alogl, alog2, u2, richl, straml, strahl, iter)
c
      call fstrat (t12, ta,  tfac, q12, qa,  z1, za, 
     >             alogu, alogu, aloga, ua, richu, stramu, strahu, iter)
c
c eliminate c/d from eq (28), tau_l/rho from (26),(27), to get
c lower-story roughness alogl. yl/bdl is (tau_l/rho)/(c+d)
c
c equation numbers correspond to lsx description section 4.e
c
      do 100 i = 1, npoi
c
        x = ((alog4(i)-alogav(i))/vonk)**2 * bdl(i)
c
        rwork = 1. / expl(i)
        yl(i) = ((x+1)*expl(i) + (x-1)*rwork)
     >        / ((x+1)*expl(i) - (x-1)*rwork)
c
        alogl(i) = alog3(i) - vonk * sqrt(yl(i)/bdl(i))
c
 100  continue 
c
c eliminate tau_l/rho from (24),(25), tau_u/rho and a/b from
c (22),(23), to get upper-story roughness alogu
c 
c yu/bdu is (tau_u/rho)/(a+b)
c
      do 110 i = 1, npoi
c          
        x = ((alog2(i)-alogl(i))/vonk)**2 * bdu(i) / straml(i)
c
        rwork = 1. / expu(i)
        yu(i) = ((x+1)*expu(i) + (x-1)*rwork)
     >        / ((x+1)*expu(i) - (x-1)*rwork)
c
        alogu(i) = alog1(i) - vonk * sqrt(yu(i)/bdu(i))
c
 110  continue
c
c define the maximum value of cd
c
      cdmax = 300.0 / (2.0 * dtime)
c
c get tauu (=tau_u/rho) from (21), a and b from (22),(23),
c taul (=tau_u/rho) from (25), c and d from (26),(27)
c
c changed the following to eliminate small errors associated with
c moving this code to single precision - affected c and d,
c which made u_ become undefined, as well as affecting some
c other variables
c
      do 200 i = 1, npoi
c
        tauu = (ua(i) * vonk/(aloga(i)-alogu(i)))**2 * stramu(i)
c
        a = 0.5 * tauu * (yu(i)+1)/bdu(i)
        b = 0.5 * tauu * (yu(i)-1)/bdu(i)
c
        taul = bdu(i) * (a/expu(i) - b*expu(i))
c
        c = 0.5 * taul * (yl(i)+1)/bdl(i)
        d = 0.5 * taul * (yl(i)-1)/bdl(i)
c
c evaluate wind speeds at various levels, keeping a minimum 
c wind speed of 0.01 m/s at all levels
c   
        u1(i)  = max (0.01, sqrt (max (0.0, (a+b))))
        u12(i) = max (0.01, sqrt (max (0.0, (a/exphu(i)+b*exphu(i)))))
        u2(i)  = max (0.01, sqrt (max (0.0, (a/expu(i) +b*expu(i)))))
        u3(i)  = max (0.01, sqrt (max (0.0, (c+d))))
        u34(i) = max (0.01, sqrt (max (0.0, (c/exphl(i)+d*exphl(i)))))
        u4(i)  = max (0.01, sqrt (max (0.0, (c/expl(i) +d*expl(i)))))
c
 200  continue
c
c compute inverse air-air transfer coeffs
c
c use of inverse individual coeffs cai, cbi, cci, cdi, cei, cfi avoids
c divide-by-zero as vegetation vanishes - combine into
c upper-region coeff cu from za to z12, and lower-region coeff
c cl from z34 to z12, and also coeffs
c
      do 300 i = 1, npoi
c
        ca = ua(i)*strahu(i)*vonk**2  /
     >       ((aloga(i)-alogu(i)) * (aloga(i)-alog1(i)))
c
        ca = min (cdmax, ca / (1. + ca * 1.0e-20))
c
        cai = 1.0 / (rhoa(i)*ca)
c
        cbi = diu(i) * (z1(i)-z12(i)) / (rhoa(i) * 0.5*(u1(i)+u12(i)))
        cci = diu(i) * (z12(i)-z2(i)) / (rhoa(i) * 0.5*(u12(i)+u2(i)))
c
        cdi = (alog2(i)-alogl(i)) * (alog2(i)-alog3(i)) /
     >        (rhoa(i)*u2(i)*strahl(i)*vonk**2)
c
        cei = dil(i) * (z3(i)-z34(i)) / (rhoa(i) * 0.5*(u3(i)+u34(i)))
        cfi = dil(i) * (z34(i)-z4(i)) / (rhoa(i) * 0.5*(u34(i)+u4(i)))
c
        cu(i) = 1.0 / (cai + cbi)
        cl(i) = 1.0 / (cci + cdi + cei)
c
c compute air-solid transfer coeffs for upper leaves, upper
c stems, lower story (su,ss,sl)
c
        su(i) = rhoa(i) * cleaf  * sqrt (u12(i) / dleaf(2))
        ss(i) = rhoa(i) * cstem  * sqrt (u12(i) / dstem(2))
        sl(i) = rhoa(i) * cgrass * sqrt (u34(i) / dleaf(1))
c
c compute air-solid transfer coeffs for soil and snow (sg,si)
c
c old technique
c
c       sg0 = rhoa(i) * u4(i) * (vonk/(alog4(i)-alogg(i)))**2
c       si0 = rhoa(i) * u4(i) * (vonk/(alog4(i)-alogi(i)))**2
c
c replace above formulations which depend on the log-wind profile
c (which may not work well below a canopy), with empirical formulation
c of Norman's. In the original LSX, turcof.f solves for the winds at
c the various levels from the momentum equations. This gives the transfer
c coefficients for heat and moisture. Heat and moisture eqns are then solved 
c in subroutine turvap. Using the empirical formulation of John Norman is 
c not consistent with the earlier solution for u4 (based on a logarithmic 
c profile just above the ground. However, this is used here because it 
c improved a lot simulations of the sensible heat flux over the 
c HAPEX-MOBILHY and FIFE sites
c
        sg0 = rhoa(i) * (0.004 + 0.012 * u4(i))
        si0 = rhoa(i) * (0.003 + 0.010 * u4(i))
c
c modify the cofficient to deal with cfi (see above)
c
        sg(i) = 1.0 / (cfi + 1.0 / sg0)
        si(i) = 1.0 / (cfi + 1.0 / si0)
c
  300 continue
c
c JAF:  not necessary 
c
c if no veg, recalculate coefficients appropriately for a
c single logarithmic profile, and 2 fictitious levels just
c above soil/snow surface. these levels are arbitrary but are
c taken as z2 and z4, preset in vegdat to a few cm height
c for bare ground and ice. use strahu from above, which used
c t12 and alogu (ok after first iteration)
c
c     do 600 i = 1, npoi
c
c       if ((fu(i).eq.0.0).and.(fl(i).eq.0.0)) then
c
c         z = rhoa(i)*ua(i)*strahu(i)*vonk**2 / (aloga(i)-alogav(i))
c
c         ca    = z / (aloga(i)-alog2(i))
c         cu(i) = rhoa(i)*min (cdmax,
c    >                          ca / (1. + ca / 1.0e+20))
c
c         cl(i) = z / (alog2(i)-alog4(i))
c
c         sg(i) = z / (alog4(i)-alogg(i))
c         si(i) = z / (alog4(i)-alogi(i))
c
c         alogu(i) = alogav(i)
c
c       endif
c
c 600 continue
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine turvap (iter, niter)
c ---------------------------------------------------------------------
c
c solves canopy system with linearized implicit sensible heat and
c moisture fluxes
c
c first, assembles matrix arr of coeffs in linearized equations
c for tu,ts,tl,t12,t34,q12,q34,tg,ti and assembles the right hand
c sides in the rhs vector
c
c then calls linsolve to solve this system, passing template mplate of
c zeros of arr 
c 
c finally calculates the implied fluxes and stores them 
c for the agcm, soil, snow models and budget calcs
c
c common blocks
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'comhyd.h'
      include 'comatm.h'
      include 'comsno.h'
      include 'comsoi.h'
      include 'comveg.h'
      include 'com1d.h'
c
c Arguments (input)
c
      integer niter,      ! total # of iteration
     >        iter        ! # of iteration
c
c local variables
c
      integer i, k
c
      real rwork, zwtot, rwork2, tgav, tiav, tuav, 
     >     tsav, tlav, quav, qsav, qlav, qgav, qiav, zwpud, zwsoi,
     >     psig, hfac, hfac2, zwopt, zwdry, betaw, emisoil, e, qs1,
     >     dqs1, xnumer, xdenom, betafac, betas
c
      real
     >  xu(npoi),
     >  xs(npoi),
     >  xl(npoi),
     >  chux(npoi), 
     >  chsx(npoi), 
     >  chlx(npoi),
     >  chgx(npoi),
     >  wlgx(npoi),
     >  wigx(npoi),
     >  fradu(npoi), 
     >  frads(npoi), 
     >  fradl(npoi),
     >  wu(npoi),
     >  ws(npoi),
     >  wl(npoi),      
     >  wg(npoi),
     >  wi(npoi), 
     >  qu(npoi),
     >  qs(npoi),  
     >  ql(npoi),  
     >  qg(npoi),  
     >  qi(npoi),
     >  dqu(npoi),
     >  dqs(npoi), 
     >  dql(npoi), 
     >  dqg(npoi), 
     >  dqi(npoi),
     >  tuold(npoi),
     >  tsold(npoi),
     >  tlold(npoi),
     >  tgold(npoi), 
     >  tiold(npoi),
     >  tupre(npoi),
     >  tspre(npoi),
     >  tlpre(npoi),
     >  tgpre(npoi),
     >  tipre(npoi),
     >  suw(npoi),
     >  ssw(npoi),
     >  slw(npoi),
     >  sut(npoi),
     >  slt(npoi),
     >  slt0(npoi),
     >  suh(npoi),
     >  ssh(npoi), 
     >  slh(npoi),
     >  cog(npoi),
     >  coi(npoi),
     >  zirg(npoi),
     >  ziri(npoi),
     >  qgfac(npoi),
     >  qgfac0(npoi)
c
      save 
     >  xu,
     >  xs,
     >  xl,
     >  chux, 
     >  chsx, 
     >  chlx,
     >  chgx,
     >  wlgx,
     >  wigx,
     >  cog,
     >  coi,
     >  zirg,
     >  ziri,
     >  wu,  
     >  ws,    
     >  wl,    
     >  wg, 
     >  wi,
     >  tuold, 
     >  tsold, 
     >  tlold, 
     >  tgold, 
     >  tiold
c
      integer nqn
c
      parameter (nqn=9)
c
      real arr(npoi,nqn,nqn),      !    
     >     rhs(npoi,nqn),          ! right hand side
     >     vec(npoi,nqn)           ! 
c
      integer  mplate(nqn,nqn)
cc
c                  tu  ts  tl t12 t34 q12 q34  tg  ti 
c                  ----------------------------------
      data mplate / 1,  0,  0,  1,  0,  1,  0,  0,  0, !tu
     >              0,  1,  0,  1,  0,  1,  0,  0,  0, !ts
     >              0,  0,  1,  0,  1,  0,  1,  0,  0, !tl
     >              1,  1,  0,  1,  1,  0,  0,  0,  0, !t12
     >              0,  0,  1,  1,  1,  0,  0,  1,  1, !t34
     >              1,  1,  0,  0,  0,  1,  1,  0,  0, !q12
     >              0,  0,  1,  0,  0,  1,  1,  1,  1, !q34
     >              0,  0,  0,  0,  1,  0,  1,  1,  0, !tg
     >              0,  0,  0,  0,  1,  0,  1,  0,  1  !ti
     >            /
c
      include 'comsat.h'
c
c if first iteration, save original canopy temps in t*old
c (can use tsoi,tsno for original soil/snow skin temps), for
c rhs heat capacity terms in matrix soln, and for adjustment
c of canopy temps after each iteration
c
c also initialize soil/snow skin temps tg, ti to top-layer temps
c
c the variables t12, t34, q12, q34, for the first iteration
c are saved via global arrays from the previous gcm timestep,
c this is worth doing only if the agcm forcing is
c smoothly varying from timestep to timestep
c
      if (iter.eq.1) then
c
c weights for canopy coverages
c
        do 10 i = 1, npoi
          xu(i) = 2.0 * lai(i,2) * fu(i)
          xs(i) = 2.0 * sai(i,2) * fu(i)
          xl(i) = 2.0 * (lai(i,1) + sai(i,1)) * fl(i) * (1.0 - fi(i))
 10     continue 
c
c specific heats per leaf/stem area
c
        do 20 i = 1, npoi
          chux(i) = chu + ch2o * wliqu(i) + cice * wsnou(i)
          chsx(i) = chs + ch2o * wliqs(i) + cice * wsnos(i)
          chlx(i) = chl + ch2o * wliql(i) + cice * wsnol(i)
 20     continue
c
        do 30 i = 1, npoi 
c
          rwork = poros(i,1) * rhow
c
          chgx(i) = ch2o * wpud(i) + cice * wipud(i)
     >              + ((1.-poros(i,1))*csoi(i,1)*rhosoi(i,1)
     >              + rwork*(1.-wisoi(i,1))*wsoi(i,1)*ch2o
     >              + rwork*wisoi(i,1)*cice
     >                ) * hsoi(1)
c
          wlgx(i) = wpud(i) +
     >              rwork * (1. - wisoi(i,1)) *
     >              wsoi(i,1) * hsoi(1)
c
          wigx(i) = wipud(i) + rwork * wisoi(i,1) * hsoi(1)
c
 30     continue 
c
c conductivity coeffs between ground skin and first layer
c
        do 40 i = 1, npoi
          cog(i) = consoi(i,1) / (0.5 * hsoi(1))
          coi(i) = consno      / (0.5 * max (hsno(i,1), hsnotop))
 40     continue
c
c d(ir emitted) / dt for soil
c
        rwork = 4. * 0.95 * stef
c
        do 50 i = 1, npoi
          zirg(i) = rwork * (tg(i)**3)
          ziri(i) = rwork * (ti(i)**3)
 50     continue
c
c updated temperature memory
c
        do 60 i = 1, npoi
          tuold(i) = tu(i)
          tsold(i) = ts(i)
          tlold(i) = tl(i)
          tgold(i) = tg(i)
          tiold(i) = ti(i)
 60     continue
c
      endif
c
c set implicit/explicit factors w* (0 to 1) for this iteration
c w* is 1 for fully implicit, 0 for fully explicit
c for first iteration, impexp and impexp2 set w* to 1
c
      call impexp (wu, tu, chux, wliqu, wsnou, iter)
      call impexp (ws, ts, chsx, wliqs, wsnos, iter)
      call impexp (wl, tl, chlx, wliql, wsnol, iter)
      call impexp (wg, tg, chgx, wlgx,  wigx,  iter)
c
c call impexp2 for snow model
c
      call impexp2 (wi, ti, tiold, iter)
c
c adjust t* for this iteration 
c
c in this routine we are free to choose them, 
c since they are just the central values about which the 
c equations are linearized - heat is conserved in the matrix
c solution because t*old are used for the rhs heat capacities
c
c here, let t* represent the previous soln if it was fully
c implicit, but weight towards t*old depending on the amount
c (1-w*) the previous soln was explicit
c
c this weighting is necessary for melting/freezing surfaces, for which t*
c is kept at t*old, presumably at or near tmelt
c
      do 80 i = 1, npoi
        tu(i) = wu(i) * tu(i) + (1.0 - wu(i)) * tuold(i)
        ts(i) = ws(i) * ts(i) + (1.0 - ws(i)) * tsold(i)
        tl(i) = wl(i) * tl(i) + (1.0 - wl(i)) * tlold(i)
        tg(i) = wg(i) * tg(i) + (1.0 - wg(i)) * tgold(i)
        ti(i) = wi(i) * ti(i) + (1.0 - wi(i)) * tiold(i)
 80   continue 
c
c save current "central" values for final flux calculations
c
      do 90 i = 1, npoi
        tupre(i) = tu(i)
        tspre(i) = ts(i)
        tlpre(i) = tl(i)
        tgpre(i) = tg(i)
        tipre(i) = ti(i)
 90   continue
c
c calculate various terms occurring in the linearized eqns,
c using values of t12, t34, q12, q34 from
c the previous iteration
c
c specific humidities for canopy and ground, and derivs wrt t
c for canopy
c
c limit derivs to avoid -ve implicit q's below,
c as long as d(temp)s in one iteration are le 10 deg k
c
      do 100 i = 1, npoi
c
        e      = esat(tu(i))
        qu(i)  = qsat (e, psurf(i))
        dqu(i) = dqsat (tu(i), qu(i))
        dqu(i) = min (dqu(i), qu(i) * 0.1)
c
        e      = esat(ts(i))
        qs(i)  = qsat (e, psurf(i))
        dqs(i) = dqsat (ts(i), qs(i))
        dqs(i) = min (dqs(i), qs(i) * 0.1)
c
        e      = esat(tl(i))
        ql(i)  = qsat (e, psurf(i))
        dql(i) = dqsat (tl(i), ql(i))
        dql(i) = min (dql(i), ql(i) * 0.1)
c
        e      = esat(tg(i))
        qg(i)  = qsat (e, psurf(i))
        dqg(i) = dqsat (tg(i), qg(i))
        dqg(i) = min (dqg(i), qg(i) * 0.1)
c
        e      = esat(ti(i))
        qi(i)  = qsat (e, psurf(i))
        dqi(i) = dqsat (ti(i), qi(i))
        dqi(i) = min (dqi(i), qi(i) * 0.1)
c
 100  continue
c
c set qgfac0, factor by which soil surface specific humidity
c is less than saturation
c
c it is important to note that the qgfac expression should
c satisfy timestep cfl criterion for upper-layer soil moisture
c for small wsoi(i,1)
c
c for each iteration, qgfac is set to qgfac0, or to 1 if
c condensation onto soil is anticipated (loop 110 in canopy.f)
c
c Evaporation from bare soil is calculated using the "beta method"
c (e.g., eqns 5 & 7 of Mahfouf and Noilhan 1991, JAM 30 1354-1365),
c but converted to the "alpha method" (eqns 2 & 3 of M&N), to match
c the structure in IBIS. The conversion from the beta to alpha
c method is through the relationship:
c   alpha * qgs - q34 = beta * (hfac * qgs - q34),
c from which one solves for alpha (which is equal to qgfac0):
c   qgfac0 = alpha = (beta * hfac) + (1 - beta)*(q34/qgs)
c
        do 105 i = 1, npoi
c
c first calculate the total saturated fraction at the soil surface
c (including puddles ... see soil.f)
c
          zwpud = max (0.0, min (0.5, 0.5*(wpud(i)+wipud(i))/wpudmax) )
          zwsoi = (1.0 - wisoi(i,1)) * wsoi(i,1) + wisoi(i,1)
          zwtot = zwpud + (1. - zwpud) * zwsoi
c
c next calculate the matric potential (from eqn 9.3 of Campbell and
c Norman), multiply by gravitational acceleration to get in units
c of J/kg, and calculate the relative humidity at the soil water
c surface (i.e., within the soil matrix), based on thermodynamic
c theory (eqn 4.13 of C&N)
c
          psig = -grav * suction(i,1) * (zwtot ** (-bex(i,1)))
          hfac = exp(psig/(rvap*tg(i)))
c
c then calculate the relative humidity of the air (relative to
c saturation at the soil temperature). Note that if hfac2 > 1
c (which would imply condensation), then qgfac is set to 1
c later in the code (to allow condensation to proceed at the
c "potential rate")
c
          hfac2 = q34(i)/qg(i)
c
c set the "beta" factor and then calculate "alpha" (i.e., qgfac0)
c as the beta-weighted average of the soil water RH and the "air RH"
c First calculate beta_w:
c
          zwopt = 1.0
          zwdry = swilt(i,1)
          betaw = max(0.0, min(1., (zwtot - zwdry)/(zwopt - zwdry)) )
c
c Next convert beta_w to beta_s (see Milly 1992, JClim 5 209-226):
c
          emisoil = 0.95
          e      = esat(t34(i))
          qs1    = qsat (e, psurf(i))
          dqs1   = dqsat (t34(i), qs1)
          xnumer = hvap * dqs1
          xdenom = cp(i) + (4.0 * emisoil * stef * (t34(i))**3) / sg(i)
          betafac = xnumer / xdenom
          betas = betaw / (1.0 + betafac * (1.0 - betaw))
c
c Combine hfac and hfac2 into qgfac0 ("alpha") using beta_s
c
          qgfac0(i) = betas * hfac + (1. - betas) * hfac2
  105   continue
c
c set fractions covered by intercepted h2o to 1 if dew forms
c
c these fwet*x are used only in turvap, and are distinct from
c the real fractions fwet* that are set in fwetcal
c
c they must be exactly 1 if q12 > qu or q34 > ql, to zero transpiration
c by the factor 1-fwet[u,l]x below, so preventing "-ve" transp
c
c similarly, set qgfac, allowing for anticipated dew formation
c to avoid excessive dew formation (which then infiltrates) onto
c dry soils
c
      do 110 i = 1, npoi
c
        fwetux(i) = fwetu(i)
        if (q12(i).gt.qu(i)) fwetux(i) = 1.0
c
        fwetsx(i) = fwets(i)
        if (q12(i).gt.qs(i)) fwetsx(i) = 1.0
c
        fwetlx(i) = fwetl(i)
        if (q34(i).gt.ql(i)) fwetlx(i) = 1.0
c
        qgfac(i) = qgfac0(i)
        if (q34(i).gt.qg(i)) qgfac(i) = 1.0
c
c set net absorbed radiative fluxes for canopy components
c
        fradu(i) = 0.0
c
        if (lai(i,2).gt.epsilon)
     >     fradu(i) = (solu(i) + firu(i)) / (2.0 * lai(i,2))
c
        frads(i) = 0.0
c
        if (sai(i,2).gt.epsilon)
     >     frads(i) = (sols(i) + firs(i)) / (2.0 * sai(i,2))
c
        fradl(i) = 0.0
c
        if ((lai(i,1)+sai(i,1)).gt.epsilon)
     >     fradl(i) = (soll(i) + firl(i)) /
     >                (2.0 * (lai(i,1) + sai(i,1)))
c
 110  continue
c
c calculate canopy-air moisture transfer coeffs for wetted
c leaf/stem areas, and for dry (transpiring) leaf areas
c
c the wetted-area coeffs suw,ssw,slw are constrained to be less
c than what would evaporate 0.8 * the intercepted h2o mass in 
c this timestep (using previous iteration's q* values)
c
c this should virtually eliminate evaporation-overshoots and the need
c for the "negative intercepted h2o"  correction in steph2o2
c        
      do 200 i = 1, npoi
c
c coefficient for evaporation from wet surfaces in the upper canopy:
c
        suw(i) = min ( fwetux(i) * su(i), 
     >                 0.8 * (wliqu(i) + wsnou(i)) /
     >                 max (dtime * (qu(i) - q12(i)), epsilon))
c
c coefficient for transpiration from average upper canopy leaves:
c
        sut(i) = (1.0 - fwetux(i)) * 0.5 *
     >           ( totcondub(i) * frac(i,1) +
     >             totcondub(i) * frac(i,2) +
     >             totcondub(i) * frac(i,3) +
     >             totconduc(i) * frac(i,4) +
     >             totcondub(i) * frac(i,5) +
     >             totconduc(i) * frac(i,6) +
     >             totcondub(i) * frac(i,7) +
     >             totcondub(i) * frac(i,8) )
c
        sut(i) = max (0.0, sut(i))
c
c coefficient for sensible heat flux from upper canopy:
c
        suh(i) = suw(i) * (rliqu(i)  * hvapf(tu(i),ta(i))  +
     >                 (1.-rliqu(i)) * hsubf(tu(i),ta(i))) +
     >           sut(i) *              hvapf(tu(i),ta(i))
c
c coefficient for evaporation from wet surfaces on the stems:
c
        ssw(i) = min (fwetsx(i) * ss(i), 
     >                0.8 * (wliqs(i) + wsnos(i))
     >                / max (dtime * (qs(i) - q12(i)), epsilon))
c
c coefficient for sensible heat flux from stems:
c
        ssh(i) = ssw(i) * (rliqs(i)  * hvapf(ts(i),ta(i)) +
     >                 (1.-rliqs(i)) * hsubf(ts(i),ta(i)))
c
c coefficient for evaporation from wet surfaces in the lower canopy:
c
        slw(i) = min (fwetlx(i) * sl(i), 
     >                0.8 * (wliql(i) + wsnol(i))
     >                / max (dtime * (ql(i) - q34(i)), epsilon))
c
c coefficient for transpiration from average lower canopy leaves:
c
        slt0(i) = (1. - fwetlx(i)) * 0.5 *
     >            ( totcondls(i) * frac(i,9)  +
     >              totcondls(i) * frac(i,10) +
     >              totcondl4(i) * frac(i,11) +
     >              totcondl3(i) * frac(i,12) )
c
        slt0(i) = max (0., slt0(i))
c
c averaged over stems and lower canopy leaves:
c 
        slt(i) = slt0(i) * lai(i,1) / max (lai(i,1)+sai(i,1), epsilon)
c
c coefficient for sensible heat flux from lower canopy:
c
        slh(i) = slw(i) * (  rliql(i)  * hvapf(tl(i),ta(i))  +
     >                   (1.-rliql(i)) * hsubf(tl(i),ta(i))) +
     >           slt(i) *                hvapf(tl(i),ta(i))
c
 200  continue
c
c set the matrix of coefficients and the right-hand sides
c of the linearized equations
c
      call const(arr, npoi*nqn*nqn, 0.0)
      call const(rhs, npoi*nqn, 0.0)
c
      rwork = 1. / dtime
c
c upper leaf temperature tu
c
      do 300 i = 1, npoi
c
        rwork2 = su(i)*cp(i)
        arr(i,1,1) = chux(i)*rwork
     >             + wu(i)*rwork2
     >             + wu(i)*suh(i)*dqu(i)
        arr(i,1,4) = -rwork2
        arr(i,1,6) = -suh(i)
        rhs(i,1) = tuold(i)*chux(i)*rwork
     >           - (1.-wu(i))*rwork2*tu(i)
     >           - suh(i) * (qu(i)-wu(i)*dqu(i)*tu(i))
     >           + fradu(i) - pfluxu(i)
c 
 300  continue
c
c upper stem temperature ts
c
      do 310 i = 1, npoi
c
        rwork2 = ss(i)*cp(i)
        arr(i,2,2) = chsx(i)*rwork
     >             + ws(i)*rwork2
     >             + ws(i)*ssh(i)*dqs(i)
        arr(i,2,4) = -rwork2
        arr(i,2,6) = -ssh(i)
        rhs(i,2) = tsold(i)*chsx(i)*rwork
     >           - (1.-ws(i))*rwork2*ts(i)
     >           - ssh(i) * (qs(i)-ws(i)*dqs(i)*ts(i))
     >           + frads(i) - pfluxs(i)
c
 310  continue
c
c lower veg temperature tl
c
      do 320 i = 1, npoi
c
        rwork2 = sl(i)*cp(i)
        arr(i,3,3) = chlx(i)*rwork
     >             + wl(i)*rwork2
     >             + wl(i)*slh(i)*dql(i)
        arr(i,3,5) = -rwork2
        arr(i,3,7) = -slh(i)
        rhs(i,3) = tlold(i)*chlx(i)*rwork
     >           - (1.-wl(i))*rwork2*tl(i)
     >           - slh(i) * (ql(i)-wl(i)*dql(i)*tl(i))
     >           + fradl(i) - pfluxl(i)
c
 320  continue
c
c upper air temperature t12
c
      do 330 i = 1, npoi
c
        rwork = xu(i)*su(i)
        rwork2 = xs(i)*ss(i)
        arr(i,4,1) = -wu(i)*rwork
        arr(i,4,2) = -ws(i)*rwork2
        arr(i,4,4) = cu(i) + cl(i) + rwork + rwork2
        arr(i,4,5) = -cl(i)
        rhs(i,4) = cu(i)*ta(i)*tfac
     >           + (1.-wu(i))*rwork*tu(i)
     >           + (1.-ws(i))*rwork2*ts(i)
c
 330  continue
c
c lower air temperature t34
c
      do 340 i = 1, npoi
c
        rwork = xl(i)*sl(i)
        rwork2 = fi(i)*si(i)
        arr(i,5,3) = -wl(i)*rwork
        arr(i,5,4) = -cl(i)
        arr(i,5,5) = cl(i) + rwork
     >             + (1.-fi(i))*sg(i) + rwork2
        arr(i,5,8) = -wg(i)*(1.-fi(i))*sg(i)
        arr(i,5,9) = -wi(i)*rwork2
        rhs(i,5) = (1.-wl(i))*rwork           *tl(i)
     >           + (1.-wg(i))*(1.-fi(i))*sg(i)*tg(i)
     >           + (1.-wi(i))*rwork2          *ti(i)
c
 340  continue
c
c upper air specific humidity q12
c
      do 350 i = 1, npoi
c
        rwork = xu(i)*(suw(i)+sut(i))
        rwork2 = xs(i)*ssw(i)
        arr(i,6,1) = -wu(i)*rwork *dqu(i)
        arr(i,6,2) = -ws(i)*rwork2*dqs(i)
        arr(i,6,6) = cu(i) + cl(i)
     >             + rwork + rwork2
        arr(i,6,7) = -cl(i)
        rhs(i,6) = cu(i)*qa(i)
     >           + rwork  * (qu(i)-wu(i)*dqu(i)*tu(i))
     >           + rwork2 * (qs(i)-ws(i)*dqs(i)*ts(i))
c
  350 continue
c
c lower air specific humidity q34
c
      do 360 i = 1, npoi
c
        rwork  = xl(i)*(slw(i)+slt(i))
        rwork2 = (1.-fi(i))*sg(i)
        arr(i,7,3) = -wl(i)*rwork*dql(i)
        arr(i,7,6) = -cl(i)
        arr(i,7,7) = cl(i) + rwork
     >             + rwork2 +fi(i)*si(i)
        arr(i,7,8) = -wg(i)*rwork2*qgfac(i)*dqg(i)
        arr(i,7,9) = -wi(i)*fi(i)*si(i)*dqi(i)
        rhs(i,7)= rwork           *(ql(i)-wl(i)*dql(i)*tl(i))
     >          + rwork2*qgfac(i) *(qg(i)-wg(i)*dqg(i)*tg(i))
     >          + fi(i) *si(i)    *(qi(i)-wi(i)*dqi(i)*ti(i))
c
  360 continue
c
c soil skin temperature
c
c (there is no wg in this eqn since it solves for a fully
c implicit tg. wg can be thought of as the fractional soil
c area using a fully implicit soln, and 1-wg as that using a
c fully explicit soln. the combined soil temperature is felt
c by the lower air, so wg occurs in the t34,q34 eqns above.)
c
      do 370 i = 1, npoi
c
        rwork  = sg(i)*cp(i)
        rwork2 = sg(i)*hvasug(i)
        arr(i,8,5) = -rwork
        arr(i,8,7) = -rwork2
        arr(i,8,8) = rwork + rwork2*qgfac(i)*dqg(i)
     >             + cog(i) + zirg(i)
        rhs(i,8) = -rwork2*qgfac(i)*(qg(i)-dqg(i)*tg(i))
     >           + cog(i)*tsoi(i,1)
     >           + solg(i) + firg(i) + zirg(i) * tgold(i)
c
  370 continue
c
c snow skin temperature
c
c (there is no wi here, for the same reason as for wg above.)
c
      do 380 i = 1, npoi
c
	rwork  = si(i)*cp(i)
        rwork2 = si(i)*hvasui(i)
        arr(i,9,5) = -rwork
        arr(i,9,7) = -rwork2
        arr(i,9,9) = rwork + rwork2*dqi(i)
     >             + coi(i) + ziri(i)
        rhs(i,9) = -rwork2*(qi(i)-dqi(i)*ti(i))
     >           + coi(i)*tsno(i,1)
     >           + soli(i) + firi(i) + ziri(i) * tiold(i)
c
  380 continue
c
c solve the systems of equations
c
      call linsolve (arr, rhs, vec, mplate, nqn)
c
c copy this iteration's solution to t*, q12, q34
c
      do 400 i = 1, npoi
c
        tu(i)  = vec(i,1)
        ts(i)  = vec(i,2)
        tl(i)  = vec(i,3)
        t12(i) = vec(i,4)
        t34(i) = vec(i,5)
        tg(i)  = vec(i,8)
        ti(i)  = vec(i,9)
c
        q12(i) = vec(i,6)
        q34(i) = vec(i,7)
c
  400 continue
c
c all done except for final flux calculations,
c so loop back for the next iteration (except the last)
c
      if (iter.lt.niter) return
c
c evaluate sensible heat and moisture fluxes (per unit
c leaf/stem/snow-free/snow-covered area as appropriate)
c
c *******************************
c diagnostic sensible heat fluxes
c *******************************
c
      do 500 i = 1, npoi
c
        fsena(i) = cp(i) * cu(i) * (ta(i)*tfac - t12(i))
c
        tgav = wg(i)*tg(i) + (1.-wg(i))*tgpre(i)
        fseng(i) = cp(i) * sg(i) * (tgav - t34(i))
c
        tiav = wi(i)*ti(i) + (1.-wi(i))*tipre(i)
        fseni(i) = cp(i) * si(i) * (tiav - t34(i))

        tuav = wu(i)*tu(i) + (1. - wu(i))*tupre(i)
        fsenu(i) = cp(i) * su(i) * (tuav - t12(i))
c
        tsav = ws(i)*ts(i) + (1. - ws(i))*tspre(i)
        fsens(i) = cp(i) * ss(i) * (tsav - t12(i))
c
        tlav = wl(i)*tl(i) + (1. - wl(i))*tlpre(i)
        fsenl(i) = cp(i) * sl(i) * (tlav - t12(i))
c
 500  continue
c
c *************************
c calculate moisture fluxes
c *************************
c
      do 510 i = 1, npoi
c
c total evapotranspiration from the entire column
c
        fvapa(i)  = cu(i) * (qa(i)-q12(i))
c
c evaporation from wet surfaces in the upper canopy
c and transpiration per unit leaf area - upper canopy
c
        quav = qu(i) + wu(i)*dqu(i)*(tu(i)-tupre(i))
        fvapuw(i) = suw(i) * (quav-q12(i))
        fvaput(i) = max (0.0, sut(i) * (quav-q12(i)))
c
c evaporation from wet surfaces on stems
c
        qsav = qs(i) + ws(i)*dqs(i)*(ts(i)-tspre(i))
        fvaps(i) = ssw(i) * (qsav-q12(i))
c
c evaporation from wet surfaces in the lower canopy
c and transpiration per unit leaf area - lower canopy
c
        qlav = ql(i) + wl(i)*dql(i)*(tl(i)-tlpre(i))
        fvaplw(i) = slw(i)  * (qlav-q34(i))
        fvaplt(i) = max (0.0, slt0(i) * (qlav-q34(i)))
c
c evaporation from the ground
c
        qgav = qg(i) + wg(i)*dqg(i)*(tg(i)-tgpre(i))
        fvapg(i) = sg(i) * (qgfac(i)*qgav - q34(i))
c
c evaporation from the snow
c
        qiav = qi(i) + wi(i)*dqi(i)*(ti(i)-tipre(i))
        fvapi(i) = si(i) * (qiav-q34(i))
c
 510  continue
c 
c adjust ir fluxes
c
      do 520 i = 1, npoi
c
        firg(i) = firg(i) - wg(i)*zirg(i)*(tg(i) - tgold(i))
        firi(i) = firi(i) - wi(i)*ziri(i)*(ti(i) - tiold(i))
        firb(i) = firb(i) + (1.-fi(i))*wg(i)*zirg(i)*(tg(i)-tgold(i))
     >                    +     fi(i) *wi(i)*ziri(i)*(ti(i)-tiold(i))
c
c impose constraint on skin temperature
c
        ti(i) = min (ti(i), tmelt)
c
 520  continue
c
c set upsoi[u,l], the actual soil water uptake rates from each
c soil layer due to transpiration in the upper and lower stories,
c for the soil model 
c
      do 600 k = 1, nsoilay
        do 610 i = 1, npoi
c
          upsoiu(i,k) = fvaput(i) * 2.0 * lai(i,2) * fu(i) *
     >                  stressu(i,k) / max (stresstu(i), epsilon)
c
          upsoil(i,k) = fvaplt(i) * 2.0 * lai(i,1) * fl(i) *
     >                  (1. - fi(i)) *
     >                  stressl(i,k) / max (stresstl(i), epsilon)
c
 610    continue
 600  continue
c
c set net evaporation from intercepted water, net evaporation
c from the surface, and net transpiration rates
c
      do 700 i = 1, npoi
c
c evaporation from intercepted water
c
        ginvap(i) = fvapuw(i) * 2.0 * lai(i,2) * fu(i) +
     >              fvaps (i) * 2.0 * sai(i,2) * fu(i) +
     >              fvaplw(i) * 2.0 * (lai(i,1) + sai(i,1)) * 
     >                                 fl(i) * (1. - fi(i))
c
c evaporation from soil and snow surfaces
c
        gsuvap(i) = fvapg(i)  * (1. - fi(i)) + fvapi(i)  * fi(i)
c
c transpiration
c
        gtrans(i) = fvaput(i) * 2.0 * lai(i,2) * fu(i) +
     >              fvaplt(i) * 2.0 * lai(i,1) * fl(i) * (1.-fi(i))
c
        gtransu(i) = fvaput(i) * 2.0 * lai(i,2) * fu(i)
        gtransl(i) = fvaplt(i) * 2.0 * lai(i,1) * fl(i) * (1.-fi(i))
c
 700  continue
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine fstrat (tb, tt, ttfac, qb, qt, zb, zt, 
     >                   albm, albh, alt, u, rich, stram, strah, iter)
c ---------------------------------------------------------------------
c
c computes mixing-length stratification correction factors
c for momentum and heat/vapor, for current 1d strip, using
c parameterizations in louis (1979),blm,17,187. first computes
c richardson numbers. sets an upper limit to richardson numbers
c so lower-veg winds don't become vanishingly small in very
c stable conditions (cf, carson and richards,1978,blm,14,68)
c
c system (i) is as in louis(1979). system (vi) is improved as
c described in louis(1982), ecmwf workshop on planetary boundary
c layer parameterizations,november 1981,59-79 (qc880.4 b65w619)
c
c common blocks
c
      include 'implicit.h'
c
      include 'compar.h'
c
c input variables
c
      integer iter   ! current iteration number
c
      real ttfac     ! pot. temp factor for ttop (relative to bottom,supplied)
c
      real
     > tb(npoi),     ! bottom temperature (supplied)
     > tt(npoi),     ! top temperature (supplied)
     > qb(npoi),     ! bottom specific humidity (supplied)
     > qt(npoi),     ! top specific humidity (supplied)
     > zb(npoi),     ! height of bottom (supplied)
     > zt(npoi),     ! height of top (supplied)
     > albm(npoi),   ! log (bottom roughness length) for momentum (supplied)
     > albh(npoi),   ! log (bottom roughness length) for heat/h2o (supplied)
     > alt(npoi),    ! log (z at top) (supplied)
     > u(npoi),      ! wind speed at top (supplied)
     > rich(npoi),   ! richardson number (returned)
     > stram(npoi),  ! stratification factor for momentum (returned)
     > strah(npoi),  ! stratification factor for heat/vap (returned)
     > stramx(npoi), !
     > strahx(npoi)  !
c
c local variables
c
      integer 
     > indp(npoi),   !
     > indq(npoi)    !
c
      integer i, j, np, nq
c
      real zht, zhb, xm, xh, rwork, ym, yh, z, w
c ---------------------------------------------------------------------
      np = 0
      nq = 0
c
c do for all points
c
      do 100 i = 1, npoi
c
c calculate richardson numbers
c
        zht = tt(i)*ttfac*(1.+.622*qt(i))
        zhb = tb(i)*      (1.+.622*qb(i))
c
        rich(i) = grav * max (zt(i)-zb(i), 0.)
     >            * (zht-zhb) / (0.5*(zht+zhb) * u(i)**2)
c
c bound richardson number between -2.0 (unstable) to 1.0 (stable)
c
        rich(i) = max (-2.0, min (rich(i), 1.0))
c
 100  continue
c
c set up indices for points with negative or positive ri
c
      do 110 i = 1, npoi
c
        if (rich(i).le.0.) then
          np = np + 1
          indp(np) = i
        else
          nq = nq + 1
          indq(nq) = i
        endif
c
  110 continue
c
c calculate momentum and heat/vapor factors for negative ri
c
      if (np.gt.0) then
c
        do 200 j = 1, np
c
          i = indp(j)
c
          xm = max (alt(i)-albm(i), .5)
          xh = max (alt(i)-albh(i), .5)
c
          rwork = sqrt(-rich(i))
c
          ym = (vonk/xm)**2 * exp (0.5*xm) * rwork
          yh = (vonk/xh)**2 * exp (0.5*xh) * rwork
c
c system (vi)
c
          stramx(i) =   1.0 - 2*5*rich(i) / (1.0 + 75*ym)
          strahx(i) =   1.0 - 3*5*rich(i) / (1.0 + 75*yh)
c
  200   continue
c
      endif
c
c calculate momentum and heat/vapor factors for positive ri
c
      if (nq.gt.0) then
c
        do 300 j=1,nq
c
          i = indq(j)
c
c system (vi)
c
          z = sqrt(1.0 + 5 * rich(i))
c
          stramx(i) = 1.0 / (1.0 + 2*5*rich(i) / z)
          strahx(i) = 1.0 / (1.0 + 3*5*rich(i) * z)
c
  300   continue
c
      endif
c
c except for the first iteration, weight results with the
c previous iteration's values. this improves convergence by
c avoiding flip-flop between stable/unstable stratif, eg,
c with cold upper air and the lower surface being heated by
c solar radiation
c
      if (iter.eq.1) then
c
        do 400 i = 1, npoi
c
          stram(i) = stramx(i)
          strah(i) = strahx(i)
c
  400   continue
c
      else
c
        w = 0.5
c
        do 410 i = 1, npoi
c
          stram(i) = w * stramx(i) + (1.0 - w) * stram(i)
          strah(i) = w * strahx(i) + (1.0 - w) * strah(i)
c
  410   continue
c
      endif
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine impexp (wimp, tveg, ch, wliq, wsno, iter)
c ---------------------------------------------------------------------
c
c sets the implicit vs explicit fraction in turvap calcs for
c upper leaves, upper stems or lower veg. this is to account for
c temperatures of freezing/melting intercepted h2o constrained
c at the melt point. if a purely implicit calc is used for such
c a surface, the predicted temperature would be nearly the atmos
c equil temp with little sensible heat input, so the amount of
c freezing or melting is underestimated. however, if a purely
c explicit calc is used with only a small amount of intercepted
c h2o, the heat exchange can melt/freeze all the h2o and cause
c an unrealistic huge change in the veg temp. the algorithm
c below attempts to avoid both pitfalls
c
c common blocks
c
      include 'implicit.h'
      include 'compar.h'
c
c input/output variables
c
      integer iter  ! current iteration number (supplied)
c
      real      
     >  wimp(npoi), ! implicit/explicit fraction (0 to 1) (returned)
     >  tveg(npoi), ! temperature of veg (previous iteration's soln) (supp)
     >  ch(npoi),   ! heat capacity of veg (supplied)
     >  wliq(npoi), ! veg intercepted liquid (supplied)
     >  wsno(npoi)  ! veg intercepted snow (supplied)
c
c local variables
c
      integer i
c
      real h, z, winew
c
c for first iteration, set wimp to fully implicit, and return
c
      if (iter.eq.1) then
        call const(wimp, npoi, 1.0)
        return
      endif
c
c for second and subsequent iterations, estimate wimp based on
c the previous iterations's wimp and its resulting tveg.
c
c calculate h, the "overshoot" heat available to melt any snow
c or freeze any liquid. then the explicit fraction is taken to
c be the ratio of h to the existing h2o's latent heat (ie, 100%
c explicit calculation if not all of the h2o would be melted or
c frozen). so winew, the implicit amount, is 1 - that ratio.
c but since we are using the previous iteration's t* results
c for the next iteration, to ensure convergence we need to damp
c the returned estimate wimp by averaging winew with the 
c previous estimate. this works reasonably well even with a
c small number of iterations (3), since for instance with large
c amounts of h2o so that wimp should be 0., a good amount of 
c h2o is melted or frozen with wimp = .25
c
      do 100 i = 1, npoi
c
        h = ch(i) * (tveg(i) - tmelt)
        z = max (abs(h), epsilon)
c
        winew = 1.0
c
        if (h.gt.epsilon)  winew = 1. - min (1., hfus * wsno(i) / z)
        if (h.lt.-epsilon) winew = 1. - min (1., hfus * wliq(i) / z)
c
        wimp(i) = 0.5 * (wimp(i) + winew)
c
  100 continue
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine impexp2 (wimp, t, told, iter)
c ---------------------------------------------------------------------
c
c sets the implicit vs explicit fraction in turvap calcs for
c seaice or snow skin temperatures, to account for temperatures
c of freezing/melting surfaces being constrained at the melt
c point
c
c unlike impexp, don't have to allow for all h2o 
c vanishing within the timestep
c
c wimp   = implicit fraction (0 to 1) (returned)
c
      include 'implicit.h'
      include 'compar.h'
c
c input variables
c
      integer iter
      real 
     >  wimp(npoi), t(npoi), told(npoi)
c
c local variables
c
      integer i    ! loop indice
c
c for first iteration, set wimp to fully implicit, and return
c
      if (iter.eq.1) then
        call const(wimp, npoi, 1.0)
        return
      endif
c
      do 100 i = 1, npoi
c
        if ((t(i)-told(i)).gt.epsilon) wimp(i) = (tmelt - told(i)) / 
     >                                           (t(i)  - told(i))
        wimp(i) = max (0.0, min (1.0, wimp(i)))
c
 100  continue
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine fwetcal
c ---------------------------------------------------------------------
c
c calculates fwet[u,s,l], the fractional areas wetted by 
c intercepted h2o (liquid and snow combined) -  the maximum value
c fmax (<1) allows some transpiration even in soaked conditions
c
c use a linear relation between fwet* and wliq*,wsno* (at least
c for small values), so that the implied "thickness" is constant
c (equal to wliq*max, wsno*max as below) and the typical amount
c evaporated in one timestep in steph2o will not make wliq*,wsno*
c negative and thus cause a spurious unrecoverable h2o loss
c
c (the max(w*max,.01) below numericaly allows w*max = 0 without
c blowup.) in fact evaporation in one timestep *does* sometimes
c exceed wliq*max (currently 1 kg/m2), so there is an additional
c safeguard in turvap that limits the wetted-area aerodynamic
c coefficients suw,ssw,slw -- if that too fails, there is an 
c ad-hoc adjustment in steph2o2 to reset negative wliq*,wsno*
c amounts to zero by taking some water vapor from the atmosphere.
c
c also sets rliq[u,s,l], the proportion of fwet[u,s,l] due to
c liquid alone. fwet,rliq are used in turvap, rliq in steph2o. 
c (so rliq*fwet, (1-rliq)*fwet are the fractional areas wetted
c by liquid and snow individually.) if fwet is 0, choose rliq
c = 1 if t[u,s,l] ge tmelt or 0 otherwize, for use by turvap and
c steph2o in case of initial dew formation on dry surface.
c 
      include 'implicit.h'
c
      include 'compar.h'
      include 'comveg.h'
      include 'com1d.h'
c
c local variables
c
      integer i           ! loop indice
c
      real fmax,          ! maximum water cover on two-sided leaf
     >     xliq,          ! fraction of wetted leaf (liquid only)
     >     xtot           ! fraction of wetted leaf (liquid and snow)
c
c maximum water cover on two-sided leaf
c
      parameter (fmax = 0.25)
c
c upper leaves
c
      do 100 i = 1, npoi
c
        xliq = wliqu(i) / max (wliqumax, 0.01)
        xtot = xliq + wsnou(i) / max (wsnoumax, 0.01)
c
        fwetu(i) = min (fmax, xtot)
        rliqu(i) = xliq / max (xtot, epsilon)
c
        if (fwetu(i).eq.0.0) then
          rliqu(i) = 1.0
          if (tu(i).lt.tmelt) rliqu(i) = 0.0
        endif
c
  100 continue
c
c upper stems
c
      do 200 i = 1, npoi
c
        xliq = wliqs(i) / max (wliqsmax, 0.01)
        xtot = xliq + wsnos(i) / max (wsnosmax, 0.01)
c
        fwets(i) = min (fmax, xtot)
        rliqs(i) = xliq / max (xtot, epsilon)
c
        if (fwets(i).eq.0.0) then
          rliqs(i) = 1.0
          if (ts(i).lt.tmelt) rliqs(i) = 0.0
        endif
c
  200 continue
c
c lower veg
c
      do 300 i = 1, npoi
c
        xliq = wliql(i) / max (wliqlmax, 0.01)
        xtot = xliq + wsnol(i) / max (wsnolmax, 0.01)
c
        fwetl(i) = min (fmax, xtot)
        rliql(i) = xliq / max (xtot, epsilon)
c
        if (fwetl(i).eq.0.) then
          rliql(i) = 1.0
          if (tl(i).lt.tmelt) rliql(i) = 0.0
        endif
c
  300 continue
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine cascade
c ---------------------------------------------------------------------
c
c steps intercepted h2o due to drip, precip, and min/max limits
c
c calls steph2o for upper leaves, upper stems and lower veg in
c iurn, adjusting precips at each level
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'comatm.h'
      include 'comveg.h'
      include 'com1d.h'
c
c local variables
c
      integer i            ! loop indice
c
      real twet3           ! Function: wet bulb temperature (K)
      real twetbulb        ! wet bulb temperature (K)
c    
      real
     >  xai(npoi),         !lai and/or sai for veg component
                           ! (allows steph2o to work on any veg component)
     >  rain(npoi),        !rainfall at appropriate level (modified by steph2o)
     >  train(npoi),       !temperature of rain (modified by steph2o)  
     >  snow(npoi),        !snowfall at appropriate level (modified by steph2o)
     >  tsnow(npoi),       !temperature of snow (modified by steph2o)
     >  x1(npoi),        ! 
     >  x2(npoi),        ! 
     >  x3(npoi),        ! 
     >  x4(npoi)         ! 
c
c adjust rainfall and snowfall rates at above-tree level
c
c set wliqmin, wsnomin -- unlike wliq*max, wsno*max, these are
c part of the lsx numerical method and not from the vegetation
c database, and they are the same for all veg components
c
c the value 0.0010 should be small compared to typical precip rates
c times dtime to allow any intercepted h2o to be initiated, but
c not too small to allow evap rates to reduce wliq*, wsno* to
c that value in a reasonable number of time steps
c
      wliqmin = 0.0010 * (dtime/3600.) * (wliqumax / 0.2)
      wsnomin = 0.0010 * (dtime/3600.) * (wsnoumax / 2.0)
c
      do 50 i=1,npoi
        rainu(i) = raina(i)
c
c set rain temperature to the wet bulb temperature
c
        if (ta(i) .gt. tmelt) then
           twetbulb = twet3( ta(i), qa(i), psurf(i) )
        else
           twetbulb = tmelt
        endif
        trainu(i) = max (twetbulb, tmelt)
        x1(i) = 0.0
        x2(i) = max (t12(i), tmelt)
   50 continue
c
      call mix (rainu,trainu, rainu,trainu, x1,x2, vzero,vzero)
c
      do 52 i=1,npoi
        snowu(i) = snowa(i)
        tsnowu(i) = min (ta(i), tmelt)
        x1(i) = 0.0
        x2(i) = min (t12(i), tmelt)
   52 continue
c
      call mix (snowu,tsnowu, snowu,tsnowu, x1,x2, vzero,vzero)
c
c set up for upper leaves
c
      do 100 i = 1, npoi
        xai(i)   = 2.0 * lai(i,2)
        rain(i)  = rainu(i)
        train(i) = trainu(i)
        snow(i)  = snowu(i)
        tsnow(i) = tsnowu(i)
  100 continue
c
c step upper leaves
c
      call steph2o
     >  (tu,  wliqu,  wsnou,  xai,  pfluxu,  rain, train, snow, tsnow,
     >   tdripu, tblowu, wliqumax, wsnoumax, wliqmin, wsnomin)
c
c set up for upper stems
c the upper stems get precip as modified by the upper leaves
c
      do 200 i=1,npoi
        xai(i) = 2.0 * sai(i,2)
  200 continue
c
c step upper stems
c
      call steph2o
     >  (ts,  wliqs,  wsnos,  xai,  pfluxs,  rain, train, snow, tsnow,
     >   tdrips, tblows, wliqsmax, wsnosmax, wliqmin, wsnomin)
c
c adjust rainfall and snowfall rates at below-tree level
c allowing for upper-veg interception/drip/belowoff
c
      do 300 i=1,npoi
        x1(i) = fu(i)*rain(i)
        x2(i) = (1.-fu(i))*rainu(i)
        x3(i) = 0.0
        x4(i) = max (t34(i), tmelt)
  300 continue
c
      call mix (rainl,trainl, x1,train, x2,trainu, x3,x4)
c
      do 310 i=1,npoi
        x1(i) = fu(i)*snow(i)
        x2(i) = (1.-fu(i))*snowu(i)
        x3(i) = 0.0
        x4(i) = min (t34(i), tmelt)
  310 continue
c
      call mix (snowl,tsnowl, x1,tsnow, x2,tsnowu, x3,x4)
c
c set up for lower veg
c
      do 400 i = 1, npoi
        xai(i)   = 2.0 * (lai(i,1) + sai(i,1))
        rain(i)  = rainl(i)
        train(i) = trainl(i)
        snow(i)  = snowl(i)
        tsnow(i) = tsnowl(i)
  400 continue
c
c step lower veg
c
      call steph2o
     >  (tl,  wliql,  wsnol,  xai,  pfluxl,  rain, train, snow, tsnow,
     >   tdripl, tblowl, wliqlmax, wsnolmax, wliqmin, wsnomin)
c
c adjust rainfall and  snowfall rates at soil level,
c allowing for lower-veg interception/drip/blowoff
c
      do 500 i=1,npoi
        x1(i) = fl(i) * rain(i)
        x2(i) = (1.-fl(i)) * rainl(i)
  500 continue
c
      call mix (raing,traing, x1,train, x2,trainl, vzero,vzero)
c
      do 510 i=1,npoi
        x1(i) = fl(i) * snow(i)
        x2(i) = (1.-fl(i)) * snowl(i)
  510 continue
c
      call mix (snowg,tsnowg, x1,tsnow, x2,tsnowl, vzero,vzero)
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine steph2o
     >  (tveg,  wliq,  wsno,  xai,  pflux,  rain, train, snow, tsnow,
     >   tdrip, tblow, wliqmax, wsnomax, wliqmin, wsnomin)
c ---------------------------------------------------------------------
c
c steps intercepted h2o for one canopy component (upper leaves, 
c upper stems, or lower veg) through one lsx time step, adjusting
c for h2o sensible heat and phase changes. also modifies precip
c due to interception and drip,blowoff
c
c 
c
      include 'implicit.h'
c
      include 'compar.h'
c
c Arguments (all arguments are supplied (unchanged) unless otherwise noted
c
      real tdrip,       ! e-folding time of liquid drip  tdrip[u,s,l]
     >     tblow,       ! e-folding time of snow blowoff tblow[u,s,l]
     >     wliqmax,     ! max amount of intercepted liquid wliq[u,s,l]max
     >     wsnomax,     ! max amount of intercepted snow   wsno[u,s,l]max
     >     wliqmin,     ! min amount of intercepted liquid (same name for u,s,l)
     >     wsnomin      ! min amount of intercepted snow (same name for u,s,l)
c
      real       
     >  tveg(npoi),     ! temperature of veg component t[u,s,l]
     >  wliq(npoi),     ! intercepted liquid amount wliq[u,s,l] (returned)
     >  wsno(npoi),     ! intercepted snow amount wsno[u,s,l] (returned)
     >  xai(npoi),      ! lai, sai, lai+sai for upper leaves/stems,lower veg
     >  pflux(npoi),    ! ht flux due to adjust of intercep precip (returned)
     >  rain(npoi),     ! rainfall rate. Input: above veg, Output: below veg
     >  train(npoi),    ! temperature of rain. (returned)
     >  snow(npoi),     ! snowfall rate. Input: above veg, output: below veg
     >  tsnow(npoi)     ! temperature of snow (returned)
c
c local variables:
c
      integer i         ! loop indice
c
      real rwork,       ! 1/dtime
     >     x,           ! work variable
     >     rwork2,      ! work variable: ch2o - cice
     >     dw           ! correction: freezing liguid or melting snow
c
      real fint(npoi),  ! precip fraction intercepted by unit leaf/stem area
     >     drip(npoi),  ! rate of liquid drip
     >     blow(npoi)   ! rate of snow blowoff
c
c ---------------------------------------------------------------------
c
c calculate fint, the intercepted precip fraction per unit
c leaf/stem area -- note 0.5 * lai or sai (similar to irrad)
c 
      do 50 i = 1, npoi
c
        if (xai(i).ge.epsilon) then
          fint(i) = ( 1.-exp(-0.5*xai(i)) )/ xai(i)
        else
          fint(i) = 0.5
        endif
c
   50 continue
c
c step intercepted liquid and snow amounts due to drip/blow,
c intercepted rainfall/snowfall, and min/max limits. also 
c adjust temperature of intercepted precip to current veg temp,
c storing the heat needed to do this in pflux for use in turvap
c 
c without these pfluxes, the implicit turvap calcs could not
c account for the heat flux associated with precip adjustments,
c especially changes of phase (see below), and so could not
c handle equilibrium situations such as intercepted snowfall
c being continuously melted by warm atmos fluxes, with the veg 
c temp somewhat lower than the equil atmos temp to supply heat
c that melts the incoming snow; (turvap would just change veg 
c temp to atmos equil, with little sensible heat storage...then
c final phase adjustment would return veg temp to melt point)
c
c the use of the current (ie, previous timestep's) veg temp 
c gives the best estimate of what this timestep's final temp
c will be, at least for steady conditions
c
      rwork = 1. / dtime
c
      do 100 i=1,npoi
c    
c liquid
c
        drip(i) = xai(i)*wliq(i)/tdrip
        wliq(i) = wliq(i) * (1.-dtime/tdrip)
c
        wliq(i) = wliq(i) + dtime*rain(i)*fint(i)
        pflux(i) = rain(i)*fint(i) * (tveg(i)-train(i))*ch2o
        rain(i) = rain(i)*(1.-xai(i)*fint(i))
c
        x = wliq(i)
        wliq(i) = min (wliq(i), wliqmax)
        if (wliq(i).lt.wliqmin) wliq(i) = 0.
        drip(i) = drip(i) + xai(i)*(x-wliq(i))*rwork
c
c snow
c
        blow(i) = xai(i)*wsno(i)/tblow
        wsno(i) = wsno(i) * (1.-dtime/tblow)
c
        wsno(i) = wsno(i) + dtime*snow(i)*fint(i)
        pflux(i) = pflux(i) + snow(i)*fint(i) * (tveg(i)-tsnow(i))*cice
        snow(i) = snow(i)*(1.-xai(i)*fint(i))
c
        x = wsno(i)
        wsno(i) = min (wsno(i), wsnomax)
        if (wsno(i).lt.wsnomin) wsno(i) = 0. 
        blow(i) = blow(i) + xai(i)*(x-wsno(i))*rwork
c
  100 continue
c
c change phase of liquid/snow below/above melt point, and add
c required heat to pflux (see comments above). this will only
c affect the precip intercepted in this timestep, since original
c wliq, wsno must have been ge/le melt point (ensured in later
c call to cascad2/steph2o2)
c
      rwork2 = ch2o - cice
c
      do 300 i=1,npoi
c
c liquid below freezing
c
        dw = 0.
        if (tveg(i).lt.tmelt)  dw = wliq(i)
c
        pflux(i) = pflux(i)
     >           + dw * (rwork2*(tmelt-tveg(i)) - hfus) * rwork
        wliq(i) = wliq(i) - dw
        wsno(i) = wsno(i) + dw
c
c snow above freezing
c
        dw = 0.
        if (tveg(i).gt.tmelt)  dw = wsno(i)
c
        pflux(i) = pflux(i)
     >           + dw * (rwork2*(tveg(i)-tmelt) + hfus) * rwork
        wsno(i) = wsno(i) - dw
        wliq(i) = wliq(i) + dw
c
  300 continue
c
c adjust rainfall, snowfall below veg for interception 
c and drip, blowoff
c
      call mix (rain,train, rain,train, drip,tveg, vzero,vzero)
      call mix (snow,tsnow, snow,tsnow, blow,tveg, vzero,vzero)
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine cascad2
c ---------------------------------------------------------------------
c
c at end of timestep, removes evaporation from intercepted h2o,
c and does final heat-conserving adjustment for any liquid/snow 
c below/above melt point. calls steph2o2 for upper leaves, 
c upper stems and lower veg in turn.
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'comsno.h'
      include 'comveg.h'
      include 'com1d.h'
c
c local variables
c
      integer i           ! loop indice
c
      real fveg(npoi),    ! fractional areal coverage of veg component
     >     xai(npoi)      ! lai and/or sai for veg component
c
c ---------------------------------------------------------------------
c
c set up for upper leaves
c
      do 100 i=1,npoi
        fveg(i) = fu(i)
        xai(i) = 2.0 * lai(i,2)
  100 continue
c
c step upper leaves
c
      call steph2o2 (tu,wliqu,wsnou,fveg,xai,rliqu,fvapuw,chu)
c
c set up for upper stems
c
      do 200 i=1,npoi
        fveg(i) = fu(i)
        xai(i) = 2.0 * sai(i,2)
  200 continue
c
c step upper stems
c
      call steph2o2 (ts,wliqs,wsnos,fveg,xai,rliqs,fvaps,chs)
c
c set up for lower veg
c
      do 400 i=1,npoi
        fveg(i) = (1.-fi(i))*fl(i)
        xai(i) = 2.0 * (lai(i,1) + sai(i,1))
  400 continue
c
c step lower veg
c
      call steph2o2 (tl,wliql,wsnol,fveg,xai,rliql,fvaplw,chl)
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine steph2o2 (tveg,wliq,wsno,fveg,xai,rliq,fvapw,cveg)
c ---------------------------------------------------------------------
c
c removes evaporation from intercepted h2o, and does final
c heat-conserving adjustment for any liquid/snow below/above
c melt point, for one veg component
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'comatm.h'
      include 'com1d.h'
c
c Arguments (all arguments are supplied unless otherwise noted)
c
      real cveg        ! specific heat of veg component ch[u,s,l] 
c
      real      
     >  tveg(npoi),    ! temperature of veg component t[u,s,l] (returned)
     >  wliq(npoi),    ! intercepted liquid amount wliq[u,s,l] (returned)
     >  wsno(npoi),    ! intercepted snow amount wsno[u,s,l] (returned)
     >  fveg(npoi),    ! fractional areal coverage, fu or (1-fi)*fl
     >  xai(npoi),     ! lai, sai, lai+sai for upper leaves/stems,lower veg
     >  rliq(npoi),    ! ratio of area wetted by liquid to total wetted area
     >  fvapw(npoi)    ! wetted evap h2o flx per leaf/stem area fvap[uw,s,lw]
c
c local variables
c
      integer i        ! loopi indice
c
      real zm,         ! to compute corrective fluxes
     >     rwork,      ! 1/specific heat of fusion 
     >     chav        ! average specific heat for veg, liw and snow
c
      real dh(npoi),   ! correct heat flux for liquid below melt point and opposite
     >     dw(npoi)    ! correct water flux for liquid below melt point and opposite
c
c
      include 'comsat.h'
c
c ---------------------------------------------------------------------
c
c step intercepted h2o due to evaporation/sublimation.
c (fvapw already has been multiplied by fwet factor in turvap,
c so it is per unit leaf/stem area.)
c
c due to linear fwet factors (see comments in fwetcal) and
c the cap on suw,ssw,slw in turvap, evaporation in one timestep
c should hardly ever make wliq or wsno negative -- but if this
c happens, compensate by increasing vapor flux from atmosphere, 
c and decreasing sensib heat flux from atmos (the former is
c dangerous since it could suck moisture out of a dry atmos,
c and both are unphysical but do fix the budget) tveg in hvapf
c and hsubf should be pre-turvap-timestep values, but are not
c
      do 100 i = 1, npoi
c
        wliq(i) = wliq(i) - dtime *     rliq(i)  * fvapw(i)
        wsno(i) = wsno(i) - dtime * (1.-rliq(i)) * fvapw(i)
c
c check to see if predicted wliq or wsno are less than zero
c
        if ((wliq(i).lt.0. or. wsno(i).lt.0.)
     >      .and. fveg(i)*xai(i).gt.0. )  then
c
c         write (*,9999) i, wliq(i), wsno(i)
c9999     format(' ***warning: wliq<0 or wsno<0 -- steph2o2 9999',
c    >           ' i, wliq, wsno:',i4, 2f12.6)
c
c calculate corrective fluxes
c
          zm = max (-wliq(i), 0.) * fveg(i) * xai(i) / dtime
          fvapa(i) = fvapa(i) + zm
          fsena(i) = fsena(i) - zm*hvapf(tveg(i),ta(i))
          wliq(i) = max (wliq(i), 0.)
c
          zm = max (-wsno(i), 0.) * fveg(i) * xai(i) / dtime
          fvapa(i) = fvapa(i) + zm
          fsena(i) = fsena(i) - zm*hsubf(tveg(i),ta(i))
          wsno(i) = max (wsno(i), 0.)
c
        endif
c
  100 continue
c
c final heat-conserving correction for liquid/snow below/above
c melting point
c
      rwork = 1. / hfus
c
      do 200 i=1,npoi
c
        chav = cveg + ch2o*wliq(i) + cice*wsno(i)
c
c correct for liquid below melt point
c
c (nb: if tveg > tmelt or wliq = 0, nothing changes.)
c
        if (tveg(i).lt.tmelt .and. wliq(i).gt.0.0) then
          dh(i) = chav*(tmelt - tveg(i))
          dw(i) = min (wliq(i), max (0., dh(i)*rwork))
          wliq(i) = wliq(i) - dw(i)
          wsno(i) = wsno(i) + dw(i) 
          chav = cveg + ch2o*wliq(i) + cice*wsno(i)
          tveg(i) = tmelt - (dh(i)-hfus*dw(i))/chav
        endif
c
c correct for snow above melt point
c
c (nb: if tveg < tmelt or wsno = 0, nothing changes.)
c
        if (tveg(i).gt.tmelt .and. wsno(i).gt.0.0) then
          dh(i) = chav*(tveg(i) - tmelt)
          dw(i) = min (wsno(i), max (0., dh(i)*rwork))
          wsno(i) = wsno(i) - dw(i)
          wliq(i) = wliq(i) + dw(i)
          chav = cveg + ch2o*wliq(i) + cice*wsno(i)
          tveg(i) = tmelt + (dh(i)-hfus*dw(i))/chav
        endif
c
  200 continue
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine mix (xm,tm, x1,t1, x2,t2, x3,t3)
c ---------------------------------------------------------------------
c
c calorimetrically mixes masses x1,x2,x3 with temperatures
c t1,t2,t3 into combined mass xm with temperature tm
c
c xm,tm may be returned into same location as one of x1,t1,..,
c so hold result temporarily in xtmp,ttmp below
c
c will work if some of x1,x2,x3 have opposite signs, but may 
c give unphysical tm's
c
      include 'implicit.h'
c
      include 'compar.h'
c
c Arguments (input except for xm, tm)
c
      real xm(npoi),     ! resulting mass  
     >     tm(npoi),     ! resulting temp
     >     x1(npoi),     ! mass 1
     >     t1(npoi),     ! temp 1
     >     x2(npoi),     ! mass 2
     >     t2(npoi),     ! temp 2
     >     x3(npoi),     ! mass 3
     >     t3(npoi)      ! temp 3
c
c local variables
c
      integer i          ! loop indice
c
      real xtmp,         ! resulting mass (storing variable)
     >     ytmp,         !  "
     >     ttmp          ! resulting temp
c
c ---------------------------------------------------------------------
c
      do 100 i=1,npoi
c
        xtmp = x1(i) + x2(i) + x3(i)
c
        ytmp = sign (max (abs(xtmp), epsilon), xtmp)
c
        if (abs(xtmp).ge.epsilon) then
          ttmp = (t1(i)*x1(i) + t2(i)*x2(i) + t3(i)*x3(i)) / ytmp
        else
          ttmp = 0.
          xtmp = 0.
        endif
c
        xm(i) = xtmp
        tm(i) = ttmp
c
  100 continue
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine noveg
c ---------------------------------------------------------------------
c
c if no veg surfaces exist, set prog vars to nominal values
c
c (sensible fluxes fsen[u,s,l], latent fluxes fvap[u,s,l]*, 
c temperature t[u,s,l], and intercepted liquid, snow amounts 
c wliq[u,s,l], wsno[u,s,l] have been calculated for a unit 
c leaf/stem surface, whether or not one exists.)
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'comsno.h'
      include 'comsoi.h'
      include 'comveg.h'
c
c local variables
c
      integer i   ! loop indice
c
      real tav,   ! average temp for soil and snow 
     >     x,     ! total lai + sai
     >     y      ! fraction of lower canopy not snow covered 
c
      do 100 i = 1, npoi
c
        tav = (1.-fi(i))*tg(i) + fi(i)*ti(i)
c
        if (lai(i,2).eq.0. .or. fu(i).eq.0.) then
          tu(i) = tav
          wliqu(i) = 0.
          wsnou(i) = 0.
        endif
c
        if (sai(i,2).eq.0. .or. fu(i).eq.0.) then
          ts(i) = tav
          wliqs(i) = 0.
          wsnos(i) = 0.
        endif 
c
        x = 2.0 * (lai(i,1) + sai(i,1))
        y = fl(i)*(1.-fi(i))
c
        if (x .eq.0. .or. y.eq.0.) then
          tl(i) = tav 
          wliql(i) = 0.
          wsnol(i) = 0.
        endif
c
  100 continue
c
      return
      end
c
c ------------------------------------------------------------------------
      real function twet3(tak, q, p)
c ------------------------------------------------------------------------
c
c twet3.f last update 8/30/2000 C Molling
c
c This function calculates the wet bulb temperature given
c air temp, specific humidity, and air pressure.  It needs the function esat
c in order to work (in comsat.h).  The function is an approximation to
c the actual wet bulb temperature relationship.  It agrees well with the
c formula in the Smithsonian Met. Tables for moderate humidities, but differs
c by as much as 1 K in extremely dry or moist environments.
c
c INPUT
c     tak - air temp in K
c     q - specific humidity in kg/kg
c     p - air pressure in Pa (Pa = 100 * mb)
c
c OUTPUT
c     twet3 - wet bulb temp in K, accuracy?
c
      include 'implicit.h'
      include 'compar.h'
c
      integer i
c
      real tak, q, p, ta, twk, twold, diff
c
      include 'comsat.h'
c
c temperatures in twet3 equation must be in C
c pressure in qsat function must be in Pa
c temperatures in esat,hvapf functions must be in K
c
c     Air temp in C
c     -------------
      ta = tak - 273.16
c
c     First guess for wet bulb temp in C, K
c     -------------------------------------
      twet3 = ta * q / qsat(esat(tak),p)
      twk = twet3 + 273.16
c
c     Iterate to converge
c     -------------------
      do 100 i = 1, 20
         twold = twk - 273.16
         twet3 = ta - (hvapf(twk,tak)/cair) * ( qsat( esat(twk),p )-q )
         diff = twet3 - twold
c
c below, the 0.2 is the relaxation parameter that works up to 40C (at least)
c
         twk = twold + 0.2 * diff + 273.16
         if (abs(twk-273.16-twold) .lt. 0.02) goto 999
 100  continue
c
      print *, 'Warning, twet3 failed to converge after 20 iterations!'
      print *, 'twet3, twetold: ', twk, twold+273.16
      print *, 'twetbulb is being set to the air temperature'
c
      twet3 = tak
c
c     Return wet bulb temperature in K
c     --------------------------------
 999  twet3 = twk
c
      return
      end
c
