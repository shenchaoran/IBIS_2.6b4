c
c #####     ##    #####      #      ##     #####     #     ####   #    #
c #    #   #  #   #    #     #     #  #      #       #    #    #  ##   #
c #    #  #    #  #    #     #    #    #     #       #    #    #  # #  #
c #####   ######  #    #     #    ######     #       #    #    #  #  # #
c #   #   #    #  #    #     #    #    #     #       #    #    #  #   ##
c #    #  #    #  #####      #    #    #     #       #     ####   #    #
c
c ---------------------------------------------------------------------
      subroutine solset
c ---------------------------------------------------------------------
c
c zeros albedos and internal absorbed solar fluxes, and sets
c index for other solar routines. the index indsol, with number
c of points nsol, points to current 1d strip arrays whose coszen 
c values are gt 0 (indsol, nsol are in com1d)
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'comatm.h'
      include 'comsno.h'
      include 'comveg.h'
      include 'com1d.h'
c
      integer i
c
c zero albedos returned just as a niceity
c
      call const (asurd, npoi*nband, 0.0)
      call const (asuri, npoi*nband, 0.0)
c
c zeros absorbed solar fluxes sol[u,s,l,g,i]1 since only points
c with +ve coszen will be set in solarf, and since
c sol[u,l,s,g,i]1 are summed over wavebands in solarf
c
c similarly zero par-related arrays set in solarf for turvap
c
      call const (solu, npoi, 0.0)
      call const (sols, npoi, 0.0)
      call const (soll, npoi, 0.0)
      call const (solg, npoi, 0.0)
      call const (soli, npoi, 0.0)
c
      call const (topparu, npoi, 0.0)
      call const (topparl, npoi, 0.0)
c
c set canopy scaling coefficients for night-time conditions
c
      call const (scalcoefl, npoi*4, 0.0)
      call const (scalcoefu, npoi*4, 0.0)
c
c set index of points with positive coszen
c
      nsol = 0
c
      do 300 i = 1, npoi
        if (coszen(i).gt.0.) then
          nsol = nsol + 1
          indsol(nsol) = i
        endif
  300 continue
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine solsur (ib)
c ---------------------------------------------------------------------
c
c sets surface albedos for soil and snow, prior to other
c solar calculations
c
c ib = waveband number
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'comatm.h'
      include 'comsno.h'
      include 'comsoi.h'
      include 'com1d.h'
c
c input variable
c
      integer ib    ! waveband number. 1 = visible, 2 = near IR
c
c local variables
c     
      integer j,    ! loop indice on number of points with >0 coszen
     >        i     ! indice of point in (1, npoi) array. 
c
      real a7svlo,  ! snow albedo at low threshold temp., visible
     >     a7snlo,  !                                   , near IR
     >     a7svhi,  !                 high              , visible
     >     a7snhi,  !                                   , near-IR
     >     t7shi,   ! high threshold temperature for snow albed
     >     t7slo,   ! low  threshold temperature for snow albedo
     >     dinc,    ! albedo correction du to soil moisture
     >     zw       ! liquid moisture content

      real x(npoi), zfac(npoi)
c
c set the "standard" snow values:
c
      data a7svlo, a7svhi /0.90, 0.70/
      data a7snlo, a7snhi /0.60, 0.40/
c
c     t7shi ... high threshold temperature for snow albedo
c     t7slo ... low  threshold temperature for snow albedo
c
      t7shi = tmelt
      t7slo = tmelt - 15.0
c
c do nothing if all points in current strip have coszen le 0
c
      if (nsol.eq.0) then
        return
      endif
c
      if (ib.eq.1) then
c
c soil albedos (visible waveband)
c
        do 100 j = 1, nsol
c
          i = indsol(j)
c
c change the soil albedo as a function of soil moisture
c
          zw = wsoi(i,1) * (1.-wisoi(i,1))
c
          dinc = 1.0 + 1.0 * min (1., max (0.0, 1. - (zw /.50) ))
c
          albsod(i) = min (albsav(i) * dinc, .80)
          albsoi(i) = albsod(i)
c
  100   continue
c
c snow albedos (visible waveband)
c
        do 110 j = 1, nsol
c
          i = indsol(j)
c
          x(i) = (a7svhi*(tsno(i,1)-t7slo) + a7svlo*(t7shi-tsno(i,1)))
     >           / (t7shi-t7slo)
c
          x(i) = min (a7svlo, max (a7svhi, x(i)))
c
          zfac(i)   = max ( 0., 1.5 / (1.0 + 4.*coszen(i)) - 0.5 )
          albsnd(i) = min (0.99, x(i) + (1.-x(i))*zfac(i))
          albsni(i) = min (1., x(i))
c
  110   continue
c
      else
c
c soil albedos (near-ir waveband)
c
        do 200 j = 1, nsol
          i = indsol(j)
c
c lsx.2 formulation (different from lsx.1)
c
          zw = wsoi(i,1) * (1. - wisoi(i,1))
c
          dinc = 1.0 + 1.0 * min (1., max (0.0, 1.0 - (zw / .50)  ))
c
          albsod(i) = min (albsan(i) * dinc, .80)
          albsoi(i) = albsod(i)
c
  200   continue
c
c snow albedos (near-ir waveband)
c
        do 210 j = 1, nsol
c
          i = indsol(j)
c
          x(i) = (a7snhi*(tsno(i,1)-t7slo) + a7snlo*(t7shi-tsno(i,1)))
     >           / (t7shi-t7slo)
          x(i) = min (a7snlo, max (a7snhi, x(i)))
c
          zfac(i) = max ( 0., 1.5/(1.+4.*coszen(i)) - 0.5 )
c
          albsnd(i) = min (0.99, x(i) + (1.-x(i))*zfac(i))
          albsni(i) = min (1., x(i))
c
  210   continue
c
      endif
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine solalb (ib)
c ---------------------------------------------------------------------
c
c calculates effective albedos of the surface system,
c separately for unit incoming direct and diffuse flux -- the 
c incoming direct zenith angles are supplied in comatm array 
c coszen, and the effective albedos are returned in comatm
c arrays asurd, asuri -- also detailed absorbed and reflected flux
c info is stored in com1d arrays, for later use by solarf
c
c the procedure is first to calculate the grass+soil albedos,
c then the tree + (grass+soil+snow) albedos. the labels
c (a) to (d) correspond to those in the description doc
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
c Arguments
c 
      integer ib     ! waveband number (1= visible, 2= near-IR)
c
c local variables
c     
      integer j,    ! loop indice on number of points with >0 coszen
     >        i     ! indice of point in (1, npoi) array. 
c
c do nothing if all points in current strip have coszen le 0
c
      if (nsol.eq.0) return
c
c (a) obtain albedos, etc, for two-stream lower veg + soil
c     system, for direct and diffuse incoming unit flux
c
      do 100 j = 1, nsol
c
        i = indsol(j)
c
        asurd(i,ib) = albsod(i)
        asuri(i,ib) = albsoi(i)
c
  100 continue
c
      call twostr (ablod, abloi,  relod, reloi,  flodd,  dummy,
     >             flodi, floii,  asurd,  asuri,    1,   coszen, ib)
c
c (b) areally average surface albedos (lower veg, soil, snow)
c
      do 200 j = 1, nsol
c
        i = indsol(j)
c
        asurd(i,ib) = fl(i)*(1.-fi(i))*relod(i)
     >              + (1.-fl(i))*(1.-fi(i))*albsod(i)
     >              + fi(i)*albsnd(i)    
c
        asuri(i,ib) = fl(i)*(1.-fi(i))*reloi(i)
     >              + (1.-fl(i))*(1.-fi(i))*albsoi(i)
     >              + fi(i)*albsni(i)    
c
  200 continue
c
c (c) obtain albedos, etc, for two-stream upper veg + surface
c     system, for direct and diffuse incoming unit flux
c
      call twostr (abupd, abupi,  reupd, reupi,  fupdd,  dummy,
     >             fupdi, fupii,  asurd,  asuri,    2,   coszen, ib)
c
c (d) calculate average overall albedos 
c
      do 300 j = 1, nsol
c
        i = indsol(j)
c
        asurd(i,ib) = fu(i)*reupd(i)
     >              + (1.-fu(i))*asurd(i,ib)
c
        asuri(i,ib) = fu(i)*reupi(i)
     >              + (1.-fu(i))*asuri(i,ib)
c
  300 continue
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine solarf (ib)
c ---------------------------------------------------------------------
c
c calculates solar fluxes absorbed by upper and lower stories,
c soil and snow
c
c zenith angles are in comatm array coszen, and must be the same
c as supplied earlier to solalb
c
c solarf uses the results obtained earlier by solalb and 
c stored in com1d arrays. the absorbed fluxes are returned in
c com1d arrays sol[u,s,l,g,i]
c
c the procedure is first to calculate the upper-story absorbed
c fluxes and fluxes below the upper story, then the lower-story
c absorbed fluxes and fluxes below the lower story, then fluxes
c absorbed by the soil and snow
c
c ib = waveband number
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
c Arguments
c 
      integer ib     ! waveband number (1= visible, 2= near-IR)
c
c local variables
c     
      integer j,    ! loop indice on number of points with >0 coszen
     >        i     ! indice of point in (1, npoi) array. 
c
      real x, y, xd, xi, 
     >     xaiu,    ! total single-sided lai+sai, upper
     >     xail     ! total single-sided lai+sai, lower
c
c do nothing if all points in current strip have coszen le 0
c
      if (nsol.eq.0) return
c
c (f) calculate fluxes absorbed by upper leaves and stems,
c     and downward fluxes below upper veg, using unit-flux
c     results of solalb(c) (apportion absorbed flux between
c     leaves and stems in proportion to their lai and sai)
c
      do 600 j=1,nsol
c
        i = indsol(j)
        x = solad(i,ib)*abupd(i) + solai(i,ib)*abupi(i)
        y = lai(i,2) / max (lai(i,2)+sai(i,2), epsilon)
        solu(i) = solu(i) + x * y
        sols(i) = sols(i) + x * (1.-y)
        sol2d(i) = solad(i,ib)*fupdd(i)
        sol2i(i) = solad(i,ib)*fupdi(i) + solai(i,ib)*fupii(i)
c
  600 continue
c
c (g) areally average fluxes to lower veg, soil, snow
c
      do 700 j=1,nsol
c
        i = indsol(j)
        sol3d(i) = fu(i)*sol2d(i) + (1.-fu(i))*solad(i,ib)
        sol3i(i) = fu(i)*sol2i(i) + (1.-fu(i))*solai(i,ib)
c
  700 continue
c
c (h,i) calculate fluxes absorbed by lower veg, snow-free soil
c       and snow, using results of (g) and unit-flux results
c       of solalb(a)
c
      do 800 j=1,nsol
c
        i = indsol(j)
        soll(i) = soll(i) + sol3d(i)*ablod(i) + sol3i(i)*abloi(i)
c
        xd = (fl(i)*flodd(i) + 1.-fl(i)) * sol3d(i)
c
        xi = fl(i)*(sol3d(i)*flodi(i) + sol3i(i)*floii(i))
     >       + (1.-fl(i)) * sol3i(i)
c
        solg(i) = solg(i)
     >            + (1.-albsod(i))*xd + (1.-albsoi(i))*xi
c
        soli(i) = soli(i) 
     >            + (1.-albsnd(i))*sol3d(i)
     >            + (1.-albsni(i))*sol3i(i)
c
  800 continue
c
c estimate absorbed pars at top of canopy, toppar[u,l] and
c some canopy scaling parameters
c
c this neglects complications due to differing values of dead vs 
c live elements, averaged into rhoveg, tauveg in vegdat, and 
c modifications of omega due to intercepted snow in twoset
c
c do only for visible band (ib=1)
c
      if (ib.eq.1) then
c
        do 900 j = 1, nsol
c
          i = indsol(j)
c
c the canopy scaling algorithm assumes that the net photosynthesis
c is proportional to absored par (apar) during the daytime. during night,
c the respiration is scaled using a 10-day running-average daytime canopy
c scaling parameter.
c
c apar(x) = A exp(-k x) + B exp(-h x) + C exp(h x)
c
c some of the required terms (i.e. term[u,l] are calculated in the subroutine 'twostr'.
c in the equations below, 
c
c   A = scalcoefu(i,1) = term[u,l](i,1) * ipardir(0)
c   B = scalcoefu(i,2) = term[u,l](i,2) * ipardir(0) + term[u,l](i,3) * ipardif(0)
c   C = scalcoefu(i,3) = term[u,l](i,4) * ipardir(0) + term[u,l](i,5) * ipardif(0)
c   A + B + C = scalcoefu(i,4) = also absorbed par at canopy of canopy by leaves & stems
c
c upper canopy:
c
c total single-sided lai+sai
c
          xaiu = max (lai(i,2)+sai(i,2), epsilon)
c
c some terms required for use in canopy scaling:
c
          scalcoefu(i,1) = termu(i,1) * solad(i,ib)
c
          scalcoefu(i,2) = termu(i,2) * solad(i,ib) +
     >                     termu(i,3) * solai(i,ib)
c
          scalcoefu(i,3) = termu(i,4) * solad(i,ib) +
     >                     termu(i,5) * solai(i,ib)
c
          scalcoefu(i,4) = scalcoefu(i,1) + 
     >                     scalcoefu(i,2) + 
     >                     scalcoefu(i,3)
c
c apar of the "top" leaves of the canopy
c
          topparu(i) = scalcoefu(i,4) * lai(i,2) / xaiu
c
c lower canopy:
c
c total single-sided lai+sai
c
          xail = max (lai(i,1)+sai(i,1), epsilon)
c
c some terms required for use in canopy scaling:
c
          scalcoefl(i,1) = terml(i,1) * sol3d(i)
c
          scalcoefl(i,2) = terml(i,2) * sol3d(i) +
     >                     terml(i,3) * sol3i(i)
c
          scalcoefl(i,3) = terml(i,4) * sol3d(i) +
     >                     terml(i,5) * sol3i(i)
c
          scalcoefl(i,4) = scalcoefl(i,1) +
     >                     scalcoefl(i,2) +
     >                     scalcoefl(i,3)
c
c apar of the "top" leaves of the canopy
c
          topparl(i) = scalcoefl(i,4) * lai(i,1) / xail
c
  900   continue
c
      endif
c
      return
      end
c
c
c ------------------------------------------------------------------------
      subroutine twostr (abvegd, abvegi, refld, refli, fbeldd, fbeldi,
     >                   fbelid, fbelii, asurd, asuri, iv, coszen, ib)
c ------------------------------------------------------------------------
c
c solves canonical radiative transfer problem of two-stream veg
c layer + underlying surface of known albedo, for unit incoming
c direct or diffuse flux. returns flux absorbed within layer,
c reflected flux, and downward fluxes below layer. note that all
c direct fluxes are per unit horizontal zrea, ie, already 
c including a factor cos (zenith angle)
c
c the solutions for the twostream approximation follow Sellers (1985),
c and Bonan (1996) (the latter being the LSM documentation)
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'comveg.h'
      include 'com1d.h'
c
c Arguments
c
      integer ib,             ! waveband number (1= visible, 2= near-IR)
     >        iv              ! 1 for lower, 2 for upper story params (supplied)
c
      real abvegd(npoi),      ! direct flux absorbed by two-stream layer (returned)
     >     abvegi(npoi),      ! diffuse flux absorbed by two-stream layer (returned)
     >     refld(npoi),       ! direct flux reflected above two-stream layer (returned)
     >     refli(npoi),       ! diffuse flux reflected above two-stream layer (returned)
     >     fbeldd(npoi),      ! downward direct  flux below two-stream layer(returned)
     >     fbeldi(npoi),      ! downward direct  flux below two-stream layer(returned)
     >     fbelid(npoi),      ! downward diffuse flux below two-stream layer(returned)
     >     fbelii(npoi),      ! downward diffuse flux below two-stream layer(returned)
     >     asurd(npoi,nband), ! direct  albedo of underlying surface (supplied)
     >     asuri(npoi,nband), ! diffuse albedo of underlying surface (supplied)
     >     coszen(npoi)       ! cosine of direct zenith angle (supplied, must be gt 0)
c
c local variables
c
      integer j,    ! loop indice on number of points with >0 coszen
     >        i     ! indice of point in (1, npoi) array. 
c
      real b, c, c0, d, f, h, k, q, p, sigma
c
      real ud1, ui1, ud2, ui2, ud3, xai, s1, s2, p1, p2, p3, p4, 
     >     rwork, dd1, di1, dd2, di2, h1, h2, h3, h4, h5, h6, h7, h8, 
     >     h9, h10, absurd, absuri
c
c [d,i] => per unit incoming direct, diffuse (indirect) flux
c
      real omega(npoi),       !
     >     betad(npoi),       !
     >     betai(npoi),       !
     >     avmu(npoi),        !
     >     gdir(npoi),        !
     >     tmp0(npoi)         !
c
c do nothing if all points in current strip have coszen le 0
c
      if (nsol.eq.0) return
c
c calculate two-stream parameters omega, betad, betai, avmu, gdir
c
      call twoset (omega, betad, betai, avmu, gdir, coszen, iv, ib)
c
      do 100 j=1,nsol
c
        i = indsol(j)
c
c the notations used here are taken from page 21 of Bonan's LSM documentation:
c Bonan, 1996: A Land Surface Model (LSM version 1.0) for ecological, hydrological,
c and atmospheric studies: Technical description and user's guide. NCAR Technical
c Note. NCAR/TN-417+STR, January 1996.
c
c some temporary variables are also introduced, which are from the original
c lsx model.
c
        b = 1. - omega(i) * (1.-betai(i))
        c = omega(i) * betai(i)
c
        tmp0(i) = b*b-c*c
c
        q = sqrt ( max(0.0, tmp0(i)) )
        k = gdir(i) / max(coszen(i), 0.01)
        p = avmu(i) * k
c
c next line perturbs p if p = q
c
        if ( abs(p-q) .lt. .001*p )
     >  p = (1.+sign(.001,p-q)) * p
c
        c0 = omega(i) * p
        d = c0 * betad(i)
        f = c0 * (1.-betad(i))
        h = q / avmu(i)
c
        sigma = p*p - tmp0(i)
c
c direct & diffuse parameters are separately calculated
c
        ud1 = b - c/asurd(i,ib)
        ui1 = b - c/asuri(i,ib)
        ud2 = b - c*asurd(i,ib)
        ui2 = b - c*asuri(i,ib)
        ud3 = f + c*asurd(i,ib)
c
        xai = max (lai(i,iv) + sai(i,iv), epsilon)
c
        s1 = exp(-1.*h*xai)
        s2 = exp(-1.*k*xai)
c
        p1 = b + q
        p2 = b - q
        p3 = b + p
        p4 = b - p
        rwork = 1./s1
c
c direct & diffuse parameters are separately calculated
c
        dd1 = p1*(ud1-q)*rwork - p2*(ud1+q)*s1
        di1 = p1*(ui1-q)*rwork - p2*(ui1+q)*s1
        dd2 = (ud2+q)*rwork - (ud2-q)*s1
        di2 = (ui2+q)*rwork - (ui2-q)*s1
        h1 = -1.*d*p4 - c*f
        rwork = s2*(d-c-h1*(ud1+p)/sigma)
        h2 = 1./dd1*( (d-h1*p3/sigma)*(ud1-q)/s1 - 
     >       p2*rwork )
        h3 = -1./dd1*( (d-h1*p3/sigma)*(ud1+q)*s1 - 
     >       p1*rwork )
        h4 = -1.*f*p3 - c*d
        rwork = s2*(ud3-h4*(ud2-p)/sigma)
        h5 = -1./dd2*( h4*(ud2+q)/(sigma*s1) +
     >       rwork )
        h6 = 1./dd2*( h4*s1*(ud2-q)/sigma +
     >       rwork )
        h7 = c*(ui1-q)/(di1*s1)
        h8 = -1.*c*s1*(ui1+q)/di1
        h9 = (ui2+q)/(di2*s1)
        h10= -1.*s1*(ui2-q)/di2
c
c save downward direct, diffuse fluxes below two-stream layer
c
        fbeldd(i) = s2
        fbeldi(i) = 0.
        fbelid(i) = h4/sigma*s2 + h5*s1 + h6/s1
        fbelii(i) = h9*s1 + h10/s1
c
c save reflected flux, and flux absorbed by two-stream layer
c
        refld(i) = h1/sigma + h2 + h3
        refli(i) = h7 + h8
        absurd = (1.-asurd(i,ib)) * fbeldd(i)
     >         + (1.-asuri(i,ib)) * fbelid(i)
        absuri = (1.-asuri(i,ib)) * fbelii(i)
c
        abvegd(i) = max (0., 1. - refld(i) - absurd)
        abvegi(i) = max (0., 1. - refli(i) - absuri)
c
c if no veg, make sure abveg (flux absorbed by veg) is exactly zero
c if this is not done, roundoff error causes small (+/-)
c sols, soll values in solarf and subsequent problems in turvap
c via stomata
c
        if (xai.lt.epsilon) abvegd(i) = 0.0
        if (xai.lt.epsilon) abvegi(i) = 0.0
c
c some terms needed in canopy scaling
c the canopy scaling algorithm assumes that the net photosynthesis
c is proportional to absored par (apar) during the daytime. during night,
c the respiration is scaled using a 10-day running-average daytime canopy
c scaling parameter.
c
c apar(x) = A exp(-k x) + B exp(-h x) + C exp(h x)
c
c in the equations below, 
c
c   k = term[u,l](i,6)
c   h = term[u,l](i,7)
c
c   A = term[u,l](i,1) * ipardir(0)
c   B = term[u,l](i,2) * ipardir(0) + term[u,l](i,3) * ipardif(0)
c   C = term[u,l](i,4) * ipardir(0) + term[u,l](i,5) * ipardif(0)
c
c calculations performed only for visible (ib=1)
c
      if (ib.eq.1) then
c
        if (iv.eq.1) then
          terml(i,1) = k * (1. + (h4-h1) / sigma)
          terml(i,2) = h * (h5 - h2)
          terml(i,3) = h * (h9 - h7)
          terml(i,4) = h * (h3 - h6)
          terml(i,5) = h * (h8 - h10)
          terml(i,6) = k
          terml(i,7) = h
        else
          termu(i,1) = k * (1. + (h4-h1) / sigma)
          termu(i,2) = h * (h5 - h2)
          termu(i,3) = h * (h9 - h7)
          termu(i,4) = h * (h3 - h6)
          termu(i,5) = h * (h8 - h10)
          termu(i,6) = k
          termu(i,7) = h
        endif
c
      end if
c
  100 continue
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine twoset (omega, betad, betai, avmu, gdir,
     >                   coszen, iv, ib)  
c ---------------------------------------------------------------------
c
c sets two-stream parameters, given single-element transmittance
c and reflectance, leaf orientation weights, and cosine of the
c zenith angle, then adjusts for amounts of intercepted snow
c
c the two-stream parameters omega,betad,betai are weighted 
c combinations of the "exact" values for the 3 orientations:
c all vertical, all horizontal, or all random (ie, spherical)
c
c the vertical, horizontal weights are in oriev,orieh (comveg)
c
c the "exact" expressions are as derived in my notes(8/6/91,p.6).
c note that values for omega*betad and omega*betai are calculated
c and then divided by the new omega, since those products are 
c actually used in twostr. also those depend *linearly* on the
c single-element transmittances and reflectances tauveg, rhoveg,
c which are themselves linear weights of leaf and stem values 
c
c for random orientation, omega*betad depends on coszen according
c to the function in array tablemu
c
c the procedure is approximate since omega*beta[d,i] and gdir
c should depend non-linearly on the complete leaf-angle
c distribution. then we should also treat leaf and stem angle
c distributions separately, and allow for the cylindrical
c shape of stems (norman and jarvis, app.b; the expressions 
c below are appropriate for flat leaves)
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'comveg.h'
      include 'com1d.h'
c
c Arguments (all quantities are returned unless otherwise note)
c
      integer ib,             ! waveband number (1= visible, 2= near-IR)
     >        iv              ! 1 for lower, 2 for upper story params (supplied)
c
      real omega(npoi),       ! fraction of intercepted radiation that is scattered
     >     betad(npoi),       ! fraction of scattered *direct* radiation that is
     >                        !  scattered into upwards hemisphere
     >     betai(npoi),       ! fraction of scattered downward *diffuse* radiation
     >                        ! that is scattered into upwards hemisphere (or fraction
     >                        ! of scattered upward diffuse rad. into downwards hemis)
     >     avmu(npoi),        ! average diffuse optical depth
     >     gdir(npoi),        ! average projected leaf area into solar direction
     >     coszen(npoi)       ! cosine of solar zenith angle (supplied)
c
c local variables
c
      integer j,    ! loop indice on number of points with >0 coszen
     >        i,    ! indice of point in (1, npoi) array. 
     >        ntmu, !
     >        itab
c
      real zrho, ztau, orand, ztab, rwork, y, o, x, betadsno, betaisno
c
      real otmp(npoi)
c
      parameter (ntmu=100)
      real tablemu(ntmu+1), omegasno(nband)
      save tablemu, omegasno, betadsno, betaisno
c
      data tablemu /
     >   0.5000, 0.4967, 0.4933, 0.4900, 0.4867, 0.4833, 0.4800, 0.4767,
     >   0.4733, 0.4700, 0.4667, 0.4633, 0.4600, 0.4567, 0.4533, 0.4500,
     >   0.4467, 0.4433, 0.4400, 0.4367, 0.4333, 0.4300, 0.4267, 0.4233,
     >   0.4200, 0.4167, 0.4133, 0.4100, 0.4067, 0.4033, 0.4000, 0.3967,
     >   0.3933, 0.3900, 0.3867, 0.3833, 0.3800, 0.3767, 0.3733, 0.3700,
     >   0.3667, 0.3633, 0.3600, 0.3567, 0.3533, 0.3500, 0.3467, 0.3433,
     >   0.3400, 0.3367, 0.3333, 0.3300, 0.3267, 0.3233, 0.3200, 0.3167,
     >   0.3133, 0.3100, 0.3067, 0.3033, 0.3000, 0.2967, 0.2933, 0.2900,
     >   0.2867, 0.2833, 0.2800, 0.2767, 0.2733, 0.2700, 0.2667, 0.2633,
     >   0.2600, 0.2567, 0.2533, 0.2500, 0.2467, 0.2433, 0.2400, 0.2367,
     >   0.2333, 0.2300, 0.2267, 0.2233, 0.2200, 0.2167, 0.2133, 0.2100,
     >   0.2067, 0.2033, 0.2000, 0.1967, 0.1933, 0.1900, 0.1867, 0.1833,
     >   0.1800, 0.1767, 0.1733, 0.1700, 0.1667 /
c
      data omegasno /0.9, 0.7/
      data betadsno, betaisno /0.5, 0.5/
c
c set two-stream parameters omega, betad, betai, gdir and avmu
c as weights of those for 100% vert,horiz,random orientations
c
      do 100 j=1,nsol
        i = indsol(j)
c
        zrho = rhoveg(ib,iv)
        ztau = tauveg(ib,iv)
c
c weight for random orientation is 1 - those for vert and horiz
c
        orand = 1. - oriev(iv) - orieh(iv)
c
        omega(i) = zrho + ztau
c
c ztab is transmittance coeff - for random-orientation omega*betad,
c given by tablemu as a function of coszen
c
        itab = nint (coszen(i)*ntmu + 1)
        ztab = tablemu(itab)
        rwork = 1./omega(i)
c
        betad(i) = (  oriev(iv) * 0.5*(zrho + ztau)
     >              + orieh(iv) * zrho
     >              + orand       * ((1.-ztab)*zrho + ztab*ztau) )
     >             * rwork
c
        betai(i) = (  oriev(iv) * 0.5*(zrho + ztau)
     >              + orieh(iv) * zrho
     >              + orand       * ((2./3.)*zrho + (1./3.)*ztau) )
     >             * rwork
c
        gdir(i)  = oriev(iv) * (2./pi) *
     >             sqrt ( max (0., 1.-coszen(i)*coszen(i)) )
     >           + orieh(iv) * coszen(i)
     >           + orand       * 0.5
c
        avmu(i) = 1.
c
  100 continue
c
c adjust omega, betad and betai for amounts of intercepted snow
c (omegasno decreases to .6 of cold values within 1 deg of tmelt)
c
      if (iv.eq.1) then
c
c lower story
c
        do 210 j=1,nsol
          i = indsol(j)
          y = fwetl(i)*(1.-rliql(i))
          o = omegasno(ib)*(.6 + .4*max(0.,min(1.,(tmelt-tl(i))/1.0)))
          otmp(i)  = omega(i)
          rwork = y * o
          omega(i) =  (1-y)*otmp(i)          + rwork
          betad(i) = ((1-y)*otmp(i)*betad(i) + rwork*betadsno) /
     >               omega(i)  
          betai(i) = ((1-y)*otmp(i)*betai(i) + rwork*betaisno) /
     >               omega(i)  
  210   continue
c
      else
c
c upper story
c
        do 220 j=1,nsol
          i = indsol(j)
          x = lai(i,iv) / max (lai(i,iv)+sai(i,iv), epsilon)
          y = x * fwetu(i)*(1.-rliqu(i)) + (1-x) *fwets(i)*(1.-rliqs(i))
          o = (     x  * min (1., max (.6, (tmelt-tu(i))/0.1))
     >         + (1-x) * min (1., max (.6, (tmelt-ts(i))/0.1)) )
     >      *  omegasno(ib) 
c
          otmp(i)  = omega(i)
          rwork = y * o
          omega(i) =  (1-y)*otmp(i)          + rwork
          betad(i) = ((1-y)*otmp(i)*betad(i) + rwork*betadsno) /
     >               omega(i)
          betai(i) = ((1-y)*otmp(i)*betai(i) + rwork*betaisno) /
     >               omega(i)
c
  220   continue
c
      endif
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine irrad
c ---------------------------------------------------------------------
c
c calculates overall emitted ir flux, and net absorbed minus
c emitted ir fluxes for upper leaves, upper stems, lower story,
c soil and snow. assumes upper leaves, upper stems and lower
c story each form a semi-transparent plane, with the upper-leaf
c plane just above the upper-stem plane. the soil and snow 
c surfaces have emissivities of 0.95.
c
c the incoming flux is supplied in comatm array fira
c
c the emitted ir flux by overall surface system is returned in
c com1d array firb - the ir fluxes absorbed by upper leaves,
c upper stems, lower veg, soil and snow are returned in com1d 
c arrays firu, firs, firl, firg and firi
c 
c other com1d arrays used are:
c
c emu, ems, eml  = emissivities of the vegetation planes
c fup, fdown     = upward and downward fluxes below tree level
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
c Local arrays:
c
      integer i           ! loop indice
c
      real emisoil,       ! soil emissivity
     >     emisnow,       ! snow emissivity
     >     avmuir         ! average diffuse optical depth
c
      real emu(npoi),     ! ir emissivity of upper-leaves veg plane
     >     ems(npoi),     ! ir emissivity of upper-stems veg plane
     >     eml(npoi),     ! ir emissivity of lower-story veg plane
     >     emg(npoi),     ! ir emissivity (gray) of soil surface
     >     emi(npoi),     ! ir emissivity (gray) of snow surface
     >     fdown(npoi),   ! downward ir flux below tree level per overall area
     >     fdowng(npoi),  ! upward   ir flux below tree level per overall area
     >     fup(npoi),     ! downward ir flux below lower-story veg
     >     fupg(npoi),    ! upward   ir flux below lower-story veg
     >     fupgb(npoi),   ! upward   ir flux above bare soil surface
     >     fupi(npoi)     ! upward   ir flux above snow surface
c
c set emissivities of soil and snow
c
      data emisoil, emisnow
     >    /0.95, 0.95/
c
c use uniform value 1.0 for average diffuse optical depth
c (although an array for solar, all values are set to 1 in twoset).
c
      save avmuir
      data avmuir /1./
c
      do 100 i=1,npoi
c
        emu(i) = 1. - exp ( -lai(i,2) / avmuir )
        ems(i) = 1. - exp ( -sai(i,2) / avmuir )
        eml(i) = 1. - exp ( -(lai(i,1)+sai(i,1)) / avmuir )
c
        emg(i) = emisoil
        emi(i) = emisnow
c
        fdown(i) =  (1.-fu(i)) * fira(i)
     >            + fu(i) * ( (1.-emu(i))*(1.-ems(i))*fira(i)
     >                       +    emu(i)* (1.-ems(i))*stef*(tu(i)**4)
     >                       +    ems(i)*stef*(ts(i)**4) )
c
        fdowng(i) = (1.-eml(i))*fdown(i)  + eml(i)*stef*(tl(i)**4)
c
        fupg(i)   = (1.-emg(i))*fdowng(i) + emg(i)*stef*(tg(i)**4)
c
        fupgb(i)  = (1.-emg(i))*fdown(i)  + emg(i)*stef*(tg(i)**4)
c
        fupi(i)   = (1.-emi(i))*fdown(i)  + emi(i)*stef*(ti(i)**4)
c
        fup(i) = (1.-fi(i))*(      fl(i)*(       eml(i) *stef*(tl(i)**4)
     >                                     + (1.-eml(i))*fupg(i) )
     >                        +(1.-fl(i))*fupgb(i)
     >                      )
     >         +     fi(i) * fupi(i)
c
        firb(i) =   (1.-fu(i)) * fup(i)
     >            + fu(i)  * ( (1.-emu(i))*(1.-ems(i))*fup(i)
     >                        +    emu(i)*stef*(tu(i)**4)
     >                        +    ems(i)*(1.-emu(i))*stef*(ts(i)**4) )
c
        firu(i) =   emu(i)*ems(i)*stef*(ts(i)**4)
     >            + emu(i)*(1.-ems(i))*fup(i)
     >            + emu(i)*fira(i)
     >            - 2*emu(i)*stef*(tu(i)**4)
c
        firs(i) =   ems(i)*emu(i)*stef*(tu(i)**4)
     >            + ems(i)*fup(i)
     >            + ems(i)*(1.-emu(i))*fira(i)
     >            - 2*ems(i)*stef*(ts(i)**4)
c
        firl(i) =   eml(i)*fdown(i)
     >            + eml(i)*fupg(i)
     >            - 2*eml(i)*stef*(tl(i)**4)
c
        firg(i) =       fl(i)  * (fdowng(i) - fupg(i))
     >            + (1.-fl(i)) * (fdown(i)  - fupgb(i))
c
        firi(i) =   fdown(i) - fupi(i)
c
  100 continue
c
      return
      end
