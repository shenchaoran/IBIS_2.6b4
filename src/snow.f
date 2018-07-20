c
c  ####   #    #   ####   #    #
c #       ##   #  #    #  #    #
c  ####   # #  #  #    #  #    #
c      #  #  # #  #    #  # ## #
c #    #  #   ##  #    #  ##  ##
c  ####   #    #   ####   #    #
c
c ---------------------------------------------------------------------
      subroutine snow
c ---------------------------------------------------------------------
c
c steps snow model through one timestep
c
      use implicit
c
      use compar
      use comsoi
      use comsno
      use comveg
      use com1d
c
c local variables:
c
      integer i, k,             ! loop indices
     >     npn                  ! index indsno, npcounter for pts with snow

      real rwork, rwork2,       ! working vaiable
     >     finew,               ! storing variable for fi
     >     zhh,                 ! 0.5*hsnomin
     >     zdh                  ! max height of snow above hsnomin (?)     

      integer indsno(npoi)      ! index of points with snow in current 1d strip
c
      real hinit(nsnolay),      ! initial layer thicknesses when snow first forms
     >     hsnoruf(npoi),       ! heigth of snow forced to cover lower canopy (?)
     >     fiold(npoi),         ! old fi at start of this timestep
     >     fhtop(npoi),         ! heat flux into upper snow surface
     >     sflo(npoi,nsnolay+2),! heat flux across snow and buried-lower-veg layer bdries
     >     zmelt(npoi),         ! liquid mass flux increments to soil, at temperature 
     >                          ! tmelt, due to processes occuring during this step
     >     zheat(npoi),         ! heat flux to soil, due to processes occuring this step
     >     dfi(npoi),           ! change in fi
     >     xl(npoi),            ! lower veg density
     >     xh(npoi),            ! temporary arrays
     >     xm(npoi),            ! "
     >     ht(npoi),            ! "
     >     x1(npoi),            ! "
     >     x2(npoi),            ! "
     >     x3(npoi)             ! "
c
      do 10 i = 1, npoi
        hsnoruf(i) =  min (0.70, max (hsnomin+.05, fl(i)*ztop(i,1)))
        xl(i) = fl(i) * 2.0 * (lai(i,1) + sai(i,1))
        x1(i) = tlsub(i)
 10   continue
c
      hinit(1) = hsnotop
c
      do 15 k = 2, nsnolay
        hinit(k) = (hsnomin - hsnotop) / (nsnolay-1)
   15 continue
c
      do 20 i = 1, npoi
        fiold(i) = fi(i)
   20 continue
c
c zero out arrays
c
      call const (sflo, npoi*(nsnolay+2), 0.0)
      call const (zmelt, npoi, 0.0)
      call const (zheat, npoi, 0.0)
c
c set up index indsno, npn for pts with snow - indsno is used
c only by vadapt - elsewhere below, just test on npn > 0
c
      npn = 0                                        
c
      do 30 i = 1, npoi 
        if (fi(i).gt.0.) then
          npn = npn + 1
          indsno(npn) = i
        endif 
   30 continue
c
c set surface heat flux fhtop and increment top layer thickness
c due to rainfall, snowfall and sublimation on existing snow
c
      if (npn.gt.0) then
c
        rwork = dtime / rhos
c
        do 40 i = 1, npoi
c
          fhtop(i) = heati(i) +
     >               rainl(i) * (ch2o * (trainl(i) - tmelt)     + hfus +
     >                           cice * (tmelt     - tsno(i,1)))       +
     >               snowl(i) *  cice * (tsnowl(i) - tsno(i,1))
c
          if (fi(i).gt.0.) hsno(i,1) = hsno(i,1) +
     >                     (rainl(i) + snowl(i) - fvapi(i)) * rwork
c
   40   continue
c
      endif
c
c step temperatures due to heat conduction, including buried
c lower-veg temperature tlsub
c
      if (npn.gt.0) then
c
        call scopy (npoi, tlsub, x1)
        call snowheat (tlsub, fhtop, sflo, xl, chl)
c
      endif
c
c put snowfall from 1-fi snow-free area onto side of existing
c snow, or create new snow if current fi = 0. also reset index.
c (assumes total depth of newly created snow = hsnomin.)
c (fi will not become gt 1 here if one timestep's snowfall
c <= hsnomin, but protect against this anyway.)
c
c if no adjacent snowfall or fi = 1, dfi = 0, so no effect
c
      call const (ht, npoi, 0.0)
      do 190 k=1,nsnolay
        do 192 i=1,npoi
          ht(i) = ht(i) + hsno(i,k)
  192   continue
  190 continue
c
      do 195 i=1,npoi
        if (ht(i).eq.0.) ht(i) = hsnomin
  195 continue
c
      rwork = dtime / rhos
      do 200 i=1,npoi
        dfi(i) = (1.-fi(i))*rwork*snowg(i) / ht(i)
        dfi(i) = min (dfi(i), 1.-fi(i))
  200 continue
c
      do 210 k=1,nsnolay
        do 212 i=1,npoi
          if (fi(i)+dfi(i).gt.0.)
     >      tsno(i,k) = (tsno(i,k)*fi(i) + tsnowg(i)*dfi(i))
     >                / (fi(i)+dfi(i))
c
c set initial thicknesses for newly created snow
c
          if (fi(i).eq.0. .and. dfi(i).gt.0.) hsno(i,k) = hinit(k)
  212   continue
  210 continue
c
      npn = 0
      do 220 i=1,npoi
        fi(i) = fi(i) + dfi(i)
        if (fi(i).gt.0.) then 
          npn = npn + 1
          indsno(npn) = i
        endif
  220 continue
c
c melt from any layer (due to implicit heat conduction, any
c layer can exceed tmelt, not just the top layer), and reduce
c thicknesses (even to zero, and give extra heat to soil)
c
c ok to do it for non-snow points, for which xh = xm = 0
c
      if (npn.gt.0) then
c
        rwork = 1. / rhos
        do 300 k=1,nsnolay
          do 302 i=1,npoi
            xh(i) = rhos*hsno(i,k)*cice * max(tsno(i,k)-tmelt, 0.)
            xm(i) = min (rhos*hsno(i,k), xh(i)/hfus)
            hsno(i,k) = hsno(i,k) - xm(i)*rwork
            tsno(i,k) = min (tsno(i,k),tmelt)
            zmelt(i) = zmelt(i) + fi(i)*xm(i)
            zheat(i) = zheat(i) + fi(i)*(xh(i)-hfus*xm(i))
  302     continue
  300   continue
c
c adjust fi and thicknesses for coverage-vs-volume relation
c ie, total thickness = hsnomin for fi < fimax, and fi <= fimax.
c (ok to do it for no-snow points, for which ht=fi=finew=0.)
c
        call const (ht, npoi, 0.0)
        do 400 k=1,nsnolay
          do 402 i=1,npoi
            ht(i) = ht(i) + hsno(i,k)
  402     continue
  400   continue
c
c linear variation  for 0 < fi < 1
c
        zhh = 0.5*hsnomin
        do 404 i=1,npoi
          zdh = hsnoruf(i)-hsnomin
          finew = ( -zhh + sqrt(zhh**2 + zdh*fi(i)*ht(i)) ) / zdh

          finew = max (0., min (fimax, finew))
          x1(i) =  fi(i) / max (finew, epsilon)
          fi(i) =  finew
  404   continue
c
        do 406 k=1,nsnolay
          do 408 i=1,npoi
            hsno(i,k) = hsno(i,k) * x1(i)
  408     continue
  406   continue
c
      endif
c
c re-adapt snow thickness profile, so top thickness = hsnotop
c and other thicknesses are equal
c
c adjust temperature to conserve sensible heat
c
      call vadapt (hsno, tsno, hsnotop, indsno, npn, nsnolay)
c
c if fi is below fimin, melt all snow and adjust soil fluxes
c
      if (npn.gt.0) then
        call scopy (npoi, fi, x1)
        do 500 k=1,nsnolay
          do 502 i=1,npoi
            if (x1(i).lt.fimin) then
              xm(i) = x1(i) * rhos * hsno(i,k)
              zmelt(i) = zmelt(i) + xm(i)
              zheat(i) = zheat(i) - xm(i)*(cice*(tmelt-tsno(i,k))+hfus)
              hsno(i,k) = 0.
              tsno(i,k) = tmelt
              fi(i) = 0.
            endif
  502     continue
  500   continue
      endif
c
c adjust buried lower veg for fi changes. if fi has increased,
c incorporate newly buried intercepted h2o into bottom-layer 
c snow, giving associated heat increment to soil, and mix the
c specific heat of newly buried veg (at tl) into tlsub. if fi
c has decreased, change temp of newly exhumed veg to tl, giving
c assoc heat increment to soil, and smear out intercepted h2o
c
      if (npn.gt.0) then
        do 600 i=1,npoi
          dfi(i) = fi(i) - fiold(i)
c
          if (dfi(i).gt.0.) then
c
c factor of xl*chl has been divided out of next line
c
            tlsub(i)= (tlsub(i)*fiold(i) + tl(i)*dfi(i)) / fi(i)
            zheat(i) = zheat(i) + dfi(i)*xl(i)
     >               * ( wliql(i) * (ch2o*(tl(i)-tmelt) + hfus
     >                              +cice*(tmelt-tsno(i,nsnolay)))
     >                   + wsnol(i) *  cice*(tl(i)-tsno(i,nsnolay)) )
c
            hsno(i,nsnolay) = hsno(i,nsnolay)
     >                      + dfi(i)*xl(i)*(wliql(i)+wsnol(i))
     >                        / (rhos*fi(i))
          endif
c
          if (dfi(i).lt.0.) then
            zheat(i) = zheat(i) - dfi(i)*xl(i)*chl*(tlsub(i)-tl(i))
            rwork = (1.-fiold(i)) / (1.-fi(i))
            wliql(i) = wliql(i) * rwork
            wsnol(i) = wsnol(i) * rwork
          endif
c
  600   continue
      endif
c
c areally average fluxes to be used by soil model. (don't use
c index due to mix call, but only need at all if npn > 0)
c
      if (npn.gt.0) then
c
        rwork = 1. / dtime
        do 700 i=1,npoi
          rwork2 = 1. - fiold(i)
          heatg(i) = rwork2*heatg(i)
     >             + fiold(i)*sflo(i,nsnolay+2)
     >             + zheat(i)*rwork
          solg(i)  = rwork2 * solg(i)
          fvapg(i) = rwork2 * fvapg(i)
          x1(i)    = rwork2 * raing(i)
          x2(i)    = zmelt(i)*rwork
          x3(i)    = tmelt
  700   continue
c
        call mix (raing,traing, x1,traing, x2,x3, vzero,vzero)
c
      endif
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine snowheat (tlsub, fhtop, sflo, xl, chl)
c ---------------------------------------------------------------------
c
c sets up call to tridia to solve implicit snow heat conduction,
c using snow temperatures in tsno (in comsno). adds an extra
c buried-lower-veg layer to the bottom of the snow with 
c conduction coefficient conbur/xl and heat capacity chl*xl
c

      use implicit
c
      use compar
      use comsno
      use comsoi
c
c Arguments
c
      real chl                   ! specific heat of lower veg per l/s area (supplied)
      
      real tlsub(npoi),          ! temperature of buried lower veg (supplied, returned)
     >     fhtop(npoi),          ! heat flux into top snow layer from atmos (supplied)
     >     sflo(npoi,nsnolay+2), ! downward heat flow across layer boundaries (returned)
     >     xl(npoi)              ! (lai(i,1)+sai(i,1))*fl(i), lower-veg density(supplied)
c
c Local variables
c
      integer k, i,              ! loop indices
     >        km1,               ! used to avoid layer 0
     >        kp1                ! used to avoid layer nsnolay+2
c
      real rimp,                 ! implicit fraction of the calculation (0 to 1)
     >     conbur,               ! conduction coeff of buried lower veg layer 
     >                           ! for unit density xl=(lai+sai)*fl, in w m-2 k-1
     >     hfake,                ! arbitrary small thickness to allow processing 
     >                           ! for zero snow. (doesn't use index since tridia
     >                           ! not set up for index.)
     >     rwork,                ! to compute matrix diagonals and right-hand side
     >     dt,                   ! '
     >     dti                   ! '
c
      real con(npoi,nsnolay+2),  ! conduction coefficents between layers
     >     temp(npoi,nsnolay+1), ! combined snow and buried-veg temperatures
     >     d1(npoi,nsnolay+1),   ! diagonals of tridiagonal systems of equations 
     >     d2(npoi,nsnolay+1),   ! '
     >     d3(npoi,nsnolay+1),   ! '
     >     rhs(npoi,nsnolay+1),  ! right-hand sides of systems of equations
     >     w1(npoi,nsnolay+1),   ! work array needed by tridia
     >     w2(npoi,nsnolay+1)    ! '
c
c conbur (for xl=1) is chosen to be equiv to 10 cm of snow
c
      data rimp, conbur, hfake /1.0, 2.0, .01/
c
c copy snow and buried-lower-veg temperatures into combined
c array temp
c
      call scopy (npoi*nsnolay, tsno,  temp             )
      call scopy (npoi,         tlsub, temp(1,nsnolay+1))
c
c set conduction coefficients between layers
c
      do 100 k=1,nsnolay+2
        if (k.eq.1) then
          call const (con(1,k), npoi, 0.0)
c
        else if (k.le.nsnolay) then
          rwork = 0.5 / consno
          do 102 i=1,npoi
            con(i,k) = 1. / (   max(hsno(i,k-1),hfake)*rwork
     >                        + max(hsno(i,k)  ,hfake)*rwork )
  102     continue
c
        else if (k.eq.nsnolay+1) then
          rwork = 0.5 / consno
          do 104 i=1,npoi
            con(i,k) = 1. / (   max(hsno(i,k-1),hfake)*rwork
     >                        + 0.5*xl(i)/conbur )
  104     continue
c
        else if (k.eq.nsnolay+2) then
          rwork = 0.5 / conbur
          do 106 i=1,npoi
            con(i,k) = 1. / (   xl(i)*rwork
     >                        + 0.5*hsoi(1) / consoi(i,1) )
  106     continue
        endif
  100 continue
c
c set matrix diagonals and right-hand side. for layer nsnolay+1
c (buried-lower-veg layer), use explicit contact with soil, and
c multiply eqn through by xl*chl/dtime to allow zero xl.
c
      do 200 k=1,nsnolay+1
        km1 = max (k-1,1)
        kp1 = min (k+1,nsnolay+1)
c
        if (k.le.nsnolay) then
          rwork = dtime /(rhos*cice)
          do 202 i=1,npoi
            dt = rwork / (max(hsno(i,k),hfake))
            d1(i,k) =    - dt*rimp* con(i,k)
            d2(i,k) = 1. + dt*rimp*(con(i,k)+con(i,k+1))
            d3(i,k) =    - dt*rimp* con(i,k+1)
c
            rhs(i,k) = temp(i,k) + dt
     >               * ( (1.-rimp)*con(i,k)  *(temp(i,km1)-temp(i,k))
     >                 + (1.-rimp)*con(i,k+1)*(temp(i,kp1)-temp(i,k)) )
  202     continue
c
          if (k.eq.1) then 
            rwork = dtime /(rhos*cice)
            do 204 i=1,npoi
              dt = rwork / (max(hsno(i,k),hfake))
              rhs(i,k) = rhs(i,k) + dt*fhtop(i)
  204       continue
          endif
c
        else if (k.eq.nsnolay+1) then
c
          rwork = chl / dtime
          do 206 i=1,npoi
            dti = xl(i)*rwork
            d1(i,k) =     -  rimp* con(i,k)
            d2(i,k) = dti +  rimp*(con(i,k)+con(i,k+1))
            d3(i,k) = 0.
            rhs(i,k) = dti*temp(i,k)
     >               + ( (1.-rimp)*con(i,k)*(temp(i,km1)-temp(i,k))
     >                 + con(i,k+1)*(tsoi(i,1)-(1.-rimp)*temp(i,k)) )
  206     continue
        endif
  200 continue
c
c solve the tridiagonal systems
c
      call tridia (npoi, npoi, nsnolay+1, d1,d2,d3, rhs, temp, w1,w2)
c
c deduce downward heat fluxes between layers
c
      call scopy (npoi, fhtop, sflo(1,1))
c
      do 400 k=1,nsnolay+1
        if (k.le.nsnolay) then
          rwork = rhos*cice/dtime
          do 402 i=1,npoi
            sflo(i,k+1) = sflo(i,k) - rwork*hsno(i,k)
     >                                *(temp(i,k)-tsno(i,k))
  402     continue
c
        else
          rwork = chl/dtime
          do 404 i=1,npoi
            sflo(i,k+1) = sflo(i,k)
     >                  - xl(i)*rwork*(temp(i,nsnolay+1)-tlsub(i))
  404     continue
        endif
  400 continue
c
c copy temperature solution to tsno and tlsub, but not for
c points with no snow
c
      do 500 k=1,nsnolay
        do 502 i=1,npoi
          if (fi(i).gt.0.) tsno(i,k) = temp(i,k) 
  502   continue
  500 continue
c
      do 510 i=1,npoi
        if (fi(i).gt.0.) tlsub(i) = temp(i,nsnolay+1)
  510 continue
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine vadapt (hcur, tcur, htop, indp, np, nlay) 
c ---------------------------------------------------------------------
c
c re-adapt snow layer thicknesses, so top thickness
c equals hsnotop and other thicknesses are equal
c
c also adjusts profile of tracer field tcur so its vertical
c integral is conserved (eg, temperature)
c
      use implicit
c
      use compar
c
c Arguments
c
      integer np,            ! number of snow pts in current strip (supplied)
     >        nlay           ! # of layer
c
      integer indp(npoi)     ! index of snow pts in current strip (supplied)
c
      real htop              ! prescribed top layer thickness (supplied)
c
      real hcur(npoi,nlay),  ! layer thicknesses (supplied and returned)     
     >     tcur(npoi,nlay)   ! tracer field (supplied and returned)
c
c local variables
c
      integer i, j, k, ko    ! loop indices
c
      real dz, rwork
c
      real ht(npoi),         ! storing variable for zold        
     >     h1(npoi),         ! to compute new layer thickness
     >     za(npoi),         ! 
     >     zb(npoi),         ! 
     >     zheat(npoi)       !
c
      real hnew(npoi,nsnolay),    ! new layer thickness
     >     tnew(npoi,nsnolay),    ! new temperatures of layers
     >     zold(npoi,nsnolay+1)   ! distances from surface to old layer boundaries
c
c if no snow or seaice points in current 1d strip, return. note
c that the index is not used below (for cray vec and efficiency)
c except in the final loop setting the returned values
c
      if (np.eq.0) return
c
c set distances zold from surface to old layer boundaries
c
      call const (zold(1,1), npoi, 0.0)
c
      do 300 k=1,nlay
        do 302 i=1,npoi
          zold(i,k+1) = zold(i,k) + hcur(i,k)
  302   continue
  300 continue
c
c set new layer thicknesses hnew (tot thickness is unchanged).
c if total thickness is less than nlay*htop (which should be
c le hsnomin), make all new layers equal including
c top one, so other layers aren't so thin. use epsilon to 
c handle zero (snow) points
c
      call scopy (npoi, zold(1,nlay+1), ht)
c
      rwork = nlay*htop
      do 304 i=1,npoi
        if (ht(i).ge.rwork) then
          h1(i) = (ht(i)-htop)/(nlay-1)
        else
          h1(i) = max (ht(i)/nlay, epsilon)
        endif
  304 continue
c
      do 306 k=1,nlay
        do 308 i=1,npoi
          hnew(i,k) = h1(i)
  308   continue
  306 continue
c
      rwork = nlay*htop
      do 310 i=1,npoi
        if (ht(i).ge.rwork) hnew(i,1) = htop
  310 continue
c
c integrate old temperature profile (loop 410) over each
c new layer (loop 400), to get new field tnew
c
      call const (zb, npoi, 0.0)
c
      do 400 k=1,nlay
c
        do 402 i=1,npoi
          za(i) = zb(i)
          zb(i) = za(i) + hnew(i,k)
  402   continue
        call const (zheat, npoi, 0.0)
c
        do 410 ko=1,nlay
          do 412 i=1,npoi
            if (za(i).lt.zold(i,ko+1) .and. zb(i).gt.zold(i,ko)) then
              dz = min(zold(i,ko+1),zb(i)) - max(zold(i,ko),za(i))
              zheat(i) = zheat(i) + tcur(i,ko)*dz
            endif
  412     continue
  410   continue
c
        do 420 i=1,npoi
          tnew(i,k) = zheat(i) / hnew(i,k)
  420   continue
c
  400 continue
c
c use index for final copy to seaice or snow arrays, to avoid
c changing soil values (when called for seaice) and to avoid
c changing nominal snow values for no-snow points (when called
c for snow)
c
      do 500 k=1,nlay
        do 502 j=1,np
          i = indp(j)
          hcur(i,k) = hnew(i,k)
          tcur(i,k) = tnew(i,k)
  502   continue
  500 continue
c
      return
      end
