c
c  ####   #          #    #    #    ##     #####  ######
c #    #  #          #    ##  ##   #  #      #    #
c #       #          #    # ## #  #    #     #    #####
c #       #          #    #    #  ######     #    #
c #    #  #          #    #    #  #    #     #    #
c  ####   ######     #    #    #  #    #     #    ######
c
c ---------------------------------------------------------------------
      subroutine climanl
c ---------------------------------------------------------------------
c
c this subsroutine is only used to initialize growing degree days,
c coldest temp, and warmest temp at very beginning - provides a
c climate 'history' based on monthly mean values
c
c common blocks
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'combcs.h'
      include 'comveg.h'
c
c Local variables
c
      integer it1w,      ! indice of previous month (interpolation)
     >        it2w,      ! indice of following month (interpolation)
     >      i,k,lda      ! loop indices
c
      real rwork,        ! work variable (1/ndaypm)
     >     dt,           ! used for interpolation
     >     dtemp         ! interpolated temperature
c
c initialize values
c
      call const (gdd0, npoi, 0.0)
      call const (gdd5, npoi, 0.0)
c
      do 100 i = 1, npoi
c
c coldest monthly temperature (year 0) in deg c
c
        tc(i) = min (xint(i,1),  xint(i,2),  xint(i,3),
     >               xint(i,4),  xint(i,5),  xint(i,6),
     >               xint(i,7),  xint(i,8),  xint(i,9),
     >               xint(i,10), xint(i,11), xint(i,12))
c
c warmest monthly temperature (year 0) in deg c
c
        tw(i) = max (xint(i,1),  xint(i,2),  xint(i,3),
     >               xint(i,4),  xint(i,5),  xint(i,6),
     >               xint(i,7),  xint(i,8),  xint(i,9),
     >               xint(i,10), xint(i,11), xint(i,12))
c
        tcmin(i) = tc(i) + deltat(i)
c
 100  continue 
c
c interpolating climatological monthly input values to daily
c
      do 200 i = 1, npoi 
c
        do 210 k = 1, 12
c
          rwork = 1. / float(ndaypm(k))
c
          do 220 lda = 1, ndaypm(k)
c
            if (float(lda).lt.float(ndaypm(k)+1)*0.5) then
              it1w = k - 1
              it2w = k
              dt   = (float(lda) - 0.5) * rwork + 0.5
            else
              it1w = k
              it2w = k + 1
              dt   = (float(lda) - 0.5) * rwork - 0.5
            end if
c
            if (it1w.lt. 1) it1w = 12
            if (it2w.gt.12) it2w = 1
c
            dtemp = xint(i,it1w) +
     >              dt * (xint(i,it2w) - xint(i,it1w))
c
c growing degree days, using deg c
c
            gdd0(i) = gdd0(i) + max(0.0, dtemp)
            gdd5(i) = gdd5(i) + max(0.0, (dtemp - 5.0))
c
 220      continue
c
 210    continue
c
 200  continue
c
c call routine to determine pft existence arrays
c
      call existence
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine climanl2
c ---------------------------------------------------------------------
c
c this subroutine updates the growing degree days, coldest temp, and
c warmest temp if monthly anomalies or daily values are used
c
c common blocks
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'combcs.h'
      include 'comveg.h'
c
c local variables
c
      integer i             ! loop indice
c
      real zweigc,          ! 30-year e-folding time-avarage
     >     zweigw,          ! 30-year e-folding time-avarage
     >     rworkc,          ! 30-year e-folding time-avarage
     >     rworkw 
c
c calculate a 30-year e-folding time-avarage
c
      zweigc = exp(-1./30.)
      zweigw = exp(-1./30.)
c
      rworkc = 1. - zweigc
      rworkw = 1. - zweigw
c
c update critical climatic parameters with running average
c
      do 100 i = 1, npoi
c
        tc(i) = zweigc * tc(i) + rworkc * tcthis(i)
        tw(i) = zweigw * tw(i) + rworkw * twthis(i)
c
        tcmin(i) = tc(i) + deltat(i)
c
        gdd0(i) = zweigc * gdd0(i) +
     >            rworkc * gdd0this(i)
c
        gdd5(i) = zweigc * gdd5(i) +
     >            rworkc * gdd5this(i)
c
 100  continue
c
      call existence
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine existence
c ---------------------------------------------------------------------
c
c this routine determines which plant functional types (pft's) are allowed
c to exist in each gridcell, based on a simple set of climatic criteria
c
c the logic here is based on the biome3 model of haxeltine and prentice
c
c plant functional types:
c
c 1)  tropical broadleaf evergreen trees
c 2)  tropical broadleaf drought-deciduous trees
c 3)  warm-temperate broadleaf evergreen trees
c 4)  temperate conifer evergreen trees
c 5)  temperate broadleaf cold-deciduous trees
c 6)  boreal conifer evergreen trees
c 7)  boreal broadleaf cold-deciduous trees
c 8)  boreal conifer cold-deciduous trees
c 9)  evergreen shrubs
c 10) deciduous shrubs
c 11) warm (c4) grasses
c 12) cool (c3) grasses
c
c
c common blocks
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'comveg.h'
      include 'compft.h'
c
c Local variables
c
      integer i      ! loop indice
c
c ---------------------------------------------------------------------
c
      do 100 i = 1, npoi
c
c determine which plant types can exist in a given gridcell
c
        exist(i,1)  = 0.
        exist(i,2)  = 0.
        exist(i,3)  = 0.
        exist(i,4)  = 0.
        exist(i,5)  = 0.
        exist(i,6)  = 0.
        exist(i,7)  = 0.
        exist(i,8)  = 0.
        exist(i,9)  = 0.
        exist(i,10) = 0.
        exist(i,11) = 0.
        exist(i,12) = 0.
c
c 1) tropical broadleaf evergreen trees
c
c  - tcmin > 0.0
c
*        if (tcmin(i).gt.0.0)           exist(i,1) = 1.0
c
c 2) tropical broadleaf drought-deciduous trees
c
c  - tcmin > 0.0
c
*        if (tcmin(i).gt.0.0)           exist(i,2) = 1.0
c
c 3) warm-temperate broadleaf evergreen trees
c
c  - tcmin <   0.0 and
c  - tcmin > -10.0
c
*        if ((tcmin(i).lt.0.0).and.
*     >      (tcmin(i).gt.-10.0))       exist(i,3) = 1.0
c
c 4) temperate conifer evergreen trees
c
c  - tcmin <    0.0 and
c  - tcmin >  -45.0 and
c  - gdd5  > 1200.0
c
*        if ((tcmin(i).lt.0.0).and.
*     >      (tcmin(i).gt.-45.0).and.
*     >      (gdd5(i).gt.1200.0))       exist(i,4) = 1.0
c
c 5) temperate broadleaf cold-deciduous trees
c
c  - tcmin <    0.0 and
c  - tcmin >  -45.0 and
c  - gdd5  > 1200.0
c
*        if ((tcmin(i).lt.0.0).and.
*     >      (tcmin(i).gt.-45.0).and.
*     >      (gdd5(i).gt.1200.0))       exist(i,5) = 1.0
c
c 6) boreal conifer evergreen trees
c
c  - tcmin <  -45.0 or gdd5 < 1200.0, and
c  - tcmin >  -57.5 and
c  - gdd5  >  350.0
c
*        if (((tcmin(i).lt.-45.0).or.(gdd5(i).lt.1200.0)).and.
*     >       (tcmin(i).gt.-57.5).and.
*     >       (gdd5(i).gt.350.0))       exist(i,6) = 1.0
c
c 7) boreal broadleaf cold-deciduous trees
c
c  - tcmin <  -45.0 or gdd5 < 1200.0, and
c  - tcmin >  -57.5 and
c  - gdd5  >  350.0
c
*        if (((tcmin(i).lt.-45.0).or.(gdd5(i).lt.1200.0)).and.
*     >       (tcmin(i).gt.-57.5).and.
*     >       (gdd5(i).gt.350.0))       exist(i,7) = 1.0
c
c 8) boreal conifer cold-deciduous trees
c
c  - tcmin <  -45.0 or gdd5 < 1200.0, and
c  - gdd5  >  350.0
c
*        if (((tcmin(i).lt.-45.0).or.(gdd5(i).lt.1200.0)).and.
*     >       (gdd5(i).gt.350.0))       exist(i,8) = 1.0
c
c 9) evergreen shrubs
c
c  - gdd0 > 100.0
c
*        if (gdd0(i).gt.100.0)          exist(i,9) = 1.0
c
c 10) deciduous shrubs
c
c  - gdd0 > 100.0
c
*        if (gdd0(i).gt.100.0)          exist(i,10) = 1.0
c
c 11) warm (c4) grasses
c
c  - tw   >  22.0 and
c  - gdd0 > 100.0
c
*        if ((tw(i).gt.22.0).and.
*     >      (gdd0(i).gt.100.0))        exist(i,11) = 1.0
c
c 12) cool (c3) grasses
c
c  - gdd0 > 100.0
c
*        if (gdd0(i).gt.100.0)          exist(i,12) = 1.0
c
c
**** DTP 2001/06/07: Modified version of above code reads in PFT
*    existence criteria from external parameter file "params.veg"
*    These are copied here for reference.... 
*------------------------------------------------------------------
*  TminL    TminU    Twarm    GDD    PFT
*------------------------------------------------------------------
*    0.0   9999.0   9999.0   9999  !   1
*    0.0   9999.0   9999.0   9999  !   2
*  -10.0      0.0   9999.0   9999  !   3
*  -45.0      0.0   9999.0   1200  !   4
*  -45.0      0.0   9999.0   1200  !   5
*  -57.5    -45.0   9999.0    350  !   6
*  -57.5    -45.0   9999.0    350  !   7
* 9999.0    -45.0   9999.0    350  !   8
* 9999.0   9999.0   9999.0    100  !   9
* 9999.0   9999.0   9999.0    100  !  10
* 9999.0   9999.0     22.0    100  !  11
* 9999.0   9999.0   9999.0    100  !  12
*------------------------------------------------------------------

c 1) tropical broadleaf evergreen trees
c
c  - tcmin > 0.0
c
        if (tcmin(i).gt.TminL(1))      exist(i,1) = 1.0
c
c 2) tropical broadleaf drought-deciduous trees
c
c  - tcmin > 0.0
c
        if (tcmin(i).gt.TminL(2))      exist(i,2) = 1.0
c
c 3) warm-temperate broadleaf evergreen trees
c
c  - tcmin <   0.0 and
c  - tcmin > -10.0
c
        if ((tcmin(i).lt.TminU(3)).and.
     >      (tcmin(i).gt.TminL(3)))    exist(i,3) = 1.0
c
c 4) temperate conifer evergreen trees
c
c  - tcmin <    0.0 and
c  - tcmin >  -45.0 and
c  - gdd5  > 1200.0
c
        if ((tcmin(i).lt.TminU(4)).and.
     >      (tcmin(i).gt.TminL(4)).and.
     >      (gdd5(i).gt.GDD(4)))       exist(i,4) = 1.0
c
c 5) temperate broadleaf cold-deciduous trees
c
c  - tcmin <    0.0 and
c  - tcmin >  -45.0 and
c  - gdd5  > 1200.0
c
        if ((tcmin(i).lt.TminU(5)).and.
     >      (tcmin(i).gt.TminL(5)).and.
     >      (gdd5(i).gt.GDD(5)))       exist(i,5) = 1.0
c
c 6) boreal conifer evergreen trees
c
c  - tcmin <  -45.0 or gdd5 < 1200.0, and
c  - tcmin >  -57.5 and
c  - gdd5  >  350.0
c
        if (((tcmin(i).lt.TminU(6)).or.
     >      (gdd5(i).lt.GDD(4))).and.
     >      (tcmin(i).gt.TminL(6)).and.
     >      (gdd5(i).gt.GDD(6)))       exist(i,6) = 1.0
c
c 7) boreal broadleaf cold-deciduous trees
c
c  - tcmin <  -45.0 or gdd5 < 1200.0, and
c  - tcmin >  -57.5 and
c  - gdd5  >  350.0
c
        if (((tcmin(i).lt.TminU(7)).or.
     >      (gdd5(i).lt.GDD(5))).and.
     >      (tcmin(i).gt.TminL(7)).and.
     >      (gdd5(i).gt.GDD(7)))       exist(i,7) = 1.0
c
c 8) boreal conifer cold-deciduous trees
c
c  - tcmin <  -45.0 or gdd5 < 1200.0, and
c  - gdd5  >  350.0
c
        if (((tcmin(i).lt.TminU(8)).or.
     >      (gdd5(i).lt.TminL(4))).and.
     >      (gdd5(i).gt.GDD(8)))       exist(i,8) = 1.0
c
c 9) evergreen shrubs
c
c  - gdd0 > 100.0
c
        if (gdd0(i).gt.GDD(9))         exist(i,9) = 1.0
c
c 10) deciduous shrubs
c
c  - gdd0 > 100.0
c
        if (gdd0(i).gt.GDD(10))        exist(i,10) = 1.0
c
c 11) warm (c4) grasses
c
c  - tw   >  22.0 and
c  - gdd0 > 100.0
c
        if ((tw(i).gt.Twarm(11)).and.
     >      (gdd0(i).gt.GDD(11)))      exist(i,11) = 1.0
c
c 12) cool (c3) grasses
c
c  - gdd0 > 100.0
c
        if (gdd0(i).gt.GDD(12))        exist(i,12) = 1.0

 100  continue
c
      return
      end
c
