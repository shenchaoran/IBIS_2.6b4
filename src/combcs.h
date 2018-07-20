c
c ------
c combcs
c ------
c
      real
     >  xintopo(npoi),    ! topography (m)
     >  xinveg(npoi),     ! fixed vegetation map
     >  deltat(npoi)      ! absolute minimum temperature -
     >                    ! temp on average of coldest month (C)
c
c     common /combcs1/ xintopo, xinveg, deltat
c
      integer
     >  lmask(nlon,nlat)  ! landmask 0=water, 1=land
c
c     common /combcs2/ lmask
c
      real
     >  xint(npoi,12),    ! climatological temp + anomaly (C)
     >  xinq(npoi,12),    ! climatological relative humidity + anomaly (%)
     >  xinprec(npoi,12), ! climatological precipition + anomaly (mm/day)
     >  xinwind(npoi,12), ! climatological wind speed + anomaly (m s-1)
     >  xincld(npoi,12),  ! climatological cloudiness + anomaly(%)
     >  xinwet(npoi,12),  ! climatological wet days + anomaly (days/month)
     >  xintrng(npoi,12)  ! climatological temp range + anomaly(C)
c
c     common /combcs3/ xint, xinq, xinprec, xinwind, xincld, xinwet, xintrng
c
      real
     >  clmt(npoi,12),    ! climatological temperature (C)
     >  clmq(npoi,12),    ! climatological relative humidity (%)
     >  clmprec(npoi,12), ! climatological precipitation (mm/day)
     >  clmw(npoi,12),    ! climatological wind speed (m s-1)
     >  clmwet(npoi,12),  ! climatological wet days (days/month)
     >  clmcld(npoi,12),  ! climatological cloudiness (%)
     >  clmtrng(npoi,12)  ! climatological temp range (C)
c
c     common /combcs4/ clmt, clmq, clmprec, clmw, clmwet, clmcld, clmtrng
c
      real
     >  xintd(npoi),      ! daily climatological temperature (C) +
     >                    !   anomaly (C) + daily anomaly (C)
     >  xinqd(npoi),      ! daily climatological relative humidity (%) +
     >                    !   anomaly (%) * daily anomaly (fraction)
     >  xinprecd(npoi),   ! daily climatological precipitation (mm/day) +
     >                    !   anomaly (mm/day) * daily anomaly (fraction)
     >  xinwindd(npoi),   ! daily climatological windspeed (m/s) +
     >                    !   anomaly (m/s) * daily anomaly (fraction)
     >  xincldd(npoi),    ! daily climatological cloud fraction (fraction) +
     >                    !   anomaly (fraction) + daily anomaly (fraction)
     >  xintrngd(npoi)    ! daily climatological temp range (C) +
     >                    !   anomaly (C) * daily anomaly (C)
c
c     common /combcs5/ xintd, xinqd, xinprecd, xinwindd, xincldd, xintrngd
