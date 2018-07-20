c
c #    #  ######    ##     #####  #    #  ######  #####
c #    #  #        #  #      #    #    #  #       #    #
c #    #  #####   #    #     #    ######  #####   #    #
c # ## #  #       ######     #    #    #  #       #####
c ##  ##  #       #    #     #    #    #  #       #   #
c #    #  ######  #    #     #    #    #  ######  #    #
c
c ---------------------------------------------------------------------
      subroutine daily (iyear, imonth, iday, seed, jdaily, iyranom, iyrlast, nrun)
c ---------------------------------------------------------------------
c
c overview
c
c this routine generates daily weather conditions from monthly-mean
c climatic parameters
c
c specifically, this routine generates daily values of
c
c  - daily total precipitation
c  - daily maximum temperature
c  - daily minimum temperature
c  - daily average cloud cover
c  - daily average relative humidity
c  - daily average wind speed
c
c in order to generate daily weather conditions, the model uses a series
c of 'weather generator' approaches, which generate random combinations of
c weather conditions based upon the climatological conditions
c
c in general, this weather generator is based upon the so-called Richardson
c weather generator
c
c appropriate references include:
c
c Geng, S., F.W.T. Penning de Vries, and L. Supit, 1985:  A simple
c method for generating rainfall data, Agricultural and Forest
c Meteorology, 36, 363-376.
c
c Richardson, C. W. and Wright, D. A., 1984: WGEN: A model for 
c generating daily weather variables: U. S. Department of
c Agriculture, Agricultural Research Service.
c
c Richardson, C., 1981: Stochastic simulation of daily
c precipitation, temperature, and solar radiation. Water Resources 
c Research 17, 182-190.
c
c common blocks
c
      use implicit
c
      use compar
      use combcs
      use comatm
      use comsum
      use comveg
c
c Arguments
c
      integer seed, 
     >        iyear,
     >        imonth,
     >        iday,
     >        jdaily   ! 1 if reading in daily weather data
     >                 ! 0 if using random/statistical weather generator      
c
c local variables
c
      integer it1w,      ! indice of previous month (interpolation)
     >        it2w,      ! indice of following month (interpolation)
     >        i,j, k     ! loop indice
c
      real rwork,          ! 
     >     omcloud,        ! cloud cover
     >     omqd,           ! humidity
     >     omtmax,         ! maximum temperature
     >     ran2,           ! function random number generator
     >     dt,             ! used for interpolation
     >     pwet,           ! monthly-average probability of rainy day
     >     pwd,            ! probability of a wet day after a dry day
     >     pww,            ! probability of a wet day after a wet day
     >     rndnum,         ! random number to decide if wet or dry day
     >     rainpwd,        ! average rainfall per wet day
     >     alpha,          ! parameter for gamma function
     >     beta,           ! parameter for gamma function
     >     aa,
     >     ab,
     >     tr1,
     >     tr2,
     >     rn1,rn2, rn3,rn,  !random numbers
     >     s1,
     >     s2,
     >     s12,
     >     z,
     >     tdm,            ! mean daily mean temperature
     >     trngm,          ! mean daily mean temperature
     >     tmaxm,          ! mean maximum temperature
     >     tminm,          ! mean minimum temperature
     >     tmaxd,          ! maximum temperatures for dry days
     >     tmaxw,          ! maximum temperatures for wet days
     >     tmaxe,          !'expected' maximum temperature for today
     >     tmine,          !'expected' minimum temperature for today
     >     tmaxs,          ! standard deviation in maximum temperature (K)
     >     tmins,          ! standard deviation in minimum temperature (K)
     >     cloudm,         ! mean cloud cover for today (fraction)
     >     cloudd,         ! dry day cloud cover
     >     cloudw,         ! wet day cloud cover
     >     cloude,         ! expected cloud cover today
     >     clouds,         ! standard deviation of cloud fraction
     >     v,
     >     tdum,           ! storage variable
     >     qdm,            ! mean relative humidity
     >     qdd,            ! dry day relative humidity
     >     qdw,            ! wet day relative humidity 
     >     qde,            ! expected relative humidity (based on wet/dry decision)
     >     qdup,           ! upper bound of humidity distribution function
     >     qdlow,          ! lower bound of humidity distribution function
     >     y,
     >     b3,
     >     b2,
     >     b1,
     >     x1,
     >     amn,
     >     eud             ! expected daily average wind speed from monthly mean
 
c
      real a(3,3),
     >     b(3,3)
c
      real
     >     ee(3), 
     >     r(3),
     >     rr(3),
     >     x(3)
c
      integer
     >     iyranom,
     >     iyrlast,
     >     nrun
c
      real
     >     precipfac,
     >     dif,
     >     humidfrac,
     >     sphumid
c
c define autocorrelation matrices for Richardson generator
c
c note that this matrix should be based upon a statistical
c analysis of regional weather patterns
c
c for global simulations, we use 'nominal' values
c
      data a /  0.600,  0.500,  0.005,
     >          0.010,  0.250,  0.005, 
     >          0.020,  0.125,  0.250 /
c
      data b /  0.500,  0.250, -0.250,
     >          0.000,  0.500,  0.250, 
     >          0.000,  0.000,  0.500 /
c
      use comsat
c
c ---------------------------------------------------------------------- 
c * * * initial setup for daily climate calculations * * *
c ---------------------------------------------------------------------- 
c
c define working variables
c
      rwork = (grav / rair / 0.0065)
c
c 'omega' parameters used to calculate differences in expected
c climatic parameters on wet and dry days
c
c following logic of weather generator used in the EPIC crop model
c
c omcloud -- cloud cover
c omqd    -- humidity
c omtmax  -- maximum temperature
c
      omcloud = 0.90    ! originally 0.90
      omqd    = 0.50    ! originally 0.50
      omtmax  = 0.75    ! originally 0.75
c
c calculate weighting factors used in interpolating climatological
c monthly-mean input values to daily-mean values
c
c this is a simple linear interpolation technique that takes into
c account the length of each month
c
      if (jdaily .eq. 0) then
        if (float(iday).lt.float(ndaypm(imonth)+1)/2.0) then
          it1w = imonth - 1
          it2w = imonth 
          dt   = (float(iday) - 0.5) / ndaypm(imonth) + 0.5
        else
          it1w = imonth 
          it2w = imonth + 1
          dt   = (float(iday) - 0.5) / ndaypm(imonth) - 0.5
        end if
c
        if (it1w.lt. 1) it1w = 12
        if (it2w.gt.12) it2w = 1
      end if
c
c initialize this year's values of gdd0, gdd5, tc, tw
c
      if (iday.eq.1 .and. imonth.eq.1) then
c
        call const (tcthis, npoi,  100.0)
        call const (twthis, npoi, -100.0)
c
        call const (gdd0this, npoi, 0.0)
        call const (gdd5this, npoi, 0.0)
c
      endif
c
c ---------------------------------------------------------------------- 
c * * * set daily climatic variables for entire domain * * *
c ---------------------------------------------------------------------- 
c
      do 200 i = 1, npoi
c
c ---------------------------------------------------------------------- 
c * * * use weather generator to create daily statistics * * *
c ---------------------------------------------------------------------- 
c
        if (jdaily .eq. 0) then
c
c ---------------------------------------------------------------------- 
c (1) determine if today will rain or not (following Geng et al.)
c ---------------------------------------------------------------------- 
c
c implement simple first-order Markov-chain precipitation generator logic
c based on Geng et al. (1986), Richardson and Wright (1984),
c and Richardson (1981) 
c
c basically, this allows for the probability of today being a wet day 
c (a day with measureable precipitation) to be a function of what
c yesterday was (wet or dry)
c
c the logic here is that it is more likely that a wet day will follow
c another wet day -- allowing for 'storm events' to persist
c
c calculate monthly-average probability of rainy day 
c
          pwet = max (1., xinwet(i,imonth)) / ndaypm(imonth)
c
c estimate the probability of a wet day after a dry day
c
          pwd = 0.75 * pwet
c
c estimate the probability of a wet day after a wet day
c
          pww = 0.25 + pwd
c
c Beginning of block of code that figures out daily precip for
c entire month on the first day of each month
c
          if (iday .eq. 1) then
c
c Verify the dataset consistency especially when using interannual anomalies
c of precipitations (negative values, or too few wet days in a rainy month)
c
            xinprec(i, imonth) = max(0., xinprec(i, imonth))
            xinwet(i, imonth) = max(1., xinwet(i, imonth))
c
 9000       continue
c
c Initialize monthly parameters back to zero
c
            iwetdaysum(i) = 0
            precipdaysum(i) = 0
c
            do 210 j = 1, 31
              iwetday(i,j) = 0
              precipday(i,j) = 0
 210        continue
c
c Loop through number of days in this month and determine precip
c
            do 220 j = 1, ndaypm(imonth)
c
c decide if today is a wet day or a dry day using a random number
c
              rndnum = ran2(seed)
c
c
c If it is the first day of the month do not look at previous day
c
              if (j .eq. 1) then
                if (rndnum .le. pwd) then
                  iwetday(i,j) = 1
                  iwetdaysum(i) = iwetdaysum(i) + 1
                else
                  iwetday(i,j) = 0
                endif
c
              else
c
c If it is not the first day, look at yesterday's wet/dry index to help
c determine if today is wet or dry
c
                if (iwetday(i,j-1) .eq. 0) then
                  if (rndnum.le.pwd) then
                    iwetday(i,j) = 1
                    iwetdaysum(i) = iwetdaysum(i) + 1
                  endif
                else
                  if (rndnum.gt.pww) iwetday(i,j) = 0
                endif
              endif
c
c ---------------------------------------------------------------------- 
c (2) determine today's precipitation amount (following Geng et al.)
c ---------------------------------------------------------------------- 
c
c if it is going to rain today
c
              if (iwetday(i,j) .eq. 1) then
c
c calculate average rainfall per wet day
c
                rainpwd = xinprec(i,imonth) * ndaypm(imonth) /
     >                max (0.1, xinwet(i,imonth))
c
c randomly select a daily rainfall amount from a probability density
c function of rainfall
c
c method i --
c
c use the following technique from Geng et al. and Richardson
c to distribute rainfall probabilities
c 
c pick a random rainfall amount from a two-parameter gamma function
c distribution function
c
c estimate two parameters for gamma function (following Geng et al.)
c 
                beta  = max (1.0, -2.16 + 1.83 * rainpwd)
                alpha = rainpwd / beta
c
c determine daily precipitation amount from gamma distribution function
c (following WGEN code of Richardson and Wright (1984))
c
                aa = 1.0 / alpha 
                ab = 1.0 / (1.0 - alpha)
c
                tr1 = exp(-18.42 / aa)
                tr2 = exp(-18.42 / ab)
c
 12             rn1 = ran2(seed)
                rn2 = ran2(seed)
c
c CD: rewrote parts of prehistoric code in fortran 77 
c
                if ((rn1 - tr1) .le. 0) then
                  s1 = 0.0
                else 
                  s1 = rn1**aa
                end if
c
                if ((rn2 - tr2) .le. 0) then 
                  s2 = 0.0
                else 
                  s2 = rn2**ab
                end if
           
c               if (rn1 - tr1) 61, 61, 62
c 61            s1 = 0.0
c               go to 63
c 62            s1 = rn1**aa
c
c 63            if (rn2 - tr2) 64, 64, 65
c 64            s2 = 0.0
c               go to 66
c 65            s2 = rn2**ab
c
c
 66             s12 = s1 + s2
c
                if (s12 - 1.0)  13, 13, 12
 13             z = s1 / s12
c
                rn3 = ran2(seed)
c
                precipday(i,j) = -z * log(rn3) * beta
c
c method ii --
c
c here we use a one-parameter Weibull distribution function
c following the analysis of Selker and Haith (1990)
c
c Selker, J.S. and D.A. Haith, 1990: Development and testing of single-
c parameter precipitation distributions, Water Resources Research, 
c 11, 2733-2740.
c
c this technique seems to have a significant advantage over other
c means of generating rainfall distribution functions
c
c by calibrating the Weibull function to U.S. precipitation records,
c Selker and Haith were able to establish the following relationship
c
c the cumulative probability of rainfall intensity x is given as:
c
c f(x) = 1.0 - exp(-(1.191 x / rainpwd)**0.75)
c 
c where x       : rainfall intensity
c       rainpwd : rainfall per wet day
c
c using transformation method, take uniform deviate and convert it to a
c random number weighted by the following Weibull function
c
c          rndnum = ran2(seed)
c
c          precip(i) = rainpwd / 1.191 * (-log(1.0 - rndnum))**1.333333
c 
c bound daily precipitation to "realistic" range
c 
c lower end is determined by definition of a 'wet day' (at least
c 0.25 mm of total precipitation)
c
c upper end is to prevent ibis from blowing up
c
                precipday(i,j) = max (precipday(i,j),0.25) ! min =   0.25 mm/day
                precipday(i,j) = min (precipday(i,j),150.00) ! max = 150.00 mm/day
c
c Back to beginning of month loop, this is the end of it
c
              endif
c
c Add today's precip to the monthly summation
c
              precipdaysum(i) = precipdaysum(i) + precipday(i,j)
c
 220        continue
c
c Adjust daily precip amounts (using precipfac) so that the monthly
c summation equals the input precip amount, when using interannual
c anomalies
c
c
              if ((precipdaysum(i) .eq. 0) .AND.
     >                (xinprec(i,imonth) .gt. 0.)) then
                    rndnum = 1.0 + (float(ndaypm(imonth)) - 1.0)
     >                   * ran2(seed)
                    iwetday(i,nint(rndnum)) = 1
                    precipday(i,nint(rndnum)) = xinprec(i,imonth)
     >                   * float(ndaypm(imonth))
                    precipdaysum(i) = precipday(i,nint(rndnum))
                    iwetdaysum(i) = 1
                 end if
c
                 precipfac = (xinprec(i,imonth)*float(ndaypm(imonth)))
     >                / max(0.01,precipdaysum(i))
c
                 do 230 j=1,ndaypm(imonth)
                    precipday(i,j) = precipday(i,j) * precipfac
c
                    if (precipday(i,j).gt.360) then
                       if (xinwet(i,imonth) .lt. ndaypm(imonth)) then
                          xinwet(i,imonth) = xinwet(i,imonth) + 1
                          pwet = xinwet(i,imonth) / ndaypm(imonth)
                          pwd = 0.75 * pwet
                          pww = 0.25 + pwd
                          print *, 'WARNING: goto 9000a', i, int(xinwet(i
     >                         ,imonth)), iwetdaysum(i), int(precipday(i,j))
                          goto 9000
                       else
                          print *, 'WARNING: goto 9000b', i, int(xinwet(i
     >                         ,imonth)), iwetdaysum(i), int(precipday(i,j))
                          goto 9000
                       end if
                    end if
c    
  230            continue
c    
c Verification of the weather generator algorithm
c
                 iwetdaysum(i) = 0
                 precipdaysum(i) = 0.
c
                 do 240 j=1,ndaypm(imonth)
                    precipdaysum(i) = precipdaysum(i) + precipday(i,j)
                    iwetdaysum(i) = iwetdaysum(i) + iwetday(i,j)
  240            continue
c
                 dif = precipdaysum(i) - xinprec(i,imonth)
     >                   * float(ndaypm(imonth))
c    
                 if ((dif.lt.-0.1).or.(dif.gt.0.1)) print *,
     >                'ERROR in DAILY:', i, precipdaysum(i),
     >                xinprec(i,imonth)* float(ndaypm(imonth)),
     >                iwetdaysum(i), xinwet(i,imonth)
c    
c end of the verification
c
           end if               !end of the iday loop
c    
c Relate today's iwetday and precipday to iwet and precip that will be
c used below
           iwet(i) = iwetday(i,iday)
           precip(i) = precipday(i,iday)
c
c ---------------------------------------------------------------------- 
c (3) estimate expected minimum and maximum temperatures
c ---------------------------------------------------------------------- 
c
c first determine the expected maximum and minimum temperatures
c (from climatological means) for this day of the year
c
c mean daily mean temperature (K)
c
          tdm = xint(i,it1w) +
     >          dt * (xint(i,it2w) - xint(i,it1w)) + 273.16
c
c mean daily temperature range (K)
c
          trngm = xintrng(i,it1w) +
     >            dt * (xintrng(i,it2w) - xintrng(i,it1w))
c
c mean minimum and maximum temperatures
c
          tmaxm = tdm + 0.56 * trngm
          tminm = tdm - 0.44 * trngm
c
c modify maximum temperatures for wet and dry days
c
          if (pwet .ne. 0.0) then
            tmaxd = tmaxm + pwet * omtmax * trngm
            tmaxw = tmaxd -        omtmax * trngm
          else
            tmaxd = tmaxm
            tmaxw = tmaxm
          endif
c
c set the 'expected' maximum and minimum temperatures for today
c
c note that the expected minimum temperatures are the same for
c both wet and dry days
c
          if (iwet(i).eq.0) tmaxe = tmaxd
          if (iwet(i).eq.1) tmaxe = tmaxw
c
          tmine = tminm
c
c estimate variability in minimum and maximum temperatures
c
c tmaxs : standard deviation in maximum temperature (K)
c tmins : standard deviation in minimum temperature (K)
c
c Regression is based on analysis of 2-m air temperature data from the
c NCEP/NCAR reanalysis (1958-1997) for 294 land points over central
c North America (24N-52N, 130W-60W, 0.5-degree resolution): Daily max
c and min temperatures were calculated for each land point from daily
c mean temperature and temperature range. Anomalies were calculated
c by subtracting similar max and min temperatures calculated from
c monthly mean temperature and range (interpolated to daily). The 40
c years of anomalies were then binned by month and the standard
c deviation calculated for each month. The 294 x 12 standard
c deviations were then regressed against the 3528 long-term monthly
c mean temperatures.
c
c note: the values are bound to be greater than 1.0 K 
c       (at the very least they must be bound so they don't go below zero)
c
          tmaxs = max (1.0, -0.0713 * (tdm - 273.16) + 4.89)
          tmins = max (1.0, -0.1280 * (tdm - 273.16) + 5.73)
c
c ---------------------------------------------------------------------- 
c (4) estimate expected cloud cover
c ---------------------------------------------------------------------- 
c
c the formulation of dry and wet cloud cover has been
c derived from the weather generator used in the epic crop model
c
c cloudm : mean cloud cover for today
c cloudd : dry day cloud cover
c cloudw : wet day cloud cover
c cloude : expected cloud cover today
c
c Verify the data set consistency when using interannual anomalies of
c cloudiness (values under 0 % or over 100 %)
c
          if (iday.eq.1) then
             xincld(i,it1w) = max (0., xincld(i,it1w))
             xincld(i,it1w) = min (100., xincld(i,it1w))
             xincld(i,it2w) = max (0., xincld(i,it2w))
             xincld(i,it2w) = min (100., xincld(i,it2w))
          end if

c monthly mean cloud cover (%)
c
          cloudm = xincld(i,it1w) +
     >             dt * (xincld(i,it2w) - xincld(i,it1w))
c 
c convert from percent to fraction
c
          cloudm = cloudm / 100.0
c
c adjust cloud cover depending on dry day / rainy day
c following logic of the EPIC weather generator code
c
          if (pwet .ne. 0.0) then
            cloudd = (cloudm - pwet * omcloud) / (1.0 - pwet * omcloud)
            cloudd = min (1.0, max (0.0, cloudd))
            cloudw = (cloudm - (1.0 - pwet) * cloudd) / pwet
          else
            cloudd = cloudm
            cloudw = cloudm
          endif
c
          if (iwet(i).eq.0) cloude = cloudd
          if (iwet(i).eq.1) cloude = cloudw
c
c estimate variability in cloud cover for wet and dry days
c following numbers proposed by Richardson
c
c clouds : standard deviation of cloud fraction
c 
          if (iwet(i).eq.0) clouds = 0.24 * cloude
          if (iwet(i).eq.1) clouds = 0.48 * cloude
c
c ---------------------------------------------------------------------- 
c (5) determine today's temperatures and cloud cover using
c     first-order serial autoregressive technique
c ---------------------------------------------------------------------- 
c
c use the Richardson (1981) weather generator approach to simulate the
c daily values of minimum / maximum temperature and cloud cover
c
c following the implementation of the Richardson WGEN weather generator
c used in the EPIC crop model
c
c this approach uses a multivariate generator, which assumes that the
c perturbation of minimum / maximum temperature and cloud cover are
c normally distributed and that the serial correlation of each
c variable may be described by a first-order autoregressive model
c
c generate standard deviates for weather generator
c
          do j = 1, 3
 31         rn1 = ran2(seed)
            rn2 = ran2(seed)
            v = sqrt (-2.0 * log(rn1)) * cos(6.283185 * rn2)
            if (abs(v) .gt. 2.5) go to 31
            ee(j) = v
          enddo
c
c zero out vectors
c
          do j = 1, 3
            r(j)  = 0.0
            rr(j) = 0.0
          enddo
c
c update working vectors
c
          do j = 1, 3
            do k = 1, 3
              r(j)  = r(j)  + b(j,k) * ee(j)
              rr(j) = rr(j) + a(j,k) * xstore(i,k)
            enddo
          enddo
c
c solve for x() perturbation vector and save current vector
c into the xim1() storage vector (saved for each point)
c 
          do j = 1, 3
            x(j) = r(j) + rr(j)
            xstore(i,j) = x(j)
          enddo
c
c determine today's minimum and maximum temperature
c
          tmax(i)  = tmaxe + tmaxs * x(1)
          tmin(i)  = tmine + tmins * x(2)
c
c if tmin > tmax, then switch the two around
c
          if (tmin(i).gt.tmax(i)) then
            tdum    = tmax(i)
            tmax(i) = tmin(i)
            tmin(i) = tdum
          endif
c
c daily average temperature
c
          td(i) = 0.44 * tmax(i) + 0.56 * tmin(i)
c
c determine today's cloud cover
c
          cloud(i) = cloude + clouds * x(3)
c
c constrain cloud cover to be between 0 and 100%
c
          cloud(i) = max (0.0, min (1.0, cloud(i)))
c
c ---------------------------------------------------------------------- 
c (6) estimate today's surface atmospheric pressure
c ---------------------------------------------------------------------- 
c
c simply a function of the daily average temperature and topographic
c height -- nothing fancy here
c
          psurf(i) = 101325.0 *
     >               (td(i) / (td(i) + 0.0065 * xintopo(i))) ** rwork
c
c ---------------------------------------------------------------------- 
c (7) estimate today's relative humidity
c ---------------------------------------------------------------------- 
c
c the formulation of dry and wet relative humidities has been
c derived from the weather generator used in the epic crop model
c
c qdm : mean relative humidity 
c qdd : dry day relative humidity
c qdw : rainy day relative humidity
c qde : expected relative humidity (based on wet/dry decision)
c
c Verify the data set consistency when using interannual anomalies of
c relative humidity (values over 100 % or under 0 %)
c
          if (iday.eq.1) then
             xinq(i,it1w) = max (0., xinq(i,it1w))
             xinq(i,it1w) = min (100., xinq(i,it1w))
             xinq(i,it2w) = max (0., xinq(i,it2w))
             xinq(i,it2w) = min (100., xinq(i,it2w))
          end if
c
c mean relative humidity (%)
c
          qdm = xinq(i,it1w) + dt * (xinq(i,it2w) - xinq(i,it1w))
c 
c convert from percent to fraction
c
          qdm = qdm / 100.0
c
c adjust humidity depending on dry day / rainy day
c following logic of the EPIC weather generator code
c
          if (pwet .ne. 0.0) then
            qdd = (qdm - pwet * omqd) / (1.0 - pwet * omqd)
            if (qdd .lt. 0.2) then
              qdd = 0.2
              if (qdd .gt. qdm) qdm = qdd
            endif
            qdd = min(1.0, qdd)
            qdw = (qdm - (1.0 - pwet) * qdd) / pwet
          else
            qdd = qdm
            qdw = qdm
          endif
c
          if (iwet(i).eq.0) qde = qdd
          if (iwet(i).eq.1) qde = qdw 
c
c estimate lower and upper bounds of humidity distribution function
c following logic of the EPIC weather generator code
c
          qdup  = qde + (1.0 - qde) * exp (qde - 1.0)
          qdlow = qde * (1.0 - exp (-qde))
c
c randomly select humidity from triangular distribution function
c following logic of the EPIC weather generator code
c
          rn = ran2(seed)
c
          y  = 2.0 / (qdup - qdlow)
c
          b3 = qde  - qdlow
          b2 = qdup - qde
          b1 = rn / y
c
          x1 = y * b3 / 2.0 
c
          if (rn.gt.x1) then
            qd(i) = qdup  - sqrt (b2 * b2 - 2.0 * b2 * (b1 - 0.5 * b3))
          else
            qd(i) = qdlow + sqrt (2.0 * b1 * b3)
          endif
c
c adjust daily humidity to conserve monthly mean values
c
c note that this adjustment sometimes gives rise to humidity
c values greater than 1.0 -- which is corrected below
c
          amn = (qdup + qde + qdlow) / 3.0
          qd(i) = qd(i) * qde / amn
c
c constrain daily average relative humidity
c
          qd(i) = max (0.30, qd(i))
          qd(i) = min (0.99, qd(i))
c
c convert from relative humidity to specific humidity at
c daily mean temperature
c
          qd(i) = qd(i) * qsat(esat(td(i)), psurf(i))
c
c ---------------------------------------------------------------------- 
c (8) estimate today's daily average wind speed
c ---------------------------------------------------------------------- 
c
c first estimate the expected daily average wind speed (from monthly
c means)
c
          eud = xinwind(i,it1w) +
     >          dt * (xinwind(i,it2w) - xinwind(i,it1w))
c
c following logic of the EPIC weather generator
c select random wind speed following this equation
c
          ud(i) = 1.13989 * eud * (-log(ran2(seed)))**0.30 
c
c constrain daily wind speeds to be between 2.5 and 10.0 m/sec
c
          ud(i) = max (2.5, min (10.0, ud(i)))
c
c ---------------------------------------------------------------------- 
c * * * use real daily climate data * * *
c ---------------------------------------------------------------------- 
c
        else
c
c use basic daily climate data, converting units
c
c daily total precipitation
c
c Here we multiply xinprecd, the daily fraction of precip calculated from
c the NCEP dataset, by xinprec, the total monthly amount of precip taken from
c the CRU05 dataset to obtain our derived daily precip amount. Also do a check
c to see if the daily precip exceeds 360mm (as is done in the daily weather 
c generator) ... no correction is made, only a warning is printed
c
          precip(i) = (xinprec(i,imonth) * ndaypm(imonth)) * xinprecd(i)
c          if (precip(i) .gt. 360) then
c            print *, 'WARNING: daily precip exceeds 360mm for'
c            print *, 'year, month, day, gridcell = '
c            print *, iyear, imonth, iday, i
c          endif
c
c daily average temperatures
c
c Here we add the NCEP temperature anomaly to the CRU05 monthly anomaly
c The trange NCEP anomaly must also be multiplied by the climatological
c CRU05 trange in order to determine tmax and tmin
c
          td(i) = xint(i,imonth) + 273.16 + xintd(i)
          trngm = min (44.0, (xintrng(i,imonth) * xintrngd(i)))
c
          tmax(i) = td(i) + 0.56 * trngm
          tmin(i) = td(i) - 0.44 * trngm
c
c daily average cloud cover
c
c Here we add the NCEP cloud anomaly to the monthly anomaly from CRU05
c before converting percentage of cover to fraction
c We also bound cloud cover fraction between 0 and 1
c
          cloud(i) = (xincld(i,imonth) + xincldd(i)) * 0.01
c
          cloud(i) = min (cloud(i), 1.)
          cloud(i) = max (0.0, cloud(i))
c
c compute surface atmospheric pressure
c
          psurf(i) = 101325.0 *
     >               (td(i) / (td(i) + 0.0065 * xintopo(i))) ** rwork
c
c daily average specific humidity
c
c First we must convert relative humidity to a fraction and then convert
c the fraction to specific humidity
c Then we can multiply the NCEP daily anomaly by the CRU05 monthly
c anomaly
c
          humidfrac = xinq(i,imonth) / 100.
          sphumid = humidfrac * qsat(esat(td(i)),psurf(i))
c
          qd(i) = sphumid * xinqd(i)
c
c daily average wind speed
c
c Here we multiply the NCEP fraction of windspeed by the CRU05
c climatological monthly windspeed
c
          ud(i) = xinwind(i,imonth) * xinwindd(i)
c
c
        end if
c
c ---------------------------------------------------------------------- 
c * * * other daily climate calculations * * *
c ---------------------------------------------------------------------- 
c
c calculated temperature extremes -- for vegetation limits (deg c)
c
c for this purpose, use the 10-day running mean temperature
c
        tcthis(i) = min (tcthis(i), (a10td(i) - 273.16))
        twthis(i) = max (twthis(i), (a10td(i) - 273.16))
c
c update this year's growing degree days
c
        gdd0this(i) = gdd0this(i) + max(0., (td(i) - 273.16))
        gdd5this(i) = gdd5this(i) + max(0., (td(i) - 278.16))
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
      subroutine diurnal (time, jday, plens, startp, endp, seed)
c ---------------------------------------------------------------------
c
c common blocks
c
      use implicit
c
      use compar
      use comatm
      use comwork
c
c Arguments
c
      integer jday,    ! current day
     >        seed     !
c
      real time,
     >     plens,      ! length of the precip event (s)
     >     startp,     ! time to start precip event (s)
     >     endp        ! time to end precip event (s)
c
c determine the length of a precipitation event (between 4 and 24 hours),
c and time to start and end precipitation event. plen is in timesteps, while
c plens, startp, and endp are in seconds
c
c
c local variables
c
      integer i,         ! loop indice
     >        jj, 
     >        ib         ! waveband number 1= visible, 2= near-IR
c
      real rtime,        ! time in hours
     >     orbit,        ! earth's orbital angle (around the sun) in radians
     >     angle,        ! solar hour angle in radians
     >     xdecl,        ! solar declination angle
     >     sw,           ! effective solar constant
     >     xlat,         ! latitude in radians
     >     trans,        ! solar transmission through the atmosphere
     >     fdiffuse,     ! fraction of indirect (diffuse) solar radiation
     >     frac,         ! fraction of energy in each waveband
     >     gamma,        !
     >     qmin,         !
     >     qmax,
     >     qsa,
     >     ran2,
     >     emb,
     >     ea,
     >     ec,
     >     dtair,
     >     dtcloud
c
      integer checkP,
     >        niter,
     >        plen,
     >        plenmin,
     >        plenmax
c
      use comsat
c

c ---------------------------------------------------------------------- 
c * * * calendar and orbital calculations * * *
c ---------------------------------------------------------------------- 
c
c calculate time in hours
c
      rtime = time / 3600.0
c
c calculate the earth's orbital angle (around the sun) in radians
c
      orbit = 2.0 * pi * float(jday) / 365.2425
c
c calculate the solar hour angle in radians
c
      angle  = 2.0 * pi * (rtime - 12.0) / 24.0
c
c calculate the current solar declination angle
c ref: global physical climatology, hartmann, appendix a
c
      xdecl =  0.006918                    -
     >         0.399912 * cos(orbit)       + 
     >         0.070257 * sin(orbit)       -
     >         0.006758 * cos(2.0 * orbit) +
     >         0.000907 * sin(2.0 * orbit) -
     >         0.002697 * cos(3.0 * orbit) +
     >         0.001480 * sin(3.0 * orbit)
c
c calculate the effective solar constant, including effects of eccentricity
c ref: global physical climatology, hartmann, appendix a
c
      sw = 1370. * (1.000110 +
     >              0.034221 * cos(orbit)        + 
     >              0.001280 * sin(orbit)        + 
     >              0.000719 * cos(2.0 * orbit)  +
     >              0.000077 * sin(2.0 * orbit)) 
c
 9001 continue
c
c do for all gridcells
c
      do 100 i = 1, npoi
c
c ---------------------------------------------------------------------- 
c * * * solar calculations * * *
c ---------------------------------------------------------------------- 
c
c calculate the latitude in radians
c
        jj = latindex(i)
c
        xlat = latscale(jj) * pi / 180.0
c
c calculate the cosine of the solar zenith angle
c
        coszen(i) = max (0.0, (sin(xlat) * sin(xdecl) +
     >                         cos(xlat) * cos(xdecl) * cos(angle)))
c
c calculate the solar transmission through the atmosphere
c using simple linear function of tranmission and cloud cover
c
c note that the 'cloud cover' data is typically obtained from
c sunshine hours -- not direct cloud observations
c
c where, cloud cover = 1 - sunshine fraction 
c
c different authors present different values for the slope and 
c intercept terms of this equation
c
c Friend, A: Parameterization of a global daily weather generator for
c terrestrial ecosystem and biogeochemical modelling, Ecological 
c Modelling
c
c Spitters et al., 1986: Separating the diffuse and direct component
c of global radiation and its implications for modeling canopy
c photosynthesis, Part I: Components of incoming radiation,
c Agricultural and Forest Meteorology, 38, 217-229.
c
c A. Friend       : trans = 0.251 + 0.509 * (1.0 - cloud(i))
c Spitters et al. : trans = 0.200 + 0.560 * (1.0 - cloud(i))
c
c we are using the values from A. Friend
c
        trans = 0.251 + 0.509 * (1.0 - cloud(i)) 
c
c calculate the fraction of indirect (diffuse) solar radiation
c based upon the cloud cover
c
c note that these relationships typically are measured for either
c monthly or daily timescales, and may not be exactly appropriate
c for hourly calculations -- however, in ibis, cloud cover is fixed
c through the entire day so it may not make much difference
c
c method i --
c
c we use a simple empirical relationships from Nikolov and Zeller (1992)
c
c Nikolov, N. and K.F. Zeller, 1992:  A solar radiation algorithm for ecosystem
c dynamics models, Ecological Modelling, 61, 149-168.
c
        fdiffuse = 1.0045 + 0.0435 * trans 
     >                    - 3.5227 * trans**2
     >                    + 2.6313 * trans**3
c
        if (trans.gt.0.75) fdiffuse = 0.166
c
c method ii --
c
c another method was suggested by Spitters et al. (1986) based on
c long-term data from the Netherlands
c
c Spitters et al., 1986: Separating the diffuse and direct component
c of global radiation and its implications for modeling canopy
c photosynthesis, Part I: Components of incoming radiation,
c Agricultural and Forest Meteorology, 38, 217-229.
c
c       if ((trans.eq.0.00).and.(trans.lt.0.07)) then
c         fdiffuse = 1.0
c       else if ((trans.ge.0.07).and.(trans.lt.0.35)) then
c         fdiffuse = 1.0 - 2.3 * (trans - 0.07)**2
c       else if ((trans.ge.0.35).and.(trans.lt.0.75)) then
c         fdiffuse = 1.33 - 1.46 * trans
c       else
c         fdiffuse = 0.23
c       endif
c
c do for each waveband
c
        do 120 ib = 1, nband
c
c calculate the fraction in each waveband
c
          frac = 0.46 + 0.08 * float(ib - 1)
c
c calculate the direct and indirect solar radiation
c
          solad(i,ib) = sw * coszen(i) * frac * trans *
     >                  (1. - fdiffuse)  
c
          solai(i,ib) = sw * coszen(i) * frac * trans * fdiffuse
c
  120   continue
c
c ---------------------------------------------------------------------- 
c * * * temperature calculations * * *
c ---------------------------------------------------------------------- 
c
c assign hourly temperatures using tmax and tmin 
c following Environmental Biophysics, by Campbell and Norman, p.23
c
c this function fits a fourier series to the diurnal temperature cycle
c note that the maximum temperature occurs at 2:00 pm local solar time
c
c note that the daily mean value of gamma is 0.44, 
c so td = 0.44 * tmax + 0.56 * tmin,  instead of
c    td = 0.50 * tmax + 0.50 * tmin
c
        gamma = 0.44 - 0.46 * sin (      pi / 12.0 * rtime + 0.9) 
     >               + 0.11 * sin (2.0 * pi / 12.0 * rtime + 0.9)
c
        ta(i) = tmax(i) * gamma + tmin(i) * (1.0 - gamma)
c
c ---------------------------------------------------------------------- 
c * * * humidity calculations * * *
c ---------------------------------------------------------------------- 
c
c adjust specific humidity against daily minimum temperatures
c
c To do this, qa is written as an approximate sine function (same as ta)
c to preserve the daily mean specific humidity, while also preventing rh
c from exceeding 99% at night
c
c Note that the daily mean RH is *not* preserved, and therefore the
c output RH will be slightly different from what was read in.
c
c first adjust so that maximum RH cannot exceed 99% at night
c
        qmin = min (qd(i), 0.99 * qsat(esat(tmin(i)), psurf(i)))
        qmax = (qd(i) - 0.56 * qmin) / 0.44
c
c if needed, adjust again to 99% at other times of the day (in which
c case the daily mean *specific* humidity is also not preserved)
c
        qsa  = 0.99 * qsat(esat(ta(i)), psurf(i))
c
c calculate the hourly specific humidity, using the above adjustments
c
        qa(i) = min (qsa, qmax * gamma + qmin * (1.0 - gamma))
c
c calculate the hourly relative humidity 
c
        rh(i) = 100.0 * qa(i) / qsat(esat(ta(i)), psurf(i))
c
c ---------------------------------------------------------------------- 
c * * * wind speed calculations * * *
c ---------------------------------------------------------------------- 
c
c following logic of the EPIC weather generator
c select random wind speed following this equation
c
        ua(i) = 1.13989 * ud(i) * (-log(ran2(seed)))**0.30 
c
c fix wind speeds to always be above 2.5 m/sec and below 10.0 m/sec
c
        ua(i) = max (2.5, min (10.0, ua(i)))
c
c ---------------------------------------------------------------------- 
c * * * ir flux calculations * * *
c ---------------------------------------------------------------------- 
c
c clear-sky emissivity as a function of water vapor pressure
c and atmospheric temperature
c
c calculate the ir emissivity of the clear sky
c using equation from idso (1981) water resources res., 17, 295-304
c
        emb = 0.01 * (psurf(i) * qa(i) / (0.622 + qa(i)))
        ea  = 0.70 + 5.95e-5 * emb * exp (1500.0 / ta(i))
c
c assign the ir emissivity of clouds (assume they are ~black in the ir)
c
        ec = 0.950
c
c assign the temperature difference of emissions (air + cloud) from
c the surface air temperature
c
        dtair   = 0.0
        dtcloud = 0.0
c
c total downward ir is equal to the sum of:
c
c (1) clear sky contribution to downward ir radiation flux
c (2) cloud contribution to downward ir radiation flux
c
        fira(i) = (1. -  cloud(i)) * ea * stef * (ta(i) - dtair  )**4 +
     >                   cloud(i)  * ec * stef * (ta(i) - dtcloud)**4
c
c ---------------------------------------------------------------------- 
c * * * snow and rain calculations * * *
c ---------------------------------------------------------------------- 
c
c reset snow and rain to zero
c
        snowa(i) = 0.0
        raina(i) = 0.0
c
c determine the number of timesteps per day
c
        niter = int (86400.0 / dtime)
c
c change the rain length when the amount of rainfall/timestep is
C too high (at the first time step)
c
        if (time.lt.dtime) then
c
           plen = plens / dtime
           plenmin = 1 +  int ((4.0 * 3600. - 1.) / dtime)
           plenmax = max (int (24.0 * 3600. / dtime), plenmin)
           checkP = 0
c
           do  while (((precip(i)/plen) .gt. 15).and.(plen.lt.plenmax))
              plen = plen + 1
              checkP = 1
           end do
c
           if (checkP.eq.1) then
c
c              print *, 'WARNING: plen changed', i,
c     $             int(precip(i)), int(plens/dtime), plen
              plens = dtime * plen
              startp = dtime * min (niter-plen,
     >             int(ran2(seed)*(niter-plen+1)))
              endp = startp + plen *dtime
              goto 9001
           end if
c
        end if
c
c if precipitation event then calculate
c
        if (time.ge.startp .and. time.lt.endp) then  
c
c for rain / snow partitioning, make it all rain if 
c ta > 2.5 C and all snow if ta <= 2.5 C
c
c reference:
c
c Auer, A. H., 1974: The rain versus snow threshold temperatures,
c Weatherwise, 27, 67.
c
          if (ta(i)-273.15 .gt. 2.5) then
            raina(i) = precip(i) / plens
          else
            snowa(i) = precip(i) / plens
          endif
c
        endif
c
  100 continue
c
      return
      end
c
