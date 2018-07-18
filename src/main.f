c
c #    #    ##       #    #    #
c ##  ##   #  #      #    ##   #
c # ## #  #    #     #    # #  #
c #    #  ######     #    #  # #
c #    #  #    #     #    #   ##
c #    #  #    #     #    #    #
c
c ---------------------------------------------------------------
      program main
c ---------------------------------------------------------------
c
c common blocks
c 
      include 'implicit.h'
c
      include 'compar.h'
      include 'comatm.h'
      include 'comdiag.h'
      include 'argvs.h'
c
c local variables
c

      integer icount,        ! number of times ran2 is called (for restart)
     >        icountdum      !
c
      integer istep,         ! timestep counter (per day)
     >        iday,          ! daily loop counter
     >        imonth,        ! monthly loop counter
     >        istyrd,        ! first yr daily means exist (from precip file)
     >        istyrm,        ! first yr monthly anoms exist (from file)
     >        iwest,         ! 1st lon index for subset
     >        iy1,           ! first year for year loop
     >        iy2,           ! last year for year loop
     >        iyear,         ! yearly loop counter
     >        iyear0,        ! first year of simulation
     >        iyranom,       ! year to start reading monthly anomalies
     >        iyrdaily,      ! year to start reading daily means
     >        iyrlast,       ! last year of previous run (for restart)
     >        idiag,         ! number of diagnostic files requested
     >        idailyout,     ! 0: no daily output 1: daily output
     >        imonthout,     ! 0: no monthly output 1: monthly output
     >        iyearout,      ! 0: no yearly output 1: yearly output
     >        irestart,      ! 0: normal mode 1: restart mode
     >        isimveg,       ! 0: static veg 1: dynam veg initialized w/ fixed
     >                       ! 2: dynam veg initialized w/ cold start
     >        isimfire,      ! 0: fixed fire  1: dynamic fire
     >        isimco2,       ! 0: fixed co2   1: changing co2
     >        jday,          ! julian day of the simulation
     >        jnorth,        ! 1st lat index for subset
     >        nanom,         ! # of years in the anomaly files
     >        niter,         ! total number of time iterations per day
     >        nday,          ! number of days since the start of the simulation
     >        nrun,          ! # of years in this run
     >        nspinsoil,     ! year of simulation that soil c reaches equilibrium 
     >        plen,          ! length of precipitation event in timesteps (see plens)
     >        plenmax,       ! upper bound on plen
     >        plenmin,       ! lower bound on plen
     >        seed,          ! value used to initialize the random # generator
     >        spin,          ! counter for iterations used to spin up soilbgc
     >        spinmax,       ! maximum iterations used to spin up soilbgc
     >        eqyears,       !
     >        soilcspin,     ! 0: no spinup procedure for soil c  1: acceleration procedure used
     >        lun, ifile,    ! file indices
     >        i, n, ivar     ! loop indices
c
      real    co2init,       ! initial co2 concentration in mol/mol
     >        o2init,        ! initial o2 concentration in mol/mol
     >        dran,          ! random number from generator
     >        plens,         ! length of precipitation event in seconds
     >        startp,        ! time to start precipitation event (seconds since midnight)
     >        endp,          ! time to end precipitation event (seconds since midnight)
     >        spincons,      ! # times soil gbc called each day during spinup
     >        spinfrac,      ! fraction of nspinsoil time with max iteration spinmax
     >        slope,         ! rate of decrease of number of iteration for soil C spinup 
     >        snorth,        ! north latitude for subsetting std grid
     >        ssouth,        ! south latitude for subsetting std grid
     >        swest,         ! west longitude for subsetting std grid
     >        seast,         ! east longitude for subsetting std grid
     >        test,          ! test on timestep 
     >        time           ! time of day since midnight (in seconds)
c
c External
c
      real ran2              ! Function : random number generator
c
c number of days per month
c
      data ndaypm /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
      data snorth, ssouth, swest, seast / 90., -90., -180., 180. /
c
c ---------------------------------------------------------------
c                       j o b   c o n t r o l 
c ---------------------------------------------------------------
c
c open local input file called 'ibis.infile'
c and read in pertinent information
c
      character dirname*100
      integer status

c
c tell user about the simulation
c
      write (*,*) ' '
      write (*,*) '****************************************'
      write (*,*) '* IBIS: Integrated BIosphere Simulator *'
      write (*,*) '* Version 2.6b3                        *'
      write (*,*) '* March 2002                           *'
      write (*,*) '****************************************'
      write (*,*) ' '


      write (*,*) ' '
      write (*,*) '****************************************'
      write (*,*) '              parsing argc'
      call parseArgv
      write(*,*) 'ini_infile = ', ini_infile_
      write(*,*) 'diag_infile = ', diag_infile_
      write(*,*) 'cld_mon = ', cld_mon_
      write(*,*) 'deltat_mon = ', deltat_mon_
      write(*,*) 'prec_mon = ', prec_mon_
      write(*,*) 'rh_mon = ', rh_mon_
      write(*,*) 'temp_mon = ', temp_mon_
      write(*,*) 'trange_mon = ', trange_mon_
      write(*,*) 'wetd_mon = ', wetd_mon_
      write(*,*) 'wspd_mon = ', wspd_mon_
      write(*,*) 'soita_sand = ', soita_sand_
      write(*,*) 'soita_clay = ', soita_clay_
      write(*,*) 'vegtype = ', vegtype_
      write(*,*) 'surta = ', surta_
      write(*,*) 'topo = ', topo_
      write(*,*) 'params_can = ', params_can_
      write(*,*) 'params_hyd = ', params_hyd_
      write(*,*) 'params_soi = ', params_soi_
      write(*,*) 'params_veg = ', params_veg_
      write(*,*) 'out_yearly_aet = ', out_yearly_aet_
      write(*,*) 'out_yearly_biomass = ', out_yearly_biomass_
      write(*,*) 'out_yearly_co2fluxes = ', out_yearly_co2fluxes_
      write(*,*) 'out_yearly_csoi = ', out_yearly_csoi_
      write(*,*) 'out_yearly_disturbf = ', out_yearly_disturbf_
      write(*,*) 'out_yearly_exist = ', out_yearly_exist_
      write(*,*) 'out_yearly_fcover = ', out_yearly_fcover_
      write(*,*) 'out_yearly_npp = ', out_yearly_npp_
      write(*,*) 'out_yearly_nsoi = ', out_yearly_nsoi_
      write(*,*) 'out_yearly_plai = ', out_yearly_plai_
      write(*,*) 'out_yearly_runoff = ', out_yearly_runoff_
      write(*,*) 'out_yearly_sens = ', out_yearly_sens_
      write(*,*) 'out_yearly_tsoi = ', out_yearly_tsoi_
      write(*,*) 'out_yearly_vegtype0 = ', out_yearly_vegtype0_
      write(*,*) 'out_yearly_wsoi = ', out_yearly_wsoi_
      write(*,*) 'out_yearly_zcanopy = ', out_yearly_zcanopy_
      write(*,*) 'out_yearly_sapfrac = ', out_yearly_sapfrac_
      write(*,*) 'out_yearly_dummyv = ', out_yearly_dummyv_
      write(*,*) 'out_yearly_solar = ', out_yearly_solar_
      write(*,*) 'out_yearly_albedo = ', out_yearly_albedo_
      write(*,*) 'out_yearly_latent = ', out_yearly_latent_
      write(*,*) 'out_yearly_totfall = ', out_yearly_totfall_
      write(*,*) 'out_yearly_clitw = ', out_yearly_clitw_
      write(*,*) 'out_yearly_csoislo = ', out_yearly_csoislo_
      write(*,*) 'out_yearly_csoipas = ', out_yearly_csoipas_
      write(*,*) 'out_monthly_aet = ', out_monthly_aet_
      write(*,*) 'out_monthly_cloud = ', out_monthly_cloud_
      write(*,*) 'out_monthly_co2ratio = ', out_monthly_co2ratio_
      write(*,*) 'out_monthly_ir = ', out_monthly_ir_
      write(*,*) 'out_monthly_lai = ', out_monthly_lai_
      write(*,*) 'out_monthly_latent = ', out_monthly_latent_
      write(*,*) 'out_monthly_npptot = ', out_monthly_npptot_
      write(*,*) 'out_monthly_qa = ', out_monthly_qa_
      write(*,*) 'out_monthly_rain = ', out_monthly_rain_
      write(*,*) 'out_monthly_rh = ', out_monthly_rh_
      write(*,*) 'out_monthly_runoff = ', out_monthly_runoff_
      write(*,*) 'out_monthly_sens = ', out_monthly_sens_
      write(*,*) 'out_monthly_snod = ', out_monthly_snod_
      write(*,*) 'out_monthly_snof = ', out_monthly_snof_
      write(*,*) 'out_monthly_snow = ', out_monthly_snow_
      write(*,*) 'out_monthly_solar = ', out_monthly_solar_
      write(*,*) 'out_monthly_temp = ', out_monthly_temp_
      write(*,*) 'out_monthly_tsoi = ', out_monthly_tsoi_
      write(*,*) 'out_monthly_wsoi = ', out_monthly_wsoi_
      write(*,*) 'out_monthly_albedo = ', out_monthly_albedo_
      write(*,*) 'out_monthly_dummyv = ', out_monthly_dummyv_
      write(*,*) 'out_daily_rain = ', out_daily_rain_
      write(*,*) 'out_daily_cloud = ', out_daily_cloud_
      write(*,*) 'out_daily_rh = ', out_daily_rh_
      write(*,*) 'out_daily_snow = ', out_daily_snow_
      write(*,*) 'out_daily_aet = ', out_daily_aet_
      write(*,*) 'out_daily_trunoff = ', out_daily_trunoff_
      write(*,*) 'out_daily_srunoff = ', out_daily_srunoff_
      write(*,*) 'out_daily_drainage = ', out_daily_drainage_
      write(*,*) 'out_daily_wsoi = ', out_daily_wsoi_
      write(*,*) 'out_daily_wisoi = ', out_daily_wisoi_
      write(*,*) 'out_daily_snod = ', out_daily_snod_
      write(*,*) 'out_daily_snof = ', out_daily_snof_
      write(*,*) 'out_daily_co2ratio = ', out_daily_co2ratio_
      write(*,*) 'out_daily_co2mic = ', out_daily_co2mic_
      write(*,*) 'out_daily_templ = ', out_daily_templ_
      write(*,*) 'out_daily_zcanopy = ', out_daily_zcanopy_
      write(*,*) 'out_daily_laicanopy = ', out_daily_laicanopy_
      write(*,*) 'out_global = ', out_global_
      write(*,*) 'out_vegtype = ', out_vegtype_
      write(*,*) 'out_yearsrun = ', out_yearsrun_
      write(*,*) 'temp_danom = ', temp_danom_
      write(*,*) 'trange_danom = ', trange_danom_
      write(*,*) 'prec_danom = ', prec_danom_
      write(*,*) 'cld_danom = ', cld_danom_
      write(*,*) 'rh_danom = ', rh_danom_
      write(*,*) 'wspd_danom = ', wspd_danom_
      write(*,*) 'wetd_danom = ', wetd_danom_
      write(*,*) 'prec_daily = ', prec_daily_
      write(*,*) 'temp_daily = ', temp_daily_
      write(*,*) 'trange_daily = ', trange_daily_
      write(*,*) 'cld_daily = ', cld_daily_
      write(*,*) 'wspd_daily = ', wspd_daily_
      write(*,*) 'sphum_daily = ', sphum_daily_
      write(*,*) 'prec_fanom = ', prec_fanom_
      write(*,*) 'trange_fanomc = ', trange_fanomc_
      write(*,*) 'wspd_fanomc = ', wspd_fanomc_
      write(*,*) 'deltat = ', deltat_
      write(*,*) 'sphum_fanom = ', sphum_fanom_
      write (*,*) 'parse argvs finished!'
      write (*,*) '****************************************'
      write (*,*) ' '

      lun = 12
      status = getcwd(dirname)
      write (*,*) dirname, status
      open (lun, status='old', file=ini_infile_)
      read (lun,*) irestart
      read (lun,*) iyear0
      read (lun,*) nrun
      read (lun,*) iyranom
      read (lun,*) nanom
      read (lun,*) iyrdaily
      read (lun,*) soilcspin
      read (lun,*) iyearout
      read (lun,*) imonthout
      read (lun,*) idailyout
      read (lun,*) isimveg
      read (lun,*) isimfire
      read (lun,*) isimco2
      read (lun,*) co2init
      read (lun,*) o2init
      read (lun,*) dtime
      read (lun,*) idiag
      read (lun,*,end=99) snorth
      read (lun,*,end=99) ssouth
      read (lun,*,end=99) swest
      read (lun,*,end=99) seast
c     read (lun,*) nlon
c     read (lun,*) nlat
c     read (lun,*) npoi
c     read (lun,*) xres
c     read (lun,*) yres
c     read (lun,*) nband
c     read (lun,*) nsoilay
c     read (lun,*) nsnolay
c     read (lun,*) npft
 99   close (lun)

      call testCommon(5)
c     allocate(garea(npoi))
c     allocate(vzero(npoi))
c
      if (irestart .eq. 1) then
        write (*,*) 'running in restart mode'
      end if
c
      write (*,*) ' '
      write (*,9000) nrun
      write (*,9005) iyranom
      write (*,*) ' '
      write (*,9010) xres, yres
      write (*,9020) nlon, nlat
      write (*,*) ' '
c
 9000 format (1x,'length of this simulation (years)   : ',i8)
 9005 format (1x,'year to begin using anomalies       : ',i8)
 9006 format (1x,'number of iterations per day        : ',i8)
 9010 format (1x,'model lon, lat resolution (degrees) : ',2f8.2)
 9020 format (1x,'model domain (nlon x nlat)          : ',i3,
     >        ' by ',i3)  
 9030 format (1x,'last year run in this sequence      : ',i8)
c
c ---------------------------------------------------------------
c                     t i m e   c o n t r o l 
c ---------------------------------------------------------------
c
c determine the number of timesteps per day
c
      niter = int (86400.0 / dtime)
c
c test on length of the time step and total number of iteration per day
c
      test = amod (86400.0, dtime)
c
      if (test .gt. 1.e-20) then
        write (*,*) 'dtime ', dtime, ' should be divisible into 86400'  
        stop
      else
        write (*,9006) niter
      end if
c
c check for restart conditions and set "iyrlast" accordingly, open
c ibis.annual in append mode
c
      if (irestart .eq. 1) then
        open (13, status='old', file= out_yearsrun_, err=9050)
        read (13,*) iyrlast
        close (13)
        goto 9059
 9050   write (*,*) 'warning: ibis.infile restart flag=1, indicating'
     >   //' restart, but yearsrun.dat not found'
        write (*,*) 'beginning from ', iyear0
        iyrlast = iyear0 - 1
 9059   open(20,file=out_global_,status='unknown',access='append')
        open(30,file=out_vegtype_,status='unknown',
     >       access='append')
      else
        iyrlast = iyear0 - 1
      end if
c
      write (*,9030) nrun + iyrlast
      write (*,*) ' '
c
c ---------------------------------------------------------------
c                   i n i t i a l i z a t i o n
c ---------------------------------------------------------------
c
c call RD_PARAM to read in all parameter files
c
      call rd_param
c
c initial concentration of co2, o2
c
      co2conc = co2init
      o2conc  = o2init
c
c read global boundary condition datasets
c
      call readit (isimveg,snorth,ssouth,swest,seast,iwest,jnorth)
c
c check if diagnostic output is requested, if so read info from 'diag.infile'
c
      if (idiag .ne. 0) call inidiag(idiag)
c
c preliminary analysis of climate data
c
      if (irestart .eq. 0) then
        call climanl
      end if
c
c initialize monthly anomalies, daily means if needed
c
      if(iyranom .le. iyrlast+nrun) then
c       call inird('input/anom/temp.danom.nc',istyrm)
        call inird(temp_danom_,istyrm)
        if (istyrm .le. 0) istyrm = 9999
      end if
c
      if (iyrdaily .le. iyrlast+nrun) then
c       call inird('input/anom/daily/prec.fanom.nc',istyrd)
        call inird(prec_fanom_,istyrd)
        if (istyrd .le. 0) istyrd = 9999
      end if
c
c if first year of this run/restart is an anomaly year, then read in 
c anomalies for month 12 of previous year and month 1 of this year
c (note: rdanom reads imonth+1) 
c
      iy2 = iyrlast + nrun
      if (iyrlast+1 .ge. iyranom) then
        call rdanom (11, iyrlast, iyranom, nanom, iy2, istyrm, iwest, jnorth)
        call rdanom (12, iyrlast, iyranom, nanom, iy2, istyrm, iwest, jnorth)
      end if
c
c initialize the model
c
      call initial (isimveg, irestart, iyrlast)
c
c initialize random number generator seed
c
      seed = -1
      dran = ran2 (seed)
c
c ---------------------------------------------------------------
c              s t a r t   o f   t i m e   l o o p s
c ---------------------------------------------------------------
c
      write (*,*) ' '
      write (*,*) '********************************'
      write (*,*) '* start of the IBIS simulation *'
      write (*,*) '********************************'
      write (*,*) ' '
c
c reset elapsed number of days, accumulating past days if restart mode
c
      nday = 0
c
      if (irestart .eq. 1) then
        do 150 iyear = iyear0, iyrlast
          nday = nday + 365
          if (mod(iyear,4).eq.0) then
            if (mod(iyear,100).ne.0) then
              nday = nday + 1
            else if (mod(iyear/100,4).eq.0) then
              nday = nday + 1
            end if
          end if
 150    continue
      end if
c
c start of yearly loop
c
      iy1 = iyrlast + 1
      iy2 = iyrlast + nrun
c
      do 200 iyear = iy1, iy2
c
c reset julian date
c
        jday = 0
c
c determine the calendar for this year
c
c leap years occur every calendar year that is divisible by 4 (such as
c 1980, 1984, etc.) *except* for years which begin centuries not divisible
c by 400
c
c for example, 1700, 1800, and 1900 are *not* leap years, but 1600
c and 2000 *are* leap years
c
        ndaypm(2) = 28
        ndaypy = 365
c
        if (mod(iyear,4).eq.0) then
          if (mod(iyear,100).ne.0) then
            ndaypm(2) = 29
            ndaypy = 366
          else if (mod(iyear/100,4).eq.0) then
          ndaypm(2) = 29
          ndaypy = 366
          end if
        end if
c
c start of monthly loop
c 
         do 210 imonth = 1, 12
c
          call rdanom (imonth, iyear, iyranom, nanom, iy2, istyrm, iwest, jnorth)
c
c start of daily loop
c
         do 220 iday = 1, ndaypm(imonth)
c
c update user on model progress once a day
c
            write (*,9100) iyear, imonth, iday
 9100       format (1x,'starting - year/month/day: ',i4,'/',i2,'/',i2)
c
c update elapsed number of days and julian date
c
            nday = nday + 1
            jday = jday + 1
c
c get daily means (if requested)
c and determine today's climatic conditions
c
            if (iyear.ge.iyrdaily) then
              call rdday(jday,imonth,iyear,istyrd,iwest,jnorth)
              call daily (iyear, imonth, iday, seed, 1, iyranom, iyrlast, nrun)
            else
              call daily (iyear, imonth, iday, seed, 0, iyranom, iyrlast, nrun)
            end if
c
c determine the daily vegetation cover characteristics
c
            call pheno
c
c call soil biogeochemistry model
c
c soil spin up procedure calls soil bgc model spincons times
c for each time step to spin up the process (accelerate C & N pools)
c
            spinfrac  = 0.75          ! fraction of nspinsoil time used to
c                                     ! spin up soil at maximum spin up rate
c
            spincons  = 40.0          ! number of times soilbgc subroutine is
c                                     ! called for each day in simulation during
c                                     ! acceleration procedure
c
            eqyears   = 50            ! years used past spin up to slowly bring
c                                     ! entire soil to a steady equilibrium
c
            nspinsoil = iyear0 + 150  ! year of simulation that soil c reaches equilibrium
c
c
            if (soilcspin .eq. 1) then

              if ((iyear - iyear0) .le.
     >            (spinfrac * (nspinsoil - iyear0 - eqyears))) then
                 spinmax = int(spincons)
c
              else if ((iyear - iyear0) .lt.
     >               (nspinsoil - iyear0 -  eqyears)) then
c
                slope   = spincons / ((nspinsoil - iyear0 - eqyears) -
     >                    (spinfrac * (nspinsoil - iyear0 - eqyears)))
c
                spinmax = int (spincons - (slope * ((iyear - iyear0) -
     >            (spinfrac * (nspinsoil - iyear0 - eqyears)))))
c
                spinmax = max(spinmax,1)
c
              else
c
                spinmax = 1
c
              endif
c
            else 
c
                spinmax = 1
c
            endif
c
            do 230 spin = 1, spinmax
              call soilbgc (iyear, iyear0, imonth, iday, nspinsoil,
     >                      spin, spinmax)
 230        continue
c
c determine the length of a precipitation event (between 4 and 24 hours),
c and time to start and end precipitation event. plen is in timesteps, while
c plens, startp, and endp are in seconds
c
            plenmin = 1 +  int ((4.0 * 3600. - 1.) / dtime)
            plenmax = max (int (24.0 * 3600. / dtime), plenmin)
c
            plen    = min (plenmax, int (
     >                     plenmin + ran2(seed) * (plenmax-plenmin+1) ))
c
            plens   = dtime * plen
            startp  = dtime * min (niter-plen, 
     >                             int(ran2(seed)*(niter-plen+1)))
            endp    = startp + plens
c
c start of hourly loop
c
            do 240 istep = 1, niter
c
c calculate the time of day since midnight (in seconds)
c
              time = (istep - 1) * dtime
c
c determine climatic conditions for the given timestep
c
              call diurnal (time, jday, plens, startp, endp, seed)
c
c call the land surface model
c
              call lsxmain
c
c accumulate some variables every timestep
c
              call sumnow
c
              call sumday   (istep,plens)           
              call summonth (istep, iday, imonth)
              call sumyear  (istep, iday, imonth)
c
c write out diagnostics
c
              if (idiag .ne. 0) then
c
                do 250 i = 1, idiag
                  if (iyear.ge.diagstart(i) .and.
     >                iyear.le.diagend(i)   .and. 
     >                ndiagpt(i).ge.1)  then  
                    if (mod(istep,nfreq(i)).eq.0) then
                      call wdiag (i, iyear, imonth, iday, istep)
                    end if
                  end if
 250            continue
c
              end if
c
c ---------------------------------------------------------------
c               e n d   o f   t i m e   l o o p s
c ---------------------------------------------------------------
c
c end of the hourly loop
c
 240        continue
c
c write out daily output
c
            if (idailyout.eq.1) then
              write (*,*) ' '
              write (*,*) 'writing daily output'
              write (*,*) ' '
              call wdaily (nday, iyear, iyear0) 
            endif

            progress_ =  ((iyear-iy1)*365 + (imonth-1) * 30 + iday)
     >                      / 3.65/(iy2-iy1+1)
            write (*,920) progress_
c
c end of the daily loop
c
 220      continue
c
c write out monthly output (use 1st day of this month because of GrADS bug)
c
          if (imonthout.eq.1) then
            write (*,*) ' '
            write (*,*) 'writing monthly output'
            write (*,*) ' '
            call wmonthly (nday-ndaypm(imonth)+1, imonth, iyear, 
     >                     iyear0)
          endif
c
c end of the monthly loop
c
 210    continue
c
c recalculate bioclimatic parameters if using anomalies
c
        if (iyear .ge. iyranom-1 .or. iyear.ge.iyrdaily) call climanl2
c
c get new o2 and co2 concentrations for this year
c
        if (isimco2.eq.1) call co2 (co2init, co2conc, iyear)
c
c perform vegetation dynamics
c
        if (isimveg.ne.0) call dynaveg (isimfire)
c
c calculate simple annual diagnostic variables for the globe
c
        call gdiag (iyear, iyear0)
c
c calculate simple annual diagnostic variables for each vegetation type
c
        call vdiag (iyear, iyear0)
c
c write out annual output
c
        if (iyearout.eq.1) then
          write (*,*) ' '
          write (*,*) 'writing annual output'
          write (*,*) ' '
          call wyearly (nday, iyear, iyear0)
        endif
c
c write restart files
c
        call wrestart (nday, iyear, iyear0)
c
c update iyrlast value and file yearsrun.dat
c
        iyrlast = iyrlast + 1
        open(12,file= out_yearsrun_, status='unknown')
        write(12,*) iyrlast
        close(12)
c
c end of the yearly loop
c
 200  continue
c
c end of the simulation
c
      write (*,*) ' '
      progress_=100.0
      write (*,920) progress_
      write (*,*) '*** end of run ***'
      write (*,*) ' '
c
      stop
      end
c
c
c ---------------------------------------------------------------
      subroutine lsxmain
c ---------------------------------------------------------------
c
c common blocks
c 
      include 'implicit.h'
c
      include 'compar.h'
      include 'com1d.h'
      include 'comatm.h'
      include 'comsoi.h'
      include 'argvs.h'
c
c Local variables
c
      integer ib,   ! waveband number (1= visible, 2= near-IR)
     >        i     ! loop indice
c
c set physical soil quantities
c
      call setsoi
c
c calculate areal fractions wetted by intercepted h2o
c
      call fwetcal
c
c set up for solar calculations
c
      call solset
c
c solar calculations for each waveband
c
      do 100 ib = 1, nband
c
c solsur sets surface albedos for soil and snow
c solalb performs the albedo calculations
c solarf uses the unit-incident-flux results from solalb
c to obtain absorbed fluxes sol[u,s,l,g,i] and 
c incident pars sunp[u,l]
c
        call solsur (ib)
        call solalb (ib)
        call solarf (ib)
c
 100  continue
c
c calculate ir fluxes
c
      call irrad
c
c step intercepted h2o
c
      call cascade
c
c re-calculate wetted fractions, changed by cascade
c
      call fwetcal
c
c step vegetation canopy temperatures implicitly
c and calculate sensible heat and moisture fluxes
c
      call canopy
c
c step intercepted h2o due to evaporation
c
      call cascad2
c
c arbitrarily set veg temps & intercepted h2o for no-veg locations
c
      call noveg
c
c set net surface heat fluxes for soil and snow models
c
      do 110 i = 1, npoi
c
        heatg(i) = solg(i) + firg(i) - fseng(i) -
     >             hvasug(i)*fvapg(i)
c
        heati(i) = soli(i) + firi(i) - fseni(i) -
     >             hvasui(i)*fvapi(i)
c
 110  continue
c
c step snow model
c
      call snow 
c
c step soil model
c
      call soilctl
c
c return to main program
c
      return
      end
c 

      subroutine testCommon(N)
      integer N
      real :: X(N)

      return
      end

      subroutine parseArgv

      include 'argvs.h'

      integer argc_, argcI_, filePathStart_
      character*100, dimension(50):: argvs_
      character*100 fileTag_, tempStr_, inputFPath_

      argc_ = iargc()
      if(argc_ < 13) then
          print 377, 'invalid argc number, you must input at lease 14 argc'
          stop
      end if

      do argcI_=1, argc_
          tempStr_ = argvs_(argcI_)
          call getarg(argcI_, tempStr_)
          filePathStart_ = index(tempStr_, '=')
          if(filePathStart_ == 0) then
              print 377, 'invalid file path argvs, the file&
     & path must as a prefix of "--tag="!'
              stop
          end if

          fileTag_ = trim(tempStr_(:filePathStart_))
          inputFPath_ = tempStr_(filePathStart_+1:)

       if(fileTag_ == '--ini_infile_=') then
            ini_infile_ = inputFPath_
       elseif(fileTag_ == '--diag_infile_=') then
            diag_infile_ = inputFPath_
       elseif(fileTag_ == '--cld_mon_=') then
            cld_mon_ = inputFPath_
       elseif(fileTag_ == '--deltat_mon_=') then
            deltat_mon_ = inputFPath_
       elseif(fileTag_ == '--prec_mon_=') then
            prec_mon_ = inputFPath_
       elseif(fileTag_ == '--rh_mon_=') then
            rh_mon_ = inputFPath_
       elseif(fileTag_ == '--temp_mon_=') then
            temp_mon_ = inputFPath_
       elseif(fileTag_ == '--trange_mon_=') then
            trange_mon_ = inputFPath_
       elseif(fileTag_ == '--wetd_mon_=') then
            wetd_mon_ = inputFPath_
       elseif(fileTag_ == '--wspd_mon_=') then
            wspd_mon_ = inputFPath_
       elseif(fileTag_ == '--soita_sand_=') then
            soita_sand_ = inputFPath_
       elseif(fileTag_ == '--soita_clay_=') then
            soita_clay_ = inputFPath_
       elseif(fileTag_ == '--vegtype_=') then
            vegtype_ = inputFPath_
       elseif(fileTag_ == '--surta_=') then
            surta_ = inputFPath_
       elseif(fileTag_ == '--topo_=') then
            topo_ = inputFPath_
       elseif(fileTag_ == '--params_can_=') then
            params_can_ = inputFPath_
       elseif(fileTag_ == '--params_hyd_=') then
            params_hyd_ = inputFPath_
       elseif(fileTag_ == '--params_soi_=') then
            params_soi_ = inputFPath_
       elseif(fileTag_ == '--params_veg_=') then
            params_veg_ = inputFPath_
       elseif(fileTag_ == '--out_yearly_aet_=') then
            out_yearly_aet_ = inputFPath_
       elseif(fileTag_ == '--out_yearly_biomass_=') then
            out_yearly_biomass_ = inputFPath_
       elseif(fileTag_ == '--out_yearly_co2fluxes_=') then
            out_yearly_co2fluxes_ = inputFPath_
       elseif(fileTag_ == '--out_yearly_csoi_=') then
            out_yearly_csoi_ = inputFPath_
       elseif(fileTag_ == '--out_yearly_disturbf_=') then
            out_yearly_disturbf_ = inputFPath_
       elseif(fileTag_ == '--out_yearly_exist_=') then
            out_yearly_exist_ = inputFPath_
       elseif(fileTag_ == '--out_yearly_fcover_=') then
            out_yearly_fcover_ = inputFPath_
       elseif(fileTag_ == '--out_yearly_npp_=') then
            out_yearly_npp_ = inputFPath_
       elseif(fileTag_ == '--out_yearly_nsoi_=') then
            out_yearly_nsoi_ = inputFPath_
       elseif(fileTag_ == '--out_yearly_plai_=') then
            out_yearly_plai_ = inputFPath_
       elseif(fileTag_ == '--out_yearly_runoff_=') then
            out_yearly_runoff_ = inputFPath_
       elseif(fileTag_ == '--out_yearly_sens_=') then
            out_yearly_sens_ = inputFPath_
       elseif(fileTag_ == '--out_yearly_tsoi_=') then
            out_yearly_tsoi_ = inputFPath_
       elseif(fileTag_ == '--out_yearly_vegtype0_=') then
            out_yearly_vegtype0_ = inputFPath_
       elseif(fileTag_ == '--out_yearly_wsoi_=') then
            out_yearly_wsoi_ = inputFPath_
       elseif(fileTag_ == '--out_yearly_zcanopy_=') then
            out_yearly_zcanopy_ = inputFPath_
       elseif(fileTag_ == '--out_yearly_sapfrac_=') then
            out_yearly_sapfrac_ = inputFPath_
       elseif(fileTag_ == '--out_yearly_dummyv_=') then
            out_yearly_dummyv_ = inputFPath_
       elseif(fileTag_ == '--out_yearly_solar_=') then
            out_yearly_solar_ = inputFPath_
       elseif(fileTag_ == '--out_yearly_albedo_=') then
            out_yearly_albedo_ = inputFPath_
       elseif(fileTag_ == '--out_yearly_latent_=') then
            out_yearly_latent_ = inputFPath_
       elseif(fileTag_ == '--out_yearly_totfall_=') then
            out_yearly_totfall_ = inputFPath_
       elseif(fileTag_ == '--out_yearly_clitw_=') then
            out_yearly_clitw_ = inputFPath_
       elseif(fileTag_ == '--out_yearly_csoislo_=') then
            out_yearly_csoislo_ = inputFPath_
       elseif(fileTag_ == '--out_yearly_csoipas_=') then
            out_yearly_csoipas_ = inputFPath_
       elseif(fileTag_ == '--out_monthly_aet_=') then
            out_monthly_aet_ = inputFPath_
       elseif(fileTag_ == '--out_monthly_cloud_=') then
            out_monthly_cloud_ = inputFPath_
       elseif(fileTag_ == '--out_monthly_co2ratio_=') then
            out_monthly_co2ratio_ = inputFPath_
       elseif(fileTag_ == '--out_monthly_ir_=') then
            out_monthly_ir_ = inputFPath_
       elseif(fileTag_ == '--out_monthly_lai_=') then
            out_monthly_lai_ = inputFPath_
       elseif(fileTag_ == '--out_monthly_latent_=') then
            out_monthly_latent_ = inputFPath_
       elseif(fileTag_ == '--out_monthly_npptot_=') then
            out_monthly_npptot_ = inputFPath_
       elseif(fileTag_ == '--out_monthly_qa_=') then
            out_monthly_qa_ = inputFPath_
       elseif(fileTag_ == '--out_monthly_rain_=') then
            out_monthly_rain_ = inputFPath_
       elseif(fileTag_ == '--out_monthly_rh_=') then
            out_monthly_rh_ = inputFPath_
       elseif(fileTag_ == '--out_monthly_runoff_=') then
            out_monthly_runoff_ = inputFPath_
       elseif(fileTag_ == '--out_monthly_sens_=') then
            out_monthly_sens_ = inputFPath_
       elseif(fileTag_ == '--out_monthly_snod_=') then
            out_monthly_snod_ = inputFPath_
       elseif(fileTag_ == '--out_monthly_snof_=') then
            out_monthly_snof_ = inputFPath_
       elseif(fileTag_ == '--out_monthly_snow_=') then
            out_monthly_snow_ = inputFPath_
       elseif(fileTag_ == '--out_monthly_solar_=') then
            out_monthly_solar_ = inputFPath_
       elseif(fileTag_ == '--out_monthly_temp_=') then
            out_monthly_temp_ = inputFPath_
       elseif(fileTag_ == '--out_monthly_tsoi_=') then
            out_monthly_tsoi_ = inputFPath_
       elseif(fileTag_ == '--out_monthly_wsoi_=') then
            out_monthly_wsoi_ = inputFPath_
       elseif(fileTag_ == '--out_monthly_albedo_=') then
            out_monthly_albedo_ = inputFPath_
       elseif(fileTag_ == '--out_monthly_dummyv_=') then
            out_monthly_dummyv_ = inputFPath_
       elseif(fileTag_ == '--out_daily_rain_=') then
            out_daily_rain_ = inputFPath_
       elseif(fileTag_ == '--out_daily_cloud_=') then
            out_daily_cloud_ = inputFPath_
       elseif(fileTag_ == '--out_daily_rh_=') then
            out_daily_rh_ = inputFPath_
       elseif(fileTag_ == '--out_daily_snow_=') then
            out_daily_snow_ = inputFPath_
       elseif(fileTag_ == '--out_daily_aet_=') then
            out_daily_aet_ = inputFPath_
       elseif(fileTag_ == '--out_daily_trunoff_=') then
            out_daily_trunoff_ = inputFPath_
       elseif(fileTag_ == '--out_daily_srunoff_=') then
            out_daily_srunoff_ = inputFPath_
       elseif(fileTag_ == '--out_daily_drainage_=') then
            out_daily_drainage_ = inputFPath_
       elseif(fileTag_ == '--out_daily_wsoi_=') then
            out_daily_wsoi_ = inputFPath_
       elseif(fileTag_ == '--out_daily_wisoi_=') then
            out_daily_wisoi_ = inputFPath_
       elseif(fileTag_ == '--out_daily_snod_=') then
            out_daily_snod_ = inputFPath_
       elseif(fileTag_ == '--out_daily_snof_=') then
            out_daily_snof_ = inputFPath_
       elseif(fileTag_ == '--out_daily_co2ratio_=') then
            out_daily_co2ratio_ = inputFPath_
       elseif(fileTag_ == '--out_daily_co2mic_=') then
            out_daily_co2mic_ = inputFPath_
       elseif(fileTag_ == '--out_daily_templ_=') then
            out_daily_templ_ = inputFPath_
       elseif(fileTag_ == '--out_daily_zcanopy_=') then
            out_daily_zcanopy_ = inputFPath_
       elseif(fileTag_ == '--out_daily_laicanopy_=') then
            out_daily_laicanopy_ = inputFPath_
       elseif(fileTag_ == '--temp_danom_=') then
            temp_danom_ = inputFPath_
       elseif(fileTag_ == '--trange_danom_=') then
            trange_danom_ = inputFPath_
       elseif(fileTag_ == '--prec_danom_=') then
            prec_danom_ = inputFPath_
       elseif(fileTag_ == '--cld_danom_=') then
            cld_danom_ = inputFPath_
       elseif(fileTag_ == '--rh_danom_=') then
            rh_danom_ = inputFPath_
       elseif(fileTag_ == '--wspd_danom_=') then
            wspd_danom_ = inputFPath_
       elseif(fileTag_ == '--wetd_danom_=') then
            wetd_danom_ = inputFPath_
       elseif(fileTag_ == '--prec_daily_=') then
            prec_daily_ = inputFPath_
       elseif(fileTag_ == '--temp_daily_=') then
            temp_daily_ = inputFPath_
       elseif(fileTag_ == '--trange_daily_=') then
            trange_daily_ = inputFPath_
       elseif(fileTag_ == '--cld_daily_=') then
            cld_daily_ = inputFPath_
       elseif(fileTag_ == '--wspd_daily_=') then
            wspd_daily_ = inputFPath_
       elseif(fileTag_ == '--sphum_daily_=') then
            sphum_daily_ = inputFPath_
       elseif(fileTag_ == '--prec_fanom_=') then
            prec_fanom_ = inputFPath_
       elseif(fileTag_ == '--trange_fanomc_=') then
            trange_fanomc_ = inputFPath_
       elseif(fileTag_ == '--wspd_fanomc_=') then
            wspd_fanomc_ = inputFPath_
       elseif(fileTag_ == '--sphum_fanom_=') then
            sphum_fanom_ = inputFPath_
       elseif(fileTag_ == '--deltat_=') then
            deltat_ = inputFPath_
       elseif(fileTag_ == '--out_global_=') then
            out_global_ = inputFPath_
       elseif(fileTag_ == '--out_vegtype_=') then
            out_vegtype_ = inputFPath_
       elseif(fileTag_ == '--out_yearsrun_=') then
            out_yearsrun_ = inputFPath_
       elseif(fileTag_ == '--out_diag_0_=') then
            out_diag_0_ = inputFPath_
       elseif(fileTag_ == '--out_diag_1_=') then
            out_diag_1_ = inputFPath_
       elseif(fileTag_ == '--out_diag_2_=') then
            out_diag_2_ = inputFPath_
       elseif(fileTag_ == '--out_diag_3_=') then
            out_diag_3_ = inputFPath_
       elseif(fileTag_ == '--out_diag_4_=') then
            out_diag_4_ = inputFPath_
       elseif(fileTag_ == '--out_diag_5_=') then
            out_diag_5_ = inputFPath_
       elseif(fileTag_ == '--out_diag_6_=') then
            out_diag_6_ = inputFPath_
       elseif(fileTag_ == '--out_diag_7_=') then
            out_diag_7_ = inputFPath_
       elseif(fileTag_ == '--out_diag_8_=') then
            out_diag_8_ = inputFPath_
       elseif(fileTag_ == '--out_diag_9_=') then
            out_diag_9_ = inputFPath_
       end if
      end do

      if(len_trim(ini_infile_) == 0 .or.
     >   len_trim(diag_infile_) == 0 .or.
     >   len_trim(cld_mon_) == 0 .or.
     >   len_trim(deltat_mon_) == 0 .or.
     >   len_trim(prec_mon_) == 0 .or.
     >   len_trim(rh_mon_) == 0 .or.
     >   len_trim(temp_mon_) == 0 .or.
     >   len_trim(trange_mon_) == 0 .or.
     >   len_trim(wetd_mon_) == 0 .or.
     >   len_trim(wspd_mon_) == 0 .or.
     >   len_trim(soita_sand_) == 0 .or.
     >   len_trim(soita_clay_) == 0 .or.
     >   len_trim(vegtype_) == 0 .or.
     >   len_trim(surta_) == 0 .or.
     >   len_trim(topo_) == 0) then
          print 377, 'invalid input list!'
          stop
      end if

      return
      end
