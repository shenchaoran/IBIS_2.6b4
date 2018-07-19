c io.f  last update 7/15/99 dtp
c
c    #     ####
c    #    #    #
c    #    #    #
c    #    #    #
c    #    #    #
c    #     ####
c This file contains these subroutines
c wrestart
c wdaily
c wmonthly
c wyearly
c readit
c restart
c coldstart
c rdanom
c inird
c rdday
c diaginit
c wdiag
c
c See end of file for information on how to add new code to read a
c file or write a file.
c
c ---------------------------------------------------------------------
      subroutine wrestart (nday, iyear, iyear0)
c ---------------------------------------------------------------------
c
c this subroutine writes the restart values of:
c
c  fsnocov = fractional snow cover
c  tsno    = temperature of snow
c  hsno    = snow depth
c  tsoi    = soil temperature
c  wisoi   = soil ice content
c  wsoi    = soil moisture content
c  cbiol   = carbon in leaf biomass pool
c  cbiow   = carbon in woody biomass pool
c  cbior   = carbon in fine root biomass pool
c  sapfrac = sapwood fraction
c  clitlm  = leaf metabolic litter
c  clitls  = leaf structural litter
c  clitll  = leaf lignin litter
c  clitrm  = root metabolic litter
c  clitrs  = root structural litter
c  clitrl  = root lignin litter
c  clitwm  = woody metabolic litter
c  clitws  = woody structural litter
c  clitwl  = woody lignin litter
c  falll   = annual leaf litterfall
c  fallr   = annual fine root turnover 
c  fallw   = annual wood litterfall
c  totcmic = total microbial carbon
c  csoislop= slow soil carbon, protected humus
c  csoislon= slow soil carbon, nonprotected humus
c  csoipas = passive soil carbon
c  gdd0    = growing degree days 0
c  gdd5    = growing degree days 5
c  tc      = coldest monthly temperature
c  tw      = warmest monthly temperature
c  wipud   = ice content of puddles per soil area
c  wpud    = liquid content of puddles per soil area
c  agddu   = annual accumulated growing degree days for bud burst, upper canopy
c  agddl   = annual accumulated growing degree days for bud burst, lower canopy
c  tempu   = cold-phenology trigger for trees
c  templ   = cold-phenology trigger for grasses/shrubs
c  a10td    = 10-day avg daily temp
c  a10ancub = 10-day average canopy photosynthesis rate - broadleaf
c  a10ancuc = 10-day average canopy photosynthesis rate - conifers
c  a10ancls = 10-day average canopy photosynthesis rate - shrubs
c  a10ancl4 = 10-day average canopy photosynthesis rate - c4 grasses
c  a10ancl3 = 10-day average canopy photosynthesis rate - c3 grasses
c  a10scalparamu = 10-day average canopy scaling parameter - upper canopy
c  a10scalparaml = 10-day average canopy scaling parameter - lower canopy
c  a10daylightu = 10-day average daylight - upper canopy
c  a10daylightl = 10-day average daylight - lower canopy
c  dropu   = drought-phenology trigger for trees
c  dropls  = drought-phenology trigger for shrubs
c  dropl4  = drought-phenology trigger for c4 grasses
c  dropl3  = drought-phenology trigger for c3 grasses
c (NOTE: a10ancuc is not used at this point, so its wrestart entry 
c is commented out)
c ---------------------------------------------------------------------
c
      include 'implicit.h'
c
c instantaneous output for restarts
c
      include 'compar.h'
      include 'comsoi.h'
      include 'comsno.h'
      include 'comsum.h'
      include 'comveg.h'
      include 'comwork.h'
      include 'argvs.h'
c
c Arguments
c
      integer nday,         ! number of days run since iyear0
     >        iyear,        ! this calendar year
     >        iyear0        ! initial year

c
c local variables
c
      integer lf,           ! number of characters in directory name
     >        n, k,         ! loop indices
     >        idies,        ! file indice (?) for netcdf
     >        istat,        ! error flag for netcdf      
     >        nyears        ! years of run iyear0
c
      integer istart(4),
     >        icount(4)         ! for writing restart vars
c
      character*21 tunits
      character*50 fdir         ! used to construct odd/even file names
      character*11 cdate        ! date to use in history attribute in files
      character*10 tdate        ! character date for time step
      character*80 dimnames(4)  ! names of dimensions for restart vars
      character*80 pftdef(npft) ! plant functional type defs (not used now)
      character*80 filen        ! file name
c
      real slayers(nsnolay),    ! index for snow layers
     >     depthsoi(nsoilay),   ! soil layer depths
     >     pindex(npft),        ! index for pfts
     >     ftime,               ! floating point time value
     >     tweight              ! time weight (# days/sample)
c
c External
c
      integer lenchr,           ! Function: Find length of character string
     > NF_PUT_ATT_TEXT,         ! netcdf function
     > NF_GLOBAL                ! netcdf function
c
c use workspace for local real variables
c
      equivalence (slayers(1),work(1)),(depthsoi(1),work(ndim4+1))
      equivalence (ftime,work(2*ndim4+1)), (tweight,work(3*ndim4+1))
      equivalence (pindex(1),work(4*ndim4+1))
c
c ---------------------------------------------------------------------
c
      data istart / 1,1,1,1 /, 
     >     icount / nlon,nlat,1,1 /
      icount(1) = nlonsub
      icount(2) = nlatsub
c
      tweight = 1.
c
c check to see if iyear is odd or even, construct appropriate file name root
c
      if (mod(iyear,2) .eq. 0) then
         fdir = 'restart/even'
      else
         fdir = 'restart/odd'
      end if
      lf = lenchr(fdir)
c
c tdate is december of this year, max year = 999
c
      tdate='DEC000'//char(0)//char(0)//char(0)//char(0)
      nyears = iyear - iyear0 + 1
      if (nyears .lt. 10) then
         write(tdate(6:6),'(i1)') nyears
       else if (nyears .lt. 100) then
         write(tdate(5:6),'(i2)') nyears
      else
         write(tdate(4:6),'(i3)') nyears
      end if
c
c initialize snow layer indicies, pft names, etc
c
      if (nyears .le. 2) then
         call simdate(cdate)
         ftime = nday
c
c time units is days since Dec 31 of the year before iyear0
c
         tunits = 'days since 0000-12-31'
         write(tunits(12:15),'(i4)') iyear0-1
c
         do 5 n = 1, nsnolay
            slayers(n) = float(n)
 5       continue
c
         depthsoi(1) = hsoi(1)
         do 10 n = 2, nsoilay
            depthsoi(n) = depthsoi(n-1)+hsoi(n)
 10      continue
c
c define pft by index
c
         do 12 n = 1, npft
            pindex(n) = n
 12      continue
c
c and by character label
c
         pftdef(1) = 'trbrevtr - tropical broadleaf evergreen trees'
     >    //char(0)
         pftdef(2) = 
     >    'trbrdetr - tropical broadleaf drought-deciduous trees'
     >    //char(0)
         pftdef(3) = 
     >    'wtbrevtr - warm-temperate broadleaf evergreen trees'
     >    //char(0)
         pftdef(4) = 'tecoevtr - temperate conifer evergreen trees'
     >    //char(0)
         pftdef(5) = 
     >    'tebrdetr - temperate broadleaf cold-deciduous trees'//char(0)
         pftdef(6) = 'bocoevtr - boreal conifer evergreen trees'
     >    //char(0)
         pftdef(7) = 'bocodetr - boreal conifer cold-deciduous trees'
     >    //char(0)
         pftdef(8) = 
     >    'bobrdetr - boreal broadleaf cold-deciduous trees'//char(0)
         pftdef(9) = 'evsh - evergreen shrubs'//char(0)
         pftdef(10) = 'desh - deciduous shrubs'//char(0)
         pftdef(11) = 'c4gr - warm (c4) grasses'//char(0)
         pftdef(12) = 'c3gr - cool (c3) grasses'//char(0)
c
         dimnames(1) = 'longitude'
         dimnames(2) = 'latitude'
c
c dimnames(3) is set for each variable seperately
c
         dimnames(4) = 'time'
      end if
cc
cc dummy variable example, 3-d - copy & modify for new variable.
cc If new variable is 4-d, see tsno to see how 4-d differs from 3-d.
cc
c      filen = fdir(1:lf)//'dummyv.nc'
c      if (nyears .le. 2) then
c         call inifile(idies,filen,
c     >    'restart file for dummyv',
c     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
c     >    'none',1,0.,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'dummyv',
c     >    'instantaneous dummyv','dummyvs-units',3,dimnames,
c     >    OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      call vec2arr (dummyv, cdummy)
c      icount(3) = 1
c      call writevar(filen,'dummyv',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wrestart, dummyv'
c         stop 1
c      end if
c
c fractional snow cover
c
      filen = fdir(1:lf)//'fsnocov.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for fractional snow cover',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'fsnocov',
     >    'instantaneous fractional snow cover','fraction',3,dimnames,
     >    OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (fi, cdummy)
      icount(3) = 1
      call writevar(filen,'fsnocov',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, fsnocov'
         stop 1
      end if
c
c temperature of snow layers
c
      filen = fdir(1:lf)//'tsno.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for snow temperature',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,latscale,
     >    'snowlayer','snow layers top to bottom','',nsnolay,slayers,
     >    'down',tunits,'gregorian',istat)
         dimnames(3) = 'snowlayer'
         call inivar(idies,'tsno',
     >    'instantaneous snow cover temperature','degK',4,dimnames,
     >    OCEAN,istat)
         call endini(idies,istat)
      end if
      do 15 k = 1, nsnolay
         call vec2arr (tsno(1,k), cdummy((k-1)*nlonsub*nlatsub + 1))
 15   continue
      icount(3) = nsnolay
      call writevar(filen,'tsno',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, tsno'
         stop 1
      end if
c
c thickness of snow layers
c
      filen = fdir(1:lf)//'hsno.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for snow layer thickness',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,latscale,
     >    'snowlayer','snow layers top to bottom','',nsnolay,slayers,
     >    'down',tunits,'gregorian',istat)
         dimnames(3) = 'snowlayer'
         call inivar(idies,'hsno','instantaneous snow layer thickness',
     >    'meters',4,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      do 20 k = 1, nsnolay
         call vec2arr (hsno(1,k), cdummy((k-1)*nlonsub*nlatsub + 1))
 20   continue
      icount(3) = nsnolay
      call writevar(filen,'hsno',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, hsno'
         stop 1
      end if
c
c temperature of soil layers
c
      filen = fdir(1:lf)//'tsoi.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for soil temperature',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,latscale,
     >    'soillayer','depth of soil layer bottom','meter',nsoilay,
     >    depthsoi,'down',tunits,'gregorian',istat)
         dimnames(3) = 'soillayer'
         call inivar(idies,'tsoi','instantaneous soil temperature',
     >    'degK',4,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      do 25 k = 1, nsoilay
         call vec2arr (tsoi(1,k), cdummy((k-1)*nlonsub*nlatsub + 1))
 25   continue
      icount(3) = nsoilay
      call writevar(filen,'tsoi',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, tsoi'
         stop 1
      end if
c
c ice content of soil
c
      filen = fdir(1:lf)//'wisoi.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for soil ice content',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,latscale,
     >    'soillayer','depth of soil layer bottom','meter',nsoilay,
     >    depthsoi,'down',tunits,'gregorian',istat)
         dimnames(3) = 'soillayer'
         call inivar(idies,'wisoi',
     >    'instantaneous fraction of soil pore space containing ice',
     >    'fraction',4,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      do 30 k = 1, nsoilay
         call vec2arr (wisoi(1,k), cdummy((k-1)*nlonsub*nlatsub + 1))
 30   continue
      icount(3) = nsoilay
      call writevar(filen,'wisoi',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, wisoi'
         stop 1
      end if
c
c water content of soil
c
      filen = fdir(1:lf)//'wsoi.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for soil water content',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,latscale,
     >    'soillayer','depth of soil layer bottom','meter',nsoilay,
     >    depthsoi,'down',tunits,'gregorian',istat)
         dimnames(3) = 'soillayer'
         call inivar(idies,'wsoi',
     >    'instantaneous fraction of soil pore space containing water',
     >    'fraction',4,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      do 35 k = 1, nsoilay
         call vec2arr (wsoi(1,k), cdummy((k-1)*nlonsub*nlatsub + 1))
 35   continue
      icount(3) = nsoilay
      call writevar(filen,'wsoi',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, wsoi'
         stop 1
      end if
c
c carbon in leaf biomass pool
c
      filen = fdir(1:lf)//'cbiol.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for carbon in leaf biomass pool',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,latscale,
     >    'pft','plant functional type','_',npft,pindex,'',
     >    tunits,'gregorian',istat)
c add global attribute to define pfts with text, use netcdf low-level command
         istat = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'pft_definition',
     >    npft*80,pftdef)
         dimnames(3) = 'pft'
         call inivar(idies,'cbiol',
     >    'instantaneous carbon in leaf biomass pool','kg/m^2',
     >    4,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      do 40 k = 1, npft
         call vec2arr (cbiol(1,k), cdummy((k-1)*nlonsub*nlatsub + 1))
 40   continue
      icount(3) = npft
      call writevar(filen,'cbiol',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, cbiol'
         stop 1
      end if
c
c carbon in wood
c
      filen = fdir(1:lf)//'cbiow.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for carbon in wood biomass pool',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,latscale,
     >    'pft','plant functional type','_',npft,pindex,'',
     >    tunits,'gregorian',istat)
c add global attribute to define pfts with text, use netcdf low-level command
         istat = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'pft_definition',
     >    npft*80,pftdef)
         dimnames(3) = 'pft'
         call inivar(idies,'cbiow',
     >    'instantaneous carbon in wood biomass pool','kg/m^2',
     >    4,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      do 45 k = 1, npft
         call vec2arr (cbiow(1,k), cdummy((k-1)*nlonsub*nlatsub + 1))
 45   continue
      icount(3) = npft
      call writevar(filen,'cbiow',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, cbiow'
         stop 1
      end if
c
c carbon in root
c
      filen = fdir(1:lf)//'cbior.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for carbon in root biomass pool',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,latscale,
     >    'pft','plant functional type','_',npft,pindex,'',
     >    tunits,'gregorian',istat)
c add global attribute to define pfts with text, use netcdf low-level command
         istat = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'pft_definition',
     >    npft*80,pftdef)
         dimnames(3) = 'pft'
         call inivar(idies,'cbior',
     >    'instantaneous carbon in root biomass pool','kg/m^2',
     >    4,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      do 50 k = 1, npft
         call vec2arr (cbior(1,k), cdummy((k-1)*nlonsub*nlatsub + 1))
 50   continue
      icount(3) = npft
      call writevar(filen,'cbior',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, cbior'
         stop 1
      end if
c
c sapwood fraction
c
      filen = fdir(1:lf)//'sapfrac.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for sapwood fraction',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'sapfrac','instantaneous sapwood fraction',
     >    'fraction',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (sapfrac, cdummy)
      icount(3) = 1
      call writevar(filen,'sapfrac',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, sapfrac'
         stop 1
      end if
c
c leaf metabolic litter
c
      filen = fdir(1:lf)//'clitlm.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for leaf metabolic litter',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'clitlm',
     >    'instantaneous leaf metabolic litter carbon','kg/m^2',
     >    3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (clitlm, cdummy)
      icount(3) = 1
      call writevar(filen,'clitlm',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, clitlm'
         stop 1
      end if
c
c leaf structural carbon litter
c
      filen = fdir(1:lf)//'clitls.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for leaf structural litter',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'clitls',
     >    'instantaneous leaf structural litter carbon','kg/m^2',
     >    3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (clitls, cdummy)
      icount(3) = 1
      call writevar(filen,'clitls',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, clitls'
         stop 1
      end if
c
c leaf lignin carbon litter
c
      filen = fdir(1:lf)//'clitll.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for leaf lignin litter',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'clitll',
     >    'instantaneous leaf lignin litter carbon','kg/m^2',
     >    3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (clitll, cdummy)
      icount(3) = 1
      call writevar(filen,'clitll',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, clitll'
         stop 1
      end if
c
c root metabolic litter
c
      filen = fdir(1:lf)//'clitrm.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for root metabolic litter',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'clitrm',
     >    'instantaneous root metabolic litter carbon','kg/m^2',
     >    3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (clitrm, cdummy)
      icount(3) = 1
      call writevar(filen,'clitrm',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, clitrm'
         stop 1
      end if
c
c root structural litter
c
      filen = fdir(1:lf)//'clitrs.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for root structural litter',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'clitrs',
     >    'instantaneous root structural litter carbon','kg/m^2',
     >    3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (clitrs, cdummy)
      icount(3) = 1
      call writevar(filen,'clitrs',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, clitrs'
         stop 1
      end if
c
c root lignin litter
c
      filen = fdir(1:lf)//'clitrl.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for root lignin litter',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'clitrl',
     >    'instantaneous root lignin litter carbon','kg/m^2',
     >    3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (clitrl, cdummy)
      icount(3) = 1
      call writevar(filen,'clitrl',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, clitrl'
         stop 1
      end if
c
c woody metabolic litter
c
      filen = fdir(1:lf)//'clitwm.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for woody metabolic litter',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'clitwm',
     >    'instantaneous woody metabolic litter carbon','kg/m^2',
     >    3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (clitwm, cdummy)
      icount(3) = 1
      call writevar(filen,'clitwm',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, clitwm'
         stop 1
      end if
c
c woody structural litter
c
      filen = fdir(1:lf)//'clitws.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for woody structural litter',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'clitws',
     >    'instantaneous woody structural litter carbon','kg/m^2',
     >    3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (clitws, cdummy)
      icount(3) = 1
      call writevar(filen,'clitws',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, clitws'
         stop 1
      end if
c
c woody lignin litter
c
      filen = fdir(1:lf)//'clitwl.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for woody lignin litter',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'clitwl',
     >    'instantaneous woody lignin litter carbon','kg/m^2',
     >    3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (clitwl, cdummy)
      icount(3) = 1
      call writevar(filen,'clitwl',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, clitwl'
         stop 1
      end if
c
c annual leaf litterfall
c
      filen = fdir(1:lf)//'falll.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for annual leaf litterfall',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'falll',
     >    'annual leaf litterfall carbon','kg/m^2',
     >    3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (falll, cdummy)
      icount(3) = 1
      call writevar(filen,'falll',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, falll'
         stop 1
      end if
c
c annual fine root turnover 
c
      filen = fdir(1:lf)//'fallr.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for annual fine root turnover',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'fallr',
     >    'annual fine root turnover carbon','kg/m^2',
     >    3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (fallr, cdummy)
      icount(3) = 1
      call writevar(filen,'fallr',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, fallr'
         stop 1
      end if
c
c annual wood turnover 
c
      filen = fdir(1:lf)//'fallw.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for annual woody turnover',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'fallw',
     >    'annual wood turnover carbon','kg/m^2',
     >    3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (fallw, cdummy)
      icount(3) = 1
      call writevar(filen,'fallw',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, fallw'
         stop 1
      end if
c
c total microbial carbon
c
      filen = fdir(1:lf)//'totcmic.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for total microbial carbon',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'totcmic',
     >    'instantaneous total microbial carbon','kg/m^2',
     >    3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (totcmic, cdummy)
      icount(3) = 1
      call writevar(filen,'totcmic',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, totcmic'
         stop 1
      end if
c
c slow soil carbon, protected humus
c
      filen = fdir(1:lf)//'csoislop.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for slow soil carbon, protected humus',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'csoislop',
     >    'instantaneous slow soil carbon protected humus',
     >    'kg/m^2',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (csoislop, cdummy)
      icount(3) = 1
      call writevar(filen,'csoislop',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, csoislop'
         stop 1
      end if
c
c slow soil carbon, nonprotected humus
c
      filen = fdir(1:lf)//'csoislon.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for slow soil carbon, nonprotected humus',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'csoislon',
     >    'instantaneous slow soil carbon nonprotected humus',
     >    'kg/m^2',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (csoislon, cdummy)
      icount(3) = 1
      call writevar(filen,'csoislon',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, csoislon'
         stop 1
      end if
c
c passive soil carbon
c
      filen = fdir(1:lf)//'csoipas.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for passive soil carbon',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'csoipas',
     >    'instantaneous passive soil carbon','kg/m^2',
     >    3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (csoipas, cdummy)
      icount(3) = 1
      call writevar(filen,'csoipas',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, csoipas'
         stop 1
      end if
c
c growing degree days
c
      filen = fdir(1:lf)//'gdd0.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for growing degree days above 0 deg_C',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'gdd0',
     >    'instantaneous growing degree days above 0 deg_C',
     >    'days degC',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (gdd0, cdummy)
      icount(3) = 1
      call writevar(filen,'gdd0',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, gdd0'
         stop 1
      end if
c
      filen = fdir(1:lf)//'gdd5.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for growing degree days above 5 deg_C',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'gdd5',
     >    'instantaneous growing degree days above 5 deg_C',
     >    'days degC',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (gdd5, cdummy)
      icount(3) = 1
      call writevar(filen,'gdd5',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, gdd5'
         stop 1
      end if
c
c coldest monthly temperature
c
      filen = fdir(1:lf)//'tc.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for coldest monthly temperature',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'tc',
     >    'instantaneous coldest monthly temperature',
     >    'degC',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (tc, cdummy)
      icount(3) = 1
      call writevar(filen,'tc',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, tc'
         stop 1
      end if
c
c warmest monthly temperature
c
      filen = fdir(1:lf)//'tw.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for warmest monthly temperature',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'tw',
     >    'instantaneous warmest monthly temperature',
     >    'degC',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (tw, cdummy)
      icount(3) = 1
      call writevar(filen,'tw',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, tw'
         stop 1
      end if
c
c ice content of puddles
c
      filen = fdir(1:lf)//'wipud.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for ice content of puddles',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'wipud',
     >    'instantaneous ice content of puddles',
     >    'kg/m^2',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (wipud, cdummy)
      icount(3) = 1
      call writevar(filen,'wipud',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, wipud'
         stop 1
      end if
c
c liquid content of puddles
c
      filen = fdir(1:lf)//'wpud.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for liquid content of puddles',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'wpud',
     >    'instantaneous liquid water content of puddles',
     >    'kg/m^2',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (wpud, cdummy)
      icount(3) = 1
      call writevar(filen,'wpud',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, wpud'
         stop 1
      end if
c
c annual accumulated growing degree days for bud burst, upper canopy
c
      filen = fdir(1:lf)//'agddu.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for upper canopy growing degree days',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'agddu',
     >    'instantaneous growing degree days uc',
     >    'days degC',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (agddu, cdummy)
      icount(3) = 1
      call writevar(filen,'agddu',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, agddu'
         stop 1
      end if
c
c annual accumulated growing degree days for bud burst, lower canopy
c
      filen = fdir(1:lf)//'agddl.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for lower canopy growing degree days',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'agddl',
     >    'instantaneous growing degree days lc',
     >    'days degC',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (agddl, cdummy)
      icount(3) = 1
      call writevar(filen,'agddl',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, agddl'
         stop 1
      end if
c
c cold-phenology trigger for trees
c
      filen = fdir(1:lf)//'tempu.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for cold phenology trigger for trees',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'tempu',
     >    'cold phenology trigger for trees',
     >    '_',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (tempu, cdummy)
      icount(3) = 1
      call writevar(filen,'tempu',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, tempu'
         stop 1
      end if
c
c cold-phenology trigger for grasses/shrubs
c
      filen = fdir(1:lf)//'templ.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for cold phenology trigger for grasses/shrubs',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'templ',
     >    'cold phenology trigger for grasses/shrubs',
     >    '_',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (templ, cdummy)
      icount(3) = 1
      call writevar(filen,'templ',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, templ'
         stop 1
      end if
c
c 10-day average daily air temperature
c
      filen = fdir(1:lf)//'a10td.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for 10-day average daily air T',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'a10td',
     >    '10-day average daily air T',
     >    'K',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (a10td, cdummy)
      icount(3) = 1
      call writevar(filen,'a10td',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, a10td'
         stop 1
      end if
c
c 10-day average canopy photosynthesis rate - broadleaf
c
      filen = fdir(1:lf)//'a10ancub.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for 10-day average canopy photosynth. rate,'
     >    //' broadleaf',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'a10ancub',
     >    '10-day average canopy photosynth. rate, broadleaf',
     >    'mol_co2 m-2 s-1',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (a10ancub, cdummy)
      icount(3) = 1
      call writevar(filen,'a10ancub',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, a10ancub'
         stop 1
      end if
ccc
ccc 10-day average canopy photosynthesis rate - conifer
ccc
cc      filen = fdir(1:lf)//'a10ancuc.nc'
cc      if (nyears .le. 2) then
cc         call inifile(idies,filen,
cc     >    'restart file for 10-day average canopy photosynth. rate,'
cc     >    //' conifer',
cc     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
cc     >    latscale,'','none',
cc     >    'none',1,0.,'',tunits,'gregorian',istat)
cc         dimnames(3) = 'time'
cc         call inivar(idies,'a10ancuc',
cc     >    '10-day average canopy photosynth. rate, conifer',
cc     >    'mol_co2 m-2 s-1',3,dimnames,OCEAN,istat)
cc         call endini(idies,istat)
cc      end if
cc      call vec2arr (a10ancuc, cdummy)
cc      icount(3) = 1
cc      call writevar(filen,'a10ancuc',istart,icount,cdummy,ftime,
cc     > tweight,tdate,istat)
cc      if (istat .ne. 0) then
cc         write(*,*) 'ERROR in wrestart, a10ancuc'
cc         stop 1
cc      end if
c
c 10-day average canopy photosynthesis rate - shrubs
c
      filen = fdir(1:lf)//'a10ancls.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for 10-day average canopy photosynth. rate,'
     >    //' shrubs',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'a10ancls',
     >    '10-day average canopy photosynth. rate, shrubs',
     >    'mol_co2 m-2 s-1',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (a10ancls, cdummy)
      icount(3) = 1
      call writevar(filen,'a10ancls',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, a10ancls'
         stop 1
      end if
c
c 10-day average canopy photosynthesis rate - c4 grasses
c
      filen = fdir(1:lf)//'a10ancl4.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for 10-day average canopy photosynth. rate,'
     >    //' c4 grasses',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'a10ancl4',
     >    '10-day average canopy photosynth. rate, c4 grasses',
     >    'mol_co2 m-2 s-1',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (a10ancl4, cdummy)
      icount(3) = 1
      call writevar(filen,'a10ancl4',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, a10ancl4'
         stop 1
      end if
c
c 10-day average canopy photosynthesis rate - c3 grasses
c
      filen = fdir(1:lf)//'a10ancl3.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for 10-day average canopy photosynth. rate,'
     >    //' c3 grasses',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'a10ancl3',
     >    '10-day average canopy photosynth. rate, c3 grasses',
     >    'mol_co2 m-2 s-1',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (a10ancl3, cdummy)
      icount(3) = 1
      call writevar(filen,'a10ancl3',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, a10ancl3'
         stop 1
      end if
c
c 10-day average canopy scaling parameter - upper canopy
c
      filen = fdir(1:lf)//'a10scalparamu.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for 10-day average canopy scaling parameter,'
     >    //' upper canopy',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'a10scalparamu',
     >    '10-day average canopy scaling parameter, upper canopy',
     >    '_',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (a10scalparamu, cdummy)
      icount(3) = 1
      call writevar(filen,'a10scalparamu',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, a10scalparamu'
         stop 1
      end if
c
c 10-day average canopy scaling parameter - lower canopy
c
      filen = fdir(1:lf)//'a10scalparaml.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for 10-day average canopy scaling parameter,'
     >    //' lower canopy',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'a10scalparaml',
     >    '10-day average canopy scaling parameter, lower canopy',
     >    '_',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (a10scalparaml, cdummy)
      icount(3) = 1
      call writevar(filen,'a10scalparaml',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, a10scalparaml'
         stop 1
      end if
c
c 10-day average daylight - upper canopy
c
      filen = fdir(1:lf)//'a10daylightu.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for 10-day average daylight,'
     >    //' upper canopy',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'a10daylightu',
     >    '10-day average daylight, upper canopy',
     >    'W m-2',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (a10daylightu, cdummy)
      icount(3) = 1
      call writevar(filen,'a10daylightu',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, a10daylightu'
         stop 1
      end if
c
c 10-day average daylight - lower canopy
c
      filen = fdir(1:lf)//'a10daylightl.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for 10-day average daylight,'
     >    //' lower canopy',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'a10daylightl',
     >    '10-day average daylight, lower canopy',
     >    'W m-2',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (a10daylightl, cdummy)
      icount(3) = 1
      call writevar(filen,'a10daylightl',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, a10daylightl'
         stop 1
      end if
c
c drought-phenology trigger for trees
c
      filen = fdir(1:lf)//'dropu.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for drought-pheno. trigger for trees',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'dropu',
     >    'drought-pheno. trigger for trees',
     >    '_',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (dropu, cdummy)
      icount(3) = 1
      call writevar(filen,'dropu',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, dropu'
         stop 1
      end if
c
c drought-phenology trigger for shrubs
c
      filen = fdir(1:lf)//'dropls.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for drought-pheno. trigger for shrubs',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'dropls',
     >    'drought-pheno. trigger for shrubs',
     >    '_',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (dropls, cdummy)
      icount(3) = 1
      call writevar(filen,'dropls',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, dropls'
         stop 1
      end if
c
c drought-phenology trigger for c4 grasses
c
      filen = fdir(1:lf)//'dropl4.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for for drought-pheno. trigger for c4 grasses',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'dropl4',
     >    'drought-pheno. trigger for c4 grasses',
     >    '_',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (dropl4, cdummy)
      icount(3) = 1
      call writevar(filen,'dropl4',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, dropl4'
         stop 1
      end if
c
c drought-phenology trigger for c3 grasses
c
      filen = fdir(1:lf)//'dropl3.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for drought-pheno. trigger for c3 grasses',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'dropl3',
     >    'drought-pheno. trigger for c3 grasses',
     >    '_',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (dropl3, cdummy)
      icount(3) = 1
      call writevar(filen,'dropl3',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, dropl3'
         stop 1
      end if
c
      return
c
      end
c
c
c ---------------------------------------------------------------------
      subroutine wdaily (nday, iyear, iyear0)
c ---------------------------------------------------------------------
c
c writes out daily files
c
c ---------------------------------------------------------------------
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'comatm.h'
      include 'comsum.h'
      include 'comwork.h'
      include 'comveg.h'
      include 'argvs.h'
c
c Arguments
c
      integer nday,   ! number of days run since iyear0
     >        iyear,  ! this calendar year
     >        iyear0  ! very first year ever for this run sequence
c
c local variables
c
      integer  mstep,  ! this time step for netcdf file
     >        i,       ! loop indice
     >        idies,   ! netcdf file indice
     >        istat    ! netcdf error flag

      integer istart(4), icount(4) ! for writing vars
      real    pindex(npft)         ! index used for pfts and canopies
c
      character*11 cdate        ! date to use in history attribute in files
      character*10 tdate        ! character date for time step
      character*13 canopies(2)  ! canopy definitions
      character*21 tunits       ! time units
      character*80 dimnames(4)  ! names of dimensions for vars
      character*80 pftdef(npft) ! plant functional type defs (not currently used)
      character*80 filen        ! file name
c
      real ftime,              ! real form of nday
     >     tweight             ! number of days in daily average = 1.
c
c use workspace for local real variables
c
      equivalence (ftime,work(1)),(tweight,work(ndim4+1))
c
c ---------------------------------------------------------------------
c
      data istart / 1,1,1,1 /,
     >     icount / nlon,nlat,1,1 /
      icount(1) = nlonsub
      icount(2) = nlatsub
c
c     current time value, step, time weight
c
      ftime = nday
      mstep = nday
      tweight = 1.
c
c tdate is julian day (3 char) followed by iyear (4 char)
c
      tdate='0000000'//char(0)//char(0)//char(0)
      if (nday .lt. 10) then
         write(tdate(3:3),'(i1)') nday
      else if (nday .lt. 100) then
         write(tdate(2:3),'(i2)') nday
      else
         write(tdate(1:3),'(i3)') nday
      end if
      if (iyear .ge. 1000) then
         write(tdate(4:7),'(i4)') iyear
      else if (iyear .lt. 10) then
         write(tdate(7:7),'(i1)') iyear
      else if (iyear .lt. 100) then
         write(tdate(6:7),'(i2)') iyear
      else
         write(tdate(5:7),'(i3)') iyear
      end if
c
c first time only
c
      if (mstep.eq.1) then
c
c initialize snow layer indices, pft names, etc
c
         call simdate(cdate)
c
c time units is days since Dec 31 of the year before iyear0
c
         tunits = 'days since 0000-12-31'
         write(tunits(12:15),'(i4)') iyear0-1
c
         dimnames(1) = 'longitude'
         dimnames(2) = 'latitude'
c        dimnames(3) is set for each variable seperately
c        dimnames(4) not used in this subr for now, but define for future use
         dimnames(4) = 'time'

c
c define plant functional types, canopies with indices
c
         do i = 1, npft
           pindex(i) = i
         enddo
c and with characters
         pftdef(1) = 'trbrevtr - tropical broadleaf evergreen trees'
     >    //char(0)
         pftdef(2) = 
     >    'trbrdetr - tropical broadleaf drought-deciduous trees'
     >    //char(0)
         pftdef(3) = 
     >    'wtbrevtr - warm-temperate broadleaf evergreen trees'
     >    //char(0)
         pftdef(4) = 'tecoevtr - temperate conifer evergreen trees'
     >    //char(0)
         pftdef(5) = 
     >    'tebrdetr - temperate broadleaf cold-deciduous trees'//char(0)
         pftdef(6) = 'bocoevtr - boreal conifer evergreen trees'
     >    //char(0)
         pftdef(7) = 
     >    'bocodetr - boreal conifer cold-deciduous trees'//char(0)
         pftdef(8) = 
     >    'bobrdetr - boreal broadleaf cold-deciduous trees'//char(0)
         pftdef(9) = 'evsh - evergreen shrubs'//char(0)
         pftdef(10) = 'desh - deciduous shrubs'//char(0)
         pftdef(11) = 'c4gr - warm (c4) grasses'//char(0)
         pftdef(12) = 'c3gr - cool (c3) grasses'//char(0)

         canopies(1) = 'lower canopy'//char(0)
         canopies(2) = 'upper canopy'//char(0)
c
      end if
cc
cc dummy variable example, 3-d - copy & modify for new variable.
cc
c      filen = 'output/daily/dummyv.nc'
c      filen = out_daily_dummyv_
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'daily average dummyv',
c     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
c     >    'none',1,0.,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'dummyv','average dummyv',
c     >    'dummyvs-units',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      call vec2arr (addummyv, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'dummyv',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wdaily, dummyv'
c         stop 1
c      end if
c
c rainfall
c
c      filen = 'output/daily/rain.nc'
      filen = out_daily_rain_
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'daily average rainfall',
     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'rain','average rainfall',
     >    'mm/day',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (adrain, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'rain',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wdaily, rain'
         stop 1
      end if
c
c cloudiness
c
c     filen = 'output/daily/cloud.nc'
      filen = out_daily_cloud_
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'daily average cloudiness',
     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'cloud','average cloudiness',
     >    '%',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (cloud, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'cloud',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wdaily, cloud'
         stop 1
      end if
c
c rh
c
c     filen = 'output/daily/rh.nc'
      filen = out_daily_rh_
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'daily average relative humidity',
     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'rh','average rh',
     >    '%',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (adrh, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'rh',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wdaily, rh'
         stop 1
      end if
c
c snowfall
c
c     filen = 'output/daily/snow.nc'
      filen = out_daily_snow_
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'daily average snowfall',
     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'snow','average snowfall',
     >    'mm/day',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (adsnow, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'snow',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wdaily, snow'
         stop 1
      end if
cc
cc aet
cc
c     filen = 'output/daily/aet.nc'
      filen = out_daily_aet_
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'daily average aet',
     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'aet','average aet',
     >    'mm/day',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (adaet, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'aet',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wdaily, aet'
         stop 1
      end if
cc
cc trunoff
cc
c     filen = 'output/daily/trunoff.nc'
      filen = out_daily_trunoff_
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'daily average total runoff',
     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'trunoff','average total runoff',
     >    'mm/day',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (adtrunoff, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'trunoff',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wdaily, trunoff'
         stop 1
      end if
c
c srunoff
c
c     filen = 'output/daily/srunoff.nc'
      filen = out_daily_srunoff_
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'daily average surface runoff',
     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'srunoff','average surface runoff',
     >    'mm/day',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (adsrunoff, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'srunoff',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wdaily, srunoff'
         stop 1
      end if
c
c drainage
c
c     filen = 'output/daily/drainage.nc'
      filen = out_daily_drainage_
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'daily average drainage',
     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'drainage','average drainage',
     >    'mm/day',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (addrainage, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'drainage',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wdaily, drainage'
         stop 1
      end if
cc
cc wsoi
cc
c     filen = 'output/daily/wsoi.nc'
      filen = out_daily_wsoi_
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'daily average soil moisture',
     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'wsoi','average soil moisture',
     >    'fraction',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (adwsoi, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'wsoi',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wdaily, wsoi'
         stop 1
      end if
cc
cc wisoi
cc
c     filen = 'output/daily/wisoi.nc'
      filen = out_daily_wisoi_
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'daily average soil ice',
     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'wisoi','average soil ice',
     >    'fraction',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (adwisoi, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'wisoi',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wdaily, wisoi'
         stop 1
      end if
cc
cc snod
cc
c     filen = 'output/daily/snod.nc'
      filen = out_daily_snod_
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'daily average snow depth',
     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'snod','average snow depth',
     >    'meters',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (adsnod, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'snod',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wdaily, snod'
         stop 1
      end if
cc
cc snof
cc
c      filen = 'output/daily/snof.nc'
      filen = out_daily_snof_
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'daily average snow fraction',
c     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
c     >    'none',1,0.,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'snof','average snow fraction',
c     >    'fraction',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      call vec2arr (adsnof, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'snof',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wdaily, snof'
c         stop 1
c      end if
cc
cc co2ratio
cc
c      filen = 'output/daily/co2ratio.nc'
      filen = out_daily_co2ratio_
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'daily average ratio of root to total soil co2 flux',
c     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
c     >    'none',1,0.,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'co2ratio','average co2 ratio',
c     >    'fraction',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      call vec2arr (adco2ratio, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'co2ratio',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wdaily, co2ratio'
c         stop 1
c      end if
cc
cc co2mic
cc
c      filen = 'output/daily/co2mic.nc'
c     filen = out_daily_co2mic_
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'daily flux of carbon due to soil microbe co2 flux',
c     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
c     >    'none',1,0.,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'co2mic','soil microbe carbon flux',
c     >    'kg/m^2',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      call vec2arr (adco2mic, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'co2mic',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wdaily, co2mic'
c         stop 1
c      end if
c
cc templ 
cc
c      filen = 'output/daily/templ.nc'
c     filen = out_daily_templ_
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'index based on gdd for growth/sensecence',
c     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
c     >    'none',1,0.,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'templ','index based on gdd for growth/sensescence',
c     >    'dimensionless',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      call vec2arr (templ, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'templ',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wdaily, templ'
c         stop 1
c      end if
c
c bottom and top height of lower and upper canopies
c
c      filen = 'output/daily/zcanopy.nc'
c     filen = out_daily_zcanopy_
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'daily height of vegetation canopies',
c     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale,
c     >    'canopy','canopy','_',2,pindex,'up',
c     >    tunits,'gregorian',istat)
c add global attribute to define canopies with text, use netcdf low-level com
c         istat = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'canopy_def',
c     >    26+5,'1='//canopies(1)//' 2='//canopies(2))
c         dimnames(3) = 'canopy'
c         call inivar(idies,'zbot',
c     >    'bottom heights of lower and upper canopies',
c     >    'meters',4,dimnames,OCEAN,istat)
c         call inivar(idies,'ztop',
c     >    'top heights of lower and upper canopies',
c     >    'meters',4,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      do 20 k = 1, 2
c         call vec2arr (zbot(1,k), cdummy((k-1)*nlonsub*nlatsub + 1))
c 20   continue
c      istart(3) = 1
c      istart(4) = mstep
c      icount(3) = 2
c      call writevar(filen,'zbot',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wdaily, zbot'
c         stop 1
c      end if
c      do 25 k = 1, 2
c         call vec2arr (ztop(1,k), cdummy((k-1)*nlonsub*nlatsub + 1))
c 25   continue
c      istart(3) = 1
c      istart(4) = mstep
c      icount(3) = 2
c      call writevar(filen,'ztop',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wdaily, ztop'
c         stop 1
c      end if
c
c upper and lower canopy daily lai
c
c     filen = 'output/daily/laicanopy.nc'
      filen = out_daily_laicanopy_
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'daily lai-leaf area index of upper and lower vegetation canopies',
     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale,
     >    '','','',1,0.0,'',
     >    tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'laiu',
     >    'daily lai of upper canopy',
     >    'm2/m2',3,dimnames,OCEAN,istat)
       call inivar(idies,'lail',
     >    'daily lai of lower canopy',
     >    'm2/m2',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
c
      call vec2arr (lai(1,1), cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'laiu',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wdaily, laiu'
         stop 1
      end if
      call vec2arr (lai(1,2), cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'lail',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wdaily, lail'
         stop 1
      end if
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine wmonthly (nday, imonth, iyear, iyear0)
c ---------------------------------------------------------------------
c
c writes out monthly files
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'comsum.h'
      include 'comwork.h'
      include 'argvs.h'
c
c Arguments
c
      integer nday,    ! number of days run since iyear0
     >        imonth,  ! this month
     >        iyear,   ! this calendar year
     >        iyear0   ! very first year ever for this run
c
c Local variables
c
      integer mstep,   ! time step in netcdf file
     >        idies,   ! netcdf file indice
     >        istat    ! netcdf error flag
c
      integer istart(4), icount(4) ! for writing vars
c
      character*11 cdate       ! date to use in history attribute in files
      character*10 tdate       ! character date for time step
      character*21 tunits      ! units for time
      character*80 dimnames(4) ! names of dimensions for vars
      character*80 filen       ! file name
      character*3  chmon(12)   ! month abbreviations
c
      real ftime,              ! real form of nday
     >     tweight             ! number of days in monthly average
c
c use workspace for local real variables
c
      equivalence (ftime,work(1)),(tweight,work(ndim4+1))
c
c ---------------------------------------------------------------------
c
      data chmon / 'JAN','FEB','MAR','APR','MAY','JUN',
     >             'JUL','AUG','SEP','OCT','NOV','DEC' /
      data istart / 1,1,1,1 /,
     >     icount / nlon,nlat,1,1 /
      icount(1) = nlonsub
      icount(2) = nlatsub
c
c     current time value, step, time weight
c
      ftime = nday
      mstep = 12*(iyear-iyear0) + imonth
      tweight = ndaypm(imonth)
c
c tdate is this month (3 char), followed by this year (4 char)
c
      tdate=chmon(imonth)//'0000'//char(0)//char(0)//char(0)
      if (iyear .ge. 1000) then
         write(tdate(4:7),'(i4)') iyear
      else if (iyear .lt. 10) then
         write(tdate(7:7),'(i1)') iyear
      else if (iyear .lt. 100) then
         write(tdate(6:7),'(i2)') iyear
      else
         write(tdate(5:7),'(i3)') iyear
      end if
c
c first time only
c
      if (mstep.eq.1) then
c
c initialize snow layer indices, pft names, etc
c
         call simdate(cdate)
c
c time units is days since Dec 31 of the year before iyear0
c
         tunits = 'days since 0000-12-31'
         write(tunits(12:15),'(i4)') iyear0-1
c
         dimnames(1) = 'longitude'
         dimnames(2) = 'latitude'
c        dimnames(3) is set for each variable seperately
c        dimnames(4) not used in this subr for now, but define for future use
         dimnames(4) = 'time'
      end if
cc
cc dummy variable example, 3-d - copy & modify for new variable.
cc
c      filen = 'output/monthly/dummyv.nc'
      filen = out_monthly_dummyv_
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'monthly average dummyv',
c     >    'ibis wmonthly',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
c     >    'none',1,0.,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'dummyv','average dummyv',
c     >    'dummyvs-units',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      call vec2arr (amdummyv, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'dummyv',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wmonthly, dummyv'
c         stop 1
c      end if
c
c temperature
c
c     filen = 'output/monthly/temp.nc'
      filen = out_monthly_temp_
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'monthly average air temperature',
     >    'ibis wmonthly',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'temp','average air temperature',
     >    'C',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (amtemp, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'temp',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wmonthly, temp'
         stop 1
      end if
c
c rainfall
c
c     filen = 'output/monthly/rain.nc'
      filen = out_monthly_rain_
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'monthly average rainfall',
     >    'ibis wmonthly',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'rain','average rainfall',
     >    'mm/day',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (amrain, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'rain',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wmonthly, rain'
         stop 1
      end if
c
c cloudiness
c
c     filen = 'output/monthly/cloud.nc'
      filen = out_monthly_cloud_
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'monthly average cloudiness',
     >    'ibis wmonthly',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'cloud','average cloudiness',
     >    '%',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (amcloud, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'cloud',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wmonthly, cloud'
         stop 1
      end if
c
c rh
c
c     filen = 'output/monthly/rh.nc'
      filen = out_monthly_rh_
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'monthly average relative humidity',
     >    'ibis wmonthly',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'rh','average rh',
     >    '%',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (amrh, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'rh',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wmonthly, rh'
         stop 1
      end if
c
c snowfall
c
c     filen = 'output/monthly/snow.nc'
      filen = out_monthly_snow_
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'monthly average snowfall',
     >    'ibis wmonthly',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'snow','average snowfall',
     >    'mm/day',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (amsnow, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'snow',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wmonthly, snow'
         stop 1
      end if
c
c specific humidity
c
c     filen = 'output/monthly/qa.nc'
      filen = out_monthly_qa_
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'monthly average specific humidity',
     >    'ibis wmonthly',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'qa','average specific humidity',
     >    'kg-h2o/kg-air',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (amqa, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'qa',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wmonthly, qa'
         stop 1
      end if
c
c evapotranspiration
c
c     filen = 'output/monthly/aet.nc'
      filen = out_monthly_aet_
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'monthly average aet',
     >    'ibis wmonthly',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'aet','average evapotranspiration',
     >    'mm/day',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (amaet, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'aet',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wmonthly, aet'
         stop 1
      end if
c
c trunoff, srunoff, drainage
c
c     filen = 'output/monthly/runoff.nc'
      filen = out_monthly_runoff_
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'monthly average total runoff',
     >    'ibis wmonthly',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'trunoff','average total runoff',
     >    'mm/day',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'srunoff','average surface runoff',
     >    'mm/day',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'drainage','average drainage',
     >    'mm/day',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (amtrunoff, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'trunoff',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wmonthly, trunoff'
         stop 1
      end if
      call vec2arr (amsrunoff, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'srunoff',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wmonthly, srunoff'
         stop 1
      end if
      call vec2arr (amdrainage, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'drainage',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wmonthly, drainage'
         stop 1
      end if
c
c soil temperature
c
c     filen = 'output/monthly/tsoi.nc'
      filen = out_monthly_tsoi_
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'monthly average soil temperature',
     >    'ibis wmonthly',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'tsoi','average soil temperature',
     >    'degC',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (amtsoi, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'tsoi',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wmonthly, tsoi'
         stop 1
      end if
c
c soil moisture, ice, volumetric water content, plant available water
c
c     filen = 'output/monthly/wsoi.nc'
      filen = out_monthly_wsoi_
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'monthly average soil moisture',
     >    'ibis wmonthly',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'wsoi','average soil moisture',
     >    'fraction',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'wisoi','average soil ice',
     >    'fraction',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'vwc','average volumetric water content',
     >    'fraction',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'awc',
     >    'average plant available water content',
     >    'cm',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (amwsoi, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'wsoi',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wmonthly, wsoi'
         stop 1
      end if
      call vec2arr (amwisoi, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'wisoi',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wmonthly, wisoi'
         stop 1
      end if
      call vec2arr (amvwc, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'vwc',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wmonthly, vwc'
         stop 1
      end if
      call vec2arr (amawc, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'awc',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wmonthly, awc'
         stop 1
      end if
c
c snow depth
c
c     filen = 'output/monthly/snod.nc'
      filen = out_monthly_snod_
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'monthly average snow depth',
     >    'ibis wmonthly',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'snod','average snow depth',
     >    'meters',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (amsnod, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'snod',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wmonthly, snod'
         stop 1
      end if
c
c snow fraction
c
c     filen = 'output/monthly/snof.nc'
      filen = out_monthly_snof_
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'monthly average snow fraction',
     >    'ibis wmonthly',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'snof','average snow fraction',
     >    'm^2/m^3',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (amsnof, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'snof',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wmonthly, snof'
         stop 1
      end if
c
c solar radiation
c
c     filen = 'output/monthly/solar.nc'
      filen = out_monthly_solar_
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'monthly average incident solar radiation',
     >    'ibis wmonthly',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'solar','average incident solar radiation',
     >    'W/m^2',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (amsolar, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'solar',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wmonthly, solar'
         stop 1
      end if
c
c albedo
c
c     filen = 'output/monthly/albedo.nc'
      filen = out_monthly_albedo_
c     if (mstep .eq. 1) then
c        call inifile(idies,filen,
c    >    'monthly average albedo',
c    >    'ibis wmonthly',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
c    >    'none',1,0.,'',tunits,'gregorian',istat)
c        dimnames(3) = 'time'
c        call inivar(idies,'albedo','average albedo',
c    >    'fraction',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c     end if
c     call vec2arr (amalbedo, cdummy)
c     istart(3) = mstep
c     icount(3) = 1
c     call writevar(filen,'albedo',istart,icount,cdummy,ftime,
c    > tweight,tdate,istat)
c     if (istat .ne. 0) then
c        write(*,*) 'ERROR in wmonthly, albedo'
c        stop 1
c     end if
c
c downward and upward infrared radiation
c
c     filen = 'output/monthly/ir.nc'
      filen = out_monthly_ir_
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'monthly average infrared radiation',
     >    'ibis wmonthly',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'irdown','average downward IR',
     >    'W/m^2',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'irup','average upward IR',
     >    'W/m^2',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (amirdown, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'irdown',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wmonthly, irdown'
         stop 1
      end if
      call vec2arr (amirup, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'irup',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wmonthly, irup'
         stop 1
      end if
c
c sensible heat flux
c
c     filen = 'output/monthly/sens.nc'
      filen = out_monthly_sens_
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'monthly average sensible heat flux',
     >    'ibis wmonthly',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'sens','average sensible heat flux',
     >    'W/m^2',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (amsens, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'sens',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wmonthly, sens'
         stop 1
      end if
c
c latent heat flux
c
c     filen = 'output/monthly/latent.nc'
      filen = out_monthly_latent_
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'monthly average latent heat flux',
     >    'ibis wmonthly',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'latent','average latent heat flux',
     >    'W/m^2',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (amlatent, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'latent',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wmonthly, latent'
         stop 1
      end if
c
c leaf area index upper and lower
c
      filen = 'output/monthly/lai.nc'
          filen = out_monthly_lai_
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >   'monthly average leaf area index for upper and lower canopies',
     >    'ibis wmonthly',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'laiu','average lai upper canopy',
     >    'fraction',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'lail','average lai lower canopy',
     >    'fraction',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (amlaiu, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'laiu',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wmonthly, laiu'
         stop 1
      end if
      call vec2arr (amlail, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'lail',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wmonthly, lail'
         stop 1
      end if
c
c total net primary productivity
c
      filen = 'output/monthly/npptot.nc'
          filen = out_monthly_npptot_
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'monthly total net primary productivity of carbon',
     >    'ibis wmonthly',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'npptot','total npp of carbon',
     >    'kg m-2 month-1',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (amnpptot, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'npptot',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wmonthly, npptot'
         stop 1
      end if
c
c co2ratio
c
      filen = 'output/monthly/co2ratio.nc'
          filen = out_monthly_co2ratio_
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'monthly ratio of root respiration to total soil respiration',
     >    'ibis wmonthly',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'co2ratio',
     >    'ratio of root to soil respiration',
     >    'fraction',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (amco2ratio, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'co2ratio',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wmonthly, co2ratio'
         stop 1
      end if
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine wyearly (nday,iyear,iyear0)
c ---------------------------------------------------------------------
c
c writes out yearly files
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'comsum.h'
      include 'comveg.h'
      include 'comwork.h'
      include 'argvs.h'
c
c Arguments
c
      integer nday,   ! number of days run since iyear0
     >        iyear,  ! this calendar year
     >        iyear0  ! very first year ever for this run
c
c local variables
c
      integer mstep,  ! time for this year step
     >        idies,  ! netcdf file indice
     >        istat,  ! netcdf error flag
     >        i,k     ! loop indices        
c
      integer istart(4), icount(4) ! for writing vars
      real    pindex(npft)         ! index used for pfts and canopies
c
      character*11 cdate        ! date to use in history attribute in files
      character*10 tdate        ! character date for time step
      character*13 canopies(2)  ! canopy definitions
      character*21 tunits       ! time units
      character*80 dimnames(4)  ! names of dimensions for vars
      character*80 pftdef(npft) ! plant functional type defs (not currently used)
      character*80 filen        ! file name
c
      real ftime,               ! real form of nday
     >     tweight              ! number of days in yearly average
c
c External
c 
      integer NF_PUT_ATT_TEXT,   ! netcdf function
     >        NF_GLOBAL          ! netcdf function
c
c use workspace for local real variables
c
      equivalence (ftime,work(1)),(tweight,work(ndim4+1))
c
c ---------------------------------------------------------------------
c
      data istart / 1,1,1,1 /,
     >     icount / nlon,nlat,1,1 /
      icount(1) = nlonsub
      icount(2) = nlatsub
c
c current time value and step: make ftime Jan 1 of this year 
c instead of Dec 31
c
      ftime = nday - ndaypy + 1
      tweight = ndaypy
      mstep = iyear - iyear0 + 1
c
c tdate is ANN (3 char) followed by this year (4 char)
c
      tdate='ANN0000'//char(0)//char(0)//char(0)
      if (iyear .ge. 1000) then
         write(tdate(4:7),'(i4)') iyear
      else if (iyear .lt. 10) then
         write(tdate(7:7),'(i1)') iyear
      else if (iyear .lt. 100) then
         write(tdate(6:7),'(i2)') iyear
      else
         write(tdate(5:7),'(i3)') iyear
      end if
c
c first time only
c
      if (mstep .eq. 1) then
c
         call simdate(cdate)
c
c time units is days since Dec 31 of the year before iyear0
c
         tunits = 'days since 0000-12-31'
         write(tunits(12:15),'(i4)') iyear0-1
c
c dimension names
c
         dimnames(1) = 'longitude'
         dimnames(2) = 'latitude'
c        dimnames(3) is set for each variable seperately
         dimnames(4) = 'time'
c
c define plant functional types, canopies with indices
c
         do i = 1, npft
           pindex(i) = i
         enddo
c and with characters
         pftdef(1) = 'trbrevtr - tropical broadleaf evergreen trees'
     >    //char(0)
         pftdef(2) = 
     >    'trbrdetr - tropical broadleaf drought-deciduous trees'
     >    //char(0)
         pftdef(3) = 
     >    'wtbrevtr - warm-temperate broadleaf evergreen trees'
     >    //char(0)
         pftdef(4) = 'tecoevtr - temperate conifer evergreen trees'
     >    //char(0)
         pftdef(5) = 
     >    'tebrdetr - temperate broadleaf cold-deciduous trees'//char(0)
         pftdef(6) = 'bocoevtr - boreal conifer evergreen trees'
     >    //char(0)
         pftdef(7) = 
     >    'bocodetr - boreal conifer cold-deciduous trees'//char(0)
         pftdef(8) = 
     >    'bobrdetr - boreal broadleaf cold-deciduous trees'//char(0)
         pftdef(9) = 'evsh - evergreen shrubs'//char(0)
         pftdef(10) = 'desh - deciduous shrubs'//char(0)
         pftdef(11) = 'c4gr - warm (c4) grasses'//char(0)
         pftdef(12) = 'c3gr - cool (c3) grasses'//char(0)

         canopies(1) = 'lower canopy'//char(0)
         canopies(2) = 'upper canopy'//char(0)
c
      end if
cc
cc dummy variable example, 4-d var whose 3rd dim is a character dim (pft)
cc copy and modify for a new variable
cc
c      filen = 'output/yearly/dummyv.nc'
       filen = out_yearly_dummyv_
c      if (mstep .eq. 1) then
c         call inifilec(idies,filen,
c     >    'annual dummyv',
c     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,latscale,
c     >    'pft','plant fuctional type','_',npft,80,pftdef,
c     >    tunits,'gregorian',istat)
c         dimnames(3) = 'pft'
c         call inivar(idies,'dummyv','dummyv for each pft',
c     >    'dummyvs-units',4,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      do 5 k = 1, npft
c         call vec2arr (aydummyv(1,k), cdummy((k-1)*nlonsub*nlatsub + 1))
c 5    continue
c      istart(3) = 1
c      istart(4) = mstep
c      icount(3) = npft
c      call writevar(filen,'dummyv',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wyearly, dummyv'
c         stop 1
c      end if
cc
cc dummy variable example, 3-d
cc
c      filen = 'output/yearly/dummyv.nc'
       filen = out_yearly_dummyv_
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'dummyv','annual dummyv',cdate,nlonsub,lonscale,
c     >    nlatsub,latscale,'','none','none',1,0.,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'dummyv','annual dummyv',
c     >    'dummyvs-units',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      call vec2arr (aydummyv, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'dummyv',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wyearly, dummyv'
c         stop 1
c      end if
c
c net primary productivity, by pft and total
c
      filen = 'output/yearly/npp.nc'
           filen = out_yearly_npp_
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'annual net primary productivity of carbon',
     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,latscale,
     >    'pft','plant fuctional type','_',npft,pindex,'',
     >    tunits,'gregorian',istat)
c add global attribute to define pfts with text, use netcdf low-level command
         istat = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'pft_definition',
     >    npft*80,pftdef)
         dimnames(3) = 'pft'
         call inivar(idies,'npp','npp of carbon for each pft',
     >    'kg m-2 year-1',4,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'npptot','total npp',
     >    'kg m-2 year-1',3,dimnames,OCEAN,istat)
         call inivar(idies,'anpptot','total above-ground npp',
     >    'kg m-2 year-1',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      do 5 k = 1, npft
         call vec2arr (aynpp(1,k), cdummy((k-1)*nlonsub*nlatsub + 1))
 5    continue
      istart(3) = 1
      istart(4) = mstep
      icount(3) = npft
      call writevar(filen,'npp',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, npp'
         stop 1
      end if
      call vec2arr (aynpptot, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'npptot',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, npptot'
         stop 1
      end if
      call vec2arr (ayanpptot, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'anpptot',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, anpptot'
         stop 1
      end if
c
c evapotranspiration
c
      filen = 'output/yearly/aet.nc'
           filen = out_yearly_aet_
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'annual average evapotranspiration',
     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'aet','average evapotranspiration',
     >    'mm/year',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (ayaet, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'aet',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, aet'
         stop 1
      end if
c
c trunoff, srunoff, drainage, rratio, tratio
c
      filen = 'output/yearly/runoff.nc'
           filen = out_yearly_runoff_
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'annual runoff',
     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'trunoff','total runoff',
     >    'mm/year',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'srunoff','surface runoff',
     >    'mm/year',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'drainage','drainage',
     >    'mm/year',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'rratio','average runoff ratio',
     >    'fraction',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'tratio','average transpiration ratio',
     >    'fraction',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (aytrunoff, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'trunoff',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, trunoff'
         stop 1
      end if
      call vec2arr (aysrunoff, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'srunoff',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, srunoff'
         stop 1
      end if
      call vec2arr (aydrainage, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'drainage',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, drainage'
         stop 1
      end if
      call vec2arr (ayrratio, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'rratio',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, rratio'
         stop 1
      end if
      call vec2arr (aytratio, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'tratio',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, tratio'
         stop 1
      end if
c
c soil moisture, soil ice, volumetric water content, plant available water
c
      filen = 'output/yearly/wsoi.nc'
           filen = out_yearly_wsoi_
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'annual average soil moisture',
     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'wsoi','average soil moisture',
     >    'fraction',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'wisoi','average soil ice',
     >    'fraction',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'vwc','average volumetric water content',
     >    'fraction',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'awc','average volumetric water content',
     >    'cm',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (aywsoi, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'wsoi',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, wsoi'
         stop 1
      end if
      call vec2arr (aywisoi, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'wisoi',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, wisoi'
         stop 1
      end if
      call vec2arr (ayvwc, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'vwc',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, vwc'
         stop 1
      endif
      call vec2arr (ayawc, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'awc',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, awc'
         stop 1
      endif
c
c soil temperature
c
      filen = 'output/yearly/tsoi.nc'
           filen = out_yearly_tsoi_
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'annual average soil temperature',
     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'tsoi','average soil temperature',
     >    'degrees C',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (aytsoi, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'tsoi',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, tsoi'
         stop 1
      endif
cc
cc solar radiation
cc
c      filen = 'output/yearly/solar.nc'
       filen = out_yearly_solar_
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'annual average solar incident radiation',
c     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
c     >    'none',1,0.,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'solar','average solar radiation',
c     >    'W/m^2',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      call vec2arr (aysolar, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'solar',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wyearly, solar'
c         stop 1
c      end if
cc
cc albedo
cc
c     filen = 'output/yearly/albedo.nc'
       filen = out_yearly_albedo_
c     if (mstep .eq. 1) then
c        call inifile(idies,filen,
c    >    'annual average albedo',
c    >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
c    >    'none',1,0.,'',tunits,'gregorian',istat)
c        dimnames(3) = 'time'
c        call inivar(idies,'albedo','average albedo',
c    >    'fraction',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c     end if
c     call vec2arr (ayalbedo, cdummy)
c     istart(3) = mstep
c     icount(3) = 1
c     call writevar(filen,'albedo',istart,icount,cdummy,ftime,
c    > tweight,tdate,istat)
c     if (istat .ne. 0) then
c        write(*,*) 'ERROR in wyearly, albedo'
c        stop 1
c     end if
cc
cc upward and downward infrared radiation
cc
c      filen = 'output/yearly/ir.nc'
c      filen = out_yearly_ir_
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'annual average solar infrared radiation',
c     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
c     >    'none',1,0.,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'irdown','average downward ir',
c     >    'W/m^2',3,dimnames,OCEAN,istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'irup','average upward ir',
c     >    'W/m^2',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      call vec2arr (ayirdown, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'irdown',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wyearly, irdown'
c         stop 1
c      end if
c      call vec2arr (ayirup, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'irup',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wyearly, irup'
c         stop 1
c      end if
c
c sensible heat flux
c
      filen = 'output/yearly/sens.nc'
           filen = out_yearly_sens_
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'annual average sensible heat flux',
     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'sens','average sensible heat flux',
     >    'W/m^2',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (aysens, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'sens',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, sens'
         stop 1
      end if
cc
cc latent heat flux
cc
c      filen = 'output/yearly/latent.nc'
       filen = out_yearly_latent_
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'annual average latent heat flux',
c     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
c     >    'none',1,0.,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'latent','average latent heat flux',
c     >    'W/m^2',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      call vec2arr (aylatent, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'latent',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wyearly, latent'
c         stop 1
c      end if
c
c lai, by pft, total upper canopy, total lower canopy
c
      filen = 'output/yearly/plai.nc'
           filen = out_yearly_plai_
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'annual leaf area index',
     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,latscale,
     >    'pft','plant fuctional type','_',npft,pindex,'',
     >    tunits,'gregorian',istat)
c add global attribute to define pfts with text, use netcdf low-level command
         istat = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'pft_definition',
     >    npft*80,pftdef)
         dimnames(3) = 'pft'
         call inivar(idies,'plai','leaf area index for each pft',
     >    'fraction',4,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'totlaiu','total lai for upper canopy',
     >    'fraction',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'totlail','total lai for lower canopy',
     >    'fraction',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      do 10 k = 1, npft
         call vec2arr (plai(1,k), cdummy((k-1)*nlonsub*nlatsub + 1))
 10   continue
      istart(3) = 1
      istart(4) = mstep
      icount(3) = npft
      call writevar(filen,'plai',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, plai'
         stop 1
      end if
      call vec2arr (totlaiu, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'totlaiu',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, totlaiu'
         stop 1
      end if
      call vec2arr (totlail, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'totlail',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, totlail'
         stop 1
      end if
c
c biomass, by pft, upper canopy, lower canopy
c
      filen = 'output/yearly/biomass.nc'
           filen = out_yearly_biomass_
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'annual biomass of carbon',
     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,latscale,
     >    'pft','plant fuctional type','_',npft,pindex,'',
     >    tunits,'gregorian',istat)
c add global attribute to define pfts with text, use netcdf low-level command
         istat = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'pft_definition',
     >    npft*80,pftdef)
         dimnames(3) = 'pft'
         call inivar(idies,'biomass','biomass for each pft',
     >    'kg/m^2',4,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'totbiou','total biomass for upper canopy',
     >    'kg/m^2',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'totbiol','total biomass for lower canopy',
     >    'kg/m^2',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      do 15 k = 1, npft
         call vec2arr (biomass(1,k), cdummy((k-1)*nlonsub*nlatsub + 1))
 15   continue
      istart(3) = 1
      istart(4) = mstep
      icount(3) = npft
      call writevar(filen,'biomass',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, biomass'
         stop 1
      end if
      call vec2arr (totbiou, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'totbiou',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, totbiou'
         stop 1
      end if
      call vec2arr (totbiol, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'totbiol',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, totbiol'
         stop 1
      end if
c
c soil carbon: rootbio, totalit, totrlit, totcsoi, totcmic
c
      filen = 'output/yearly/csoi.nc'
           filen = out_yearly_csoi_
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'annual total soil carbon',
     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'rootbio','total live root biomass carbon',
     >    'kg/m^2',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'totalit','total above ground litter carbon',
     >    'kg/m^2',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'totrlit','total below ground litter carbon',
     >    'kg/m^2',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'totcsoi','total soil carbon w/o litter',
     >    'kg/m^2',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'totcmic','total microbial carbon',
     >    'kg/m^2',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (ayrootbio, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'rootbio',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, rootbio'
         stop 1
      end if
      call vec2arr (ayalit, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'totalit',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, totalit'
         stop 1
      end if
      call vec2arr (ayblit, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'totrlit',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, totrlit'
         stop 1
      end if
      call vec2arr (aycsoi, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'totcsoi',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, totcsoi'
         stop 1
      end if
      call vec2arr (aycmic, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'totcmic',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, totcmic'
         stop 1
      end if
c
c soil nitrogen: totanlit, totrnlit, totnsoi, nmintot
c
      filen = 'output/yearly/nsoi.nc'
           filen = out_yearly_nsoi_
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'annual total soil nitrogen',
     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'totanlit',
     >    'total above ground litter nitrogen',
     >    'kg/m^2',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'totrnlit',
     >    'total below ground litter nitrogen',
     >    'kg/m^2',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'totnsoi','total soil nitrogen w/o litter',
     >    'kg/m^2',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'nmintot','total nitrogen mineralization',
     >    'kg/m^2',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (ayanlit, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'totanlit',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, totanlit'
         stop 1
      end if
      call vec2arr (aybnlit, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'totrnlit',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, totrnlit'
         stop 1
      end if
      call vec2arr (aynsoi, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'totnsoi',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, totnsoi'
         stop 1
      end if
      call vec2arr (aynmintot, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'nmintot',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, nmintot'
         stop 1
      end if
cc
cc total litter
cc
c      filen = 'output/yearly/totlit.nc'
c      filen = out_yearly_totlit_
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'annual total litter carbon',
c     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
c     >    'none',1,0.,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'totlit','total litter carbon',
c     >    'kg/m^2',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      call vec2arr (totlit, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'totlit',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wyearly, totlit'
c         stop 1
c      end if
cc
cc total wood litter
cc
c      filen = 'output/yearly/clitw.nc'
       filen = out_yearly_clitw_
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'annual total wood litter carbon',
c     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
c     >    'none',1,0.,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'clitw','total wood litter carbon',
c     >    'kg/m^2',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      call vec2arr (clitw, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'clitw',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wyearly, clitw'
c         stop 1
c      end if
cc
cc total litterfall
cc
c      filen = 'output/yearly/totfall.nc'
       filen = out_yearly_totfall_
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'annual total litterfall carbon',
c     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
c     >    'none',1,0.,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'totfall','total litterfall',
c     >    'kg/m^2',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      call vec2arr (totfall, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'totfall',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wyearly, totfall'
c         stop 1
c      end if
cc
cc total soil carbon in slow pool
cc
c      filen = 'output/yearly/csoislo.nc'
       filen = out_yearly_csoislo_
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'annual total soil carbon in slow pool',
c     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
c     >    'none',1,0.,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'csoislo','total soil carbon in slow pool',
c     >    'kg/m^2',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      call vec2arr (csoislo, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'csoislo',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wyearly, csoislo'
c         stop 1
c      end if
cc
cc total soil carbon in passive pool
cc
c      filen = 'output/yearly/csoipas.nc'
       filen = out_yearly_csoipas_
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'annual total soil carbon in pasive pool',
c     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
c     >    'none',1,0.,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'csoipas','total soil carbon in passive pool',
c     >    'kg/m^2',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      call vec2arr (csoipas, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'csoipas',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wyearly, csoipas'
c         stop 1
c      end if
c
c co2 carbon exchange: net ecosystem, microbial resp, root resp, soil resp
c
      filen = 'output/yearly/co2fluxes.nc'
           filen = out_yearly_co2fluxes_
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'annual total carbon from exchange of co2',
     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'neetot',
     >    'total net ecosystem echange carbon',
     >    'kg/m^2',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'co2mic','total microbe respiration carbon',
     >    'kg/m^2',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'co2root','total root respiration carbon',
     >    'kg/m^2',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'co2soi','total soil respiration carbon',
     >    'kg/m^2',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (ayneetot, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'neetot',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, neetot'
         stop 1
      end if
      call vec2arr (ayco2mic, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'co2mic',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, co2mic'
         stop 1
      end if
      call vec2arr (ayco2root, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'co2root',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, co2root'
         stop 1
      end if
      call vec2arr (ayco2soi, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'co2soi',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, co2soi'
         stop 1
      end if
c
c fire disturbance regime
c
      filen = 'output/yearly/disturbf.nc'
           filen = out_yearly_disturbf_
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'annual fire disturbance regime',
     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'disturbf','fire disturbance regime',
     >    'fraction/year',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (disturbf, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'disturbf',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, disturbf'
         stop 1
      end if
c
c vegetation type
c
      filen = 'output/yearly/vegtype0.nc'
      filen = out_yearly_vegtype0_
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'annual vegetation type - ibis classification',
     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'vegtype0','vegetation type',
     >    '_',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (vegtype0, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'vegtype0',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, vegtype0'
         stop 1
      end if
c
c fractional cover of upper and lower canopies
c
      filen = 'output/yearly/fcover.nc'
      filen = out_yearly_fcover_
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'annual fractional cover of canopies','ibis wyearly',
     >    cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
     >    'none',1,0.,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'fu','fractional cover of upper canopy',
     >    'fraction',3,dimnames,OCEAN,istat)
         call inivar(idies,'fl','fractional cover of lower canopy',
     >    'fraction',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (fu, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'fu',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, fu'
         stop 1
      end if
      call vec2arr (fl, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'fl',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, fl'
         stop 1
      end if
cc
cc sapwood fraction
cc
c      filen = 'output/yearly/sapfrac.nc'
c      filen = out_yearly_sapfrac
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'annual sapwood fraction',
c     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
c     >    'none',1,0.,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'sapfrac','sapwood fraction',
c     >    'fraction',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      call vec2arr (sapfrac, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'sapfrac',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wyearly, sapfrac'
c         stop 1
c      end if
c
c bottom and top of vegetation canopies
c
      filen = 'output/yearly/zcanopy.nc'
      filen = out_yearly_zcanopy_
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'annual height of vegetation canopies',
     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,latscale,
     >    'canopy','canopy','_',2,pindex,'up',
     >    tunits,'gregorian',istat)
c add global attribute to define canopies with text, use netcdf low-level com
         istat = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'canopy_def',
     >    26+5,'1='//canopies(1)//' 2='//canopies(2))
         dimnames(3) = 'canopy'
         call inivar(idies,'zbot',
     >    'bottom heights of lower and upper canopies',
     >    'meters',4,dimnames,OCEAN,istat)
         call inivar(idies,'ztop',
     >    'top heights of lower and upper canopies',
     >    'meters',4,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      do 20 k = 1, 2
         call vec2arr (zbot(1,k), cdummy((k-1)*nlonsub*nlatsub + 1))
 20   continue
      istart(3) = 1
      istart(4) = mstep
      icount(3) = 2
      call writevar(filen,'zbot',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, zbot'
         stop 1
      end if
      do 25 k = 1, 2
         call vec2arr (ztop(1,k), cdummy((k-1)*nlonsub*nlatsub + 1))
 25   continue
      istart(3) = 1
      istart(4) = mstep
      icount(3) = 2
      call writevar(filen,'ztop',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, ztop'
         stop 1
      end if
c
c existence of pfts
c
      filen = 'output/yearly/exist.nc'
      filen = out_yearly_exist_
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'annual existence of each plant functional type',
     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,latscale,
     >    'pft','plant fuctional type','_',npft,pindex,'',
     >    tunits,'gregorian',istat)
c add global attribute to define pfts with text, use netcdf low-level command
         istat = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'pft_definition',
     >    npft*80,pftdef)
         dimnames(3) = 'pft'
         call inivar(idies,'exist','existence for each pft',
     >    '_',4,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      do 30 k = 1, npft
         call vec2arr (exist(1,k), cdummy((k-1)*nlonsub*nlatsub + 1))
 30   continue
      istart(3) = 1
      istart(4) = mstep
      icount(3) = npft
      call writevar(filen,'exist',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, exist'
         stop 1
      end if
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine readit (isimveg,snorth,ssouth,swest,seast,iwest,jnorth)
c ---------------------------------------------------------------------
c
c reads in initialization files and initializes some fields
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'combcs.h'
      include 'comsoi.h'
      include 'comveg.h'
      include 'comwork.h'
      include 'argvs.h'
c
c Arguments
c
      integer isimveg,  ! dynamic vegetation (1) or static (0)
     >        iwest,    !
     >        jnorth    !       
c
      real snorth, ssouth, swest, seast
c
c Local variables
c

      integer  istat,    ! netcdf error flag
     > i, j, ntime,      ! loop indices
     > jj, ii, 
     > ndim,             ! number of dimensions
     > nlpoints          ! number of land points
c
      real xlat
c
      integer jsouth, ieast
      real xmask(nlon,nlat)
      equivalence ( xmask(1,1), work(1) )
c
      character*80 filen
c
      integer istart(4), icount(4) ! for reading vars
c
c ---------------------------------------------------------------------
c
      data istart / 1,1,1,1 /, icount / nlon,nlat,1,1 /
c
c 2-d surface and vegetation arrays
c
      icount(3) = 1
c
c land mask, latitudes, and longitudes
c
      filen = 'input/surta.nc'
      filen = surta_
      aname = 'surta'
      call readvar(filen,aname,'level',istart,icount,
     > xmask,lonscale,latscale,cdummy(1),cdummy(ndim4),istat)
      if (istat.lt.0) then
         write (*,9000)
         print *, 'while reading surta'
         stop 1
      end if
c
      if (abs(abs(lonscale(1)-lonscale(2))-xres).gt.0.001 .or.
     >    abs(abs(latscale(1)-latscale(2))-yres).gt.0.001) then
         write (*,9000)
         write(*,*) 'resolution mismatch!'
         write(*,*) 'xres, yres in compar.h = ', xres, yres
         write(*,*) 'xres, yres from input/surta.nc = ',
     >   abs(lonscale(1)-lonscale(2)), abs(latscale(1)-latscale(2))
         stop 1
      end if
c
c subset the grid if not default whole grid
c
       if (snorth.lt.latscale(1) .or. ssouth.gt.latscale(nlat) .or.
     >     swest.gt.lonscale(1) .or.  seast.lt.lonscale(nlon)) then
        jnorth = 0
        if (snorth .lt. (latscale(nlat)+latscale(nlat-1))/2.) then
          jnorth = nlat
        else if (snorth .ge. (latscale(1)+latscale(2))/2.) then
          jnorth = 1
        else
          do 1 j = nlat-1,1,-1
            if (snorth .ge. (latscale(j)+latscale(j+1))/2.) jnorth = j
 1        continue
        end if
        jsouth = 0
        if (ssouth .lt. (latscale(nlat)+latscale(nlat-1))/2.) then
          jsouth = nlat
        else if (ssouth .ge. (latscale(1)+latscale(2))/2.) then
          jsouth = 1
        else
          do 2 j = nlat-1,1,-1
            if (ssouth .ge. (latscale(j)+latscale(j+1))/2.) jsouth = j
 2        continue
        end if
c
        iwest = 0
        if (swest .lt. (lonscale(1)+lonscale(2))/2.) then
          iwest = 1
        else if (swest .ge. (lonscale(nlon)+lonscale(nlon-1))/2.) then
          iwest = nlon
        else
          do 3 i = 2, nlon
            if(swest .ge. (lonscale(i)+lonscale(i-1))/2.) iwest=i
 3        continue
        end if
        ieast = 0
        if (seast .lt. (lonscale(1)+lonscale(2))/2.) then
          ieast = 1
        else if (seast .ge. (lonscale(nlon)+lonscale(nlon-1))/2.) then
          ieast = nlon
        else
          do 4 i = 2, nlon
            if(seast .ge. (lonscale(i)+lonscale(i-1))/2.) ieast=i
 4        continue
        end if
        nlonsub = ieast - iwest + 1
        nlatsub = jsouth - jnorth + 1
        istart(1) = iwest
        icount(1) = nlonsub
        istart(2) = jnorth
        icount(2) = nlatsub
      else
        iwest = 1
        ieast = nlon
        jnorth = 1
        jsouth = nlat
        nlonsub = nlon
        nlatsub = nlat
      end if
cc      print *, iwest, ieast, jnorth, jsouth
cc      print *, swest, seast, snorth, ssouth
c
c
c initialize lonindex, latindex for use in arr2vec, vec2arr, etc.
c and calculate the approximate the area of the land gridcells
c
      nlpoints = 0
c
c here, i/j refers to entire grid (1 to nlon/nlat), 
c ii/jj refers to subgrid (1 to nlonsub/nlatsub)
c
      do 5 j = jnorth, jsouth
c
        jj = j - jnorth + 1
        do 6 i = iwest, ieast
c
          ii = i - iwest + 1
          lmask(ii,jj) = nint(xmask(i,j))
c
          if (lmask(ii,jj).eq.1) then
c
            nlpoints = nlpoints + 1
            lonindex(nlpoints) = ii
            latindex(nlpoints) = jj
            xlat = latscale(j) * pi / 180.0
            garea(nlpoints) = yres * 111400.0 * xres * 111400.0 *
     >                        cos(xlat)
c
          end if
c
 6      continue
 5    continue
c
cc      print *, jnorth, jsouth, iwest, ieast
cc      print *, latscale(jnorth),latscale(jsouth),lonscale(iwest),lonscale(ieast)
cc      print *, ((lmask(i,j),i=1,nlonsub),j=1,nlatsub)
cc      print *, istart
cc      print *, icount
c
      do 7 j = jnorth, jsouth
        jj = j - jnorth + 1
        latscale(jj) = latscale(j)
 7    continue
c
      do 8 i = iwest, ieast
        ii = i - iwest + 1
        lonscale(ii) = lonscale(i)
 8    continue
c
cc      print *, (lonscale(i),i=1,nlonsub)
cc      print *, (latscale(i),i=1,nlatsub)
cc      print *, lonindex
cc      print *, latindex
c
      if (nlpoints .ne. npoi) then
         write(*,9000)
         write(*,*) 'number of land points in input/surta.nc'
         write(*,*) 'does not match number of land points in compar.h'
         write(*,*) 'in surta =',nlpoints,' in compar.h =',npoi
         stop 1
      else
         write (*,9010) 
         write (*,9020) nlpoints
         write (*,9010) 
      end if
cc
cc dummy variable example, 4-d, but 3rd dim ('level') = 1
cc copy and chanve for a new variable
cc
c      filen = 'input/dummyv.nc'
c      filen = dummyv_
c      aname = 'dummyv'
c      call readvar(filen,aname,'level',istart,icount,
c     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
c      if (istat.lt.0) then
c         write(*,9000)
c         print *, 'while reading dummyv'
c         stop 1
c      end if
c      call arr2vec (cdummy, xindummy)
c
c topography
c
      filen = 'input/topo.nc'
      filen = topo_
      aname = 'topo'
      call readvar(filen,aname,'level',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat.lt.0) then
         write(*,9000)
         print *, 'while reading topo'
         stop 1
      end if
      call arr2vec (cdummy, xintopo)
c
c fixed vegetation map
c
      if (isimveg .le. 1) then
         filen = 'input/vegtype.nc'
         filen = vegtype_
         aname = 'vegtype'
         call readvar(filen,aname,'level',istart,icount,cdummy,
     >    work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
         if (istat.lt.0) then
            write(*,9000)
            print *, 'while reading vegtype'
            stop 1
         end if
         call arr2vec (cdummy, xinveg)
      end if
cc
cc 2-d soil array
cc
c     filen = 'input/soil.nc'
c     filen = soil_
c     aname = 'soil'
c     call readvar(filen,aname,'level',istart,icount,
c    > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
c     if (istat.lt.0) then
c        write(*,9000)
c        print *, 'while reading soil'
c        stop 1
c     end if
c     call arr2vec (cdummy, soita)
c
c delta t
c
c     filen = 'input/deltat.nc'
      filen = deltat_mon_
      aname = 'deltat'
      call readvar(filen,aname,'level',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat.lt.0) then
         write(*,9000)
         print *, 'while reading deltat'
         stop 1
      end if
      call arr2vec (cdummy, deltat)
c
c 3-d soil texture array
c
c icount(3) is the 6 layers used in soita.sand.nc soita.clay.nc
c
      icount(3) = 6 
      icount(4) = 1
      filen = 'input/soita.sand.nc'
      filen = soita_sand_
      aname = 'sandpct'
      call readvar(filen,aname,'layer',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat.lt.0) then
         write(*,9000)
         print *, 'while reading soita.sand'
         stop 1
      end if
      do 12 j = 1, nsoilay
        call arr2vec (cdummy((j-1)*nlonsub*nlatsub + 1), sand(1,j))
 12   continue
c
      icount(3) = 6 
      icount(4) = 1
      filen = 'input/soita.clay.nc'
      filen = soita_clay_
      aname = 'claypct'
      call readvar(filen,aname,'layer',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat.lt.0) then
         write(*,9000)
         print *, 'while reading soita.clay'
         stop 1
      end if
      do 13 j = 1, nsoilay
        call arr2vec (cdummy((j-1)*nlonsub*nlatsub + 1), clay(1,j))
 13   continue
c
c 3-d climate arrays
c
      icount(3) = 1
      icount(4) = 12
c
      filen = 'input/wetd.mon.nc'
      filen = wetd_mon_
      aname = 'wetd'
      call readvar(filen,aname,'level',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat.lt.0) then
         write(*,9000)
         print *, 'while reading wetd'
         stop 1
      end if
      do 15 ntime = 1,12
         call arr2vec (cdummy((ntime-1)*nlonsub*nlatsub + 1),
     >    clmwet(1,ntime))
 15   continue
c
      filen = 'input/temp.mon.nc'
      filen = temp_mon_
      aname = 'temp'
      call readvar(filen,aname,'level',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat.lt.0) then
         write(*,9000)
         print *, 'while reading temp'
         stop 1
      end if
      do 20 ntime = 1,12
         call arr2vec (cdummy((ntime-1)*nlonsub*nlatsub + 1),
     >    clmt(1,ntime))
 20   continue
c
      filen = 'input/trange.mon.nc'
      filen = trange_mon_
      aname = 'trange'
      call readvar(filen,aname,'level',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat.lt.0) then
         write(*,9000)
         print *, 'while reading trange'
         stop 1
      end if
      do 25 ntime = 1,12
         call arr2vec (cdummy((ntime-1)*nlonsub*nlatsub + 1),
     >    clmtrng(1,ntime))
 25   continue
c
      filen = 'input/prec.mon.nc'
      filen = prec_mon_
      aname = 'prec'
      call readvar(filen,aname,'level',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat.lt.0) then
         write(*,9000)
         print *, 'while reading prec'
         stop 1
      end if
      do 30 ntime = 1,12
         call arr2vec (cdummy((ntime-1)*nlonsub*nlatsub + 1),
     >    clmprec(1,ntime))
 30   continue
c
      filen = 'input/wspd.mon.nc'
      filen = wspd_mon_
      aname = 'wspd'
      call readvar(filen,aname,'level',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat.lt.0) then
         write(*,9000)
         print *, 'while reading wspd'
         stop 1
      end if
      do 35 ntime = 1,12
         call arr2vec (cdummy((ntime-1)*nlonsub*nlatsub + 1),
     >    xinwind(1,ntime))
 35   continue
c
      filen = 'input/cld.mon.nc'
      filen = cld_mon_
      aname = 'cld'
      call readvar(filen,aname,'level',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat.lt.0) then
         write(*,9000)
         print *, 'while reading cld'
         stop 1
      end if
      do 40 ntime = 1,12
         call arr2vec (cdummy((ntime-1)*nlonsub*nlatsub + 1),
     >    clmcld(1,ntime))
 40   continue
c
      filen = 'input/rh.mon.nc'
      filen = rh_mon_
      aname = 'rh'
      call readvar(filen,aname,'level',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat.lt.0) then
         write(*,9000)
         print *, 'while reading rh'
         stop 1
      end if
      do 45 ntime = 1,12
         call arr2vec (cdummy((ntime-1)*nlonsub*nlatsub + 1),
     >    clmq(1,ntime))
 45   continue
c
c copy all 5 climatology fields to clm+anom fields for spin up
c
      ndim = npoi*12
      call scopy (ndim, clmt, xint)
      call scopy (ndim, clmtrng, xintrng)
      call scopy (ndim, clmprec, xinprec)
      call scopy (ndim, clmcld, xincld)
      call scopy (ndim, clmq, xinq)
      call scopy (ndim, clmwet, xinwet)
c
 9000 format (1x,'ERROR in subroutine readit')
 9010 format (1x,' ')
 9020 format (1x,'number of land points: ', i10)
c
c return to main program
c
      return
      end
c 
c
c ---------------------------------------------------------------------
      subroutine restart (iyrlast)
c ---------------------------------------------------------------------
c
c reads in restart files, initializes some variables
c
c this subroutine reads the restart values of:
c
c  fsnocov = fractional snow cover
c  tsno    = temperature of snow
c  hsno    = snow depth
c  tsoi    = soil temperature
c  wisoi   = soil ice content
c  wsoi    = soil moisture content
c  cbiol   = carbon in leaf biomass pool
c  cbiow   = carbon in woody biomass pool
c  cbior   = carbon in fine root biomass pool
c  sapfrac = sapwood fraction
c  clitlm  = leaf metabolic litter
c  clitls  = leaf structural litter
c  clitll  = leaf lignin litter
c  clitrm  = root metabolic litter
c  clitrs  = root structural litter
c  clitrl  = root lignin litter
c  clitwm  = woody metabolic litter
c  clitws  = woody structural litter
c  clitwl  = woody lignin litter
c  falll   = annual leaf litterfall 
c  fallr   = annual fine root turnover
c  fallw   = annual woody turnover
c  totcmic = total microbial carbon
c  csoislop= slow soil carbon, protected humus
c  csoislon= slow soil carbon, nonprotected humus
c  csoipas = passive soil carbon
c  gdd0    = growing degree days 0
c  gdd5    = growing degree days 5
c  tc      = coldest monthly temperature
c  tw      = warmest monthly temperature
c  wipud   = ice content of puddles per soil area
c  wpud    = liquid content of puddles per soil area
c  agddu   = annual accumulated growing degree days for bud burst, upper canopy
c  agddl   = annual accumulated growing degree days for bud burst, lower canopy
c  tempu   = cold-phenology trigger for trees
c  templ   = cold-phenology trigger for grasses/shrubs
c  a10td    = 10-day avg daily temp
c  a10ancub = 10-day average canopy photosynthesis rate - broadleaf
c  a10ancuc = 10-day average canopy photosynthesis rate - conifer
c  a10ancls = 10-day average canopy photosynthesis rate - shrubs
c  a10ancl4 = 10-day average canopy photosynthesis rate - c4 grasses
c  a10ancl3 = 10-day average canopy photosynthesis rate - c3 grasses
c  a10scalparamu = 10-day average canopy scaling parameter - upper canopy
c  a10scalparaml = 10-day average canopy scaling parameter - lower canopy
c  a10daylightu = 10-day average daylight - upper canopy
c  a10daylightl = 10-day average daylight - lower canopy
c  dropu   = drought-phenology trigger for trees
c  dropls  = drought-phenology trigger for shrubs
c  dropl4  = drought-phenology trigger for c4 grasses
c  dropl3  = drought-phenology trigger for c3 grasses
c (NOTE: a10ancuc is not used at this point, so its restart entry 
c is commented out)
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'comsoi.h'
      include 'comsno.h'
      include 'combcs.h'
      include 'comsum.h'
      include 'comveg.h'
      include 'comwork.h'
      include 'argvs.h'
c
c Arguments
c
      integer iyrlast
c
c Local variables
c
      integer istart(4), icount(4) ! for reading restart vars
c
      character*20 fdir  ! used to construct odd/even file names
      character*80 filen ! file name
c
      integer  lf,           ! number of characters in directory name
     >         i,            ! loop indices
     >         istat,        ! error flag for netcdf
     >         nlevel        ! loop indice on the soil layer
c
c External
c
      integer lenchr,       ! Function: Find length of character string
     >  NF_PUT_ATT_TEXT,    ! netcdf function
     >  NF_GLOBAL           ! netcdf function
c
c ---------------------------------------------------------------------
c
      data istart / 1,1,1,1 /, icount / nlon,nlat,1,1 /
      icount(1) = nlonsub
      icount(2) = nlatsub
c
c check to see if iyrlast is odd or even and read from the appropriate
c files for that year
c
      if (mod(iyrlast,2) .eq. 0) then
         fdir = 'restart/even'
      else
         fdir = 'restart/odd'
      end if
      lf = lenchr(fdir)
cc
cc dummy variable example, 3-d - copy and change for ne variable
cc
c      icount(3) = 1
c      filen = fdir(1:lf)//'dummyv.nc'
c      call readvar(filen,'dummyv','',istart,icount,
c     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in restart, dummyvv'
c         stop 1
c      end if
c      call arr2vec (cdummy, dummyv)
cc
cc dummy variable example, 4-d, 3rd dim (snowlayer) = nsnolay
cc
c      icount(3) = nsnolay
c      filen = fdir(1:lf)//'dummyv.nc'
c      call readvar(filen,'dummyv','snowlayer',istart,icount,
c     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in restart, dummyv'
c         stop 1
c      end if
c      do 5 nlevel = 1,nsnolay
c         call arr2vec (cdummy((nlevel-1)*nlonsub*nlatsub + 1),
c     >    dummyv(1,nlevel))
c 5    continue
c
c fsnocov
c
      filen = fdir(1:lf)//'fsnocov.nc'
      icount(3) = 1
      call readvar(filen,'fsnocov','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, fsnocov'
         stop 1
      end if
      call arr2vec (cdummy, fi)
c
c nsnolay variables: tsno and hsno
c
      icount(3) = nsnolay
c
      filen = fdir(1:lf)//'tsno.nc'
      call readvar(filen,'tsno','snowlayer',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, tsno'
         stop 1
      end if
      do 5 nlevel = 1,nsnolay
         call arr2vec (cdummy((nlevel-1)*nlonsub*nlatsub + 1),
     >    tsno(1,nlevel))
 5    continue
c
      filen = fdir(1:lf)//'hsno.nc'
      call readvar(filen,'hsno','snowlayer',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, hsno'
         stop 1
      end if
      do 10 nlevel = 1,nsnolay
         call arr2vec (cdummy((nlevel-1)*nlonsub*nlatsub + 1),
     >    hsno(1,nlevel))
 10   continue
c
c nsoilay variables: tsoi, wisoi, wsoi
c
      icount(3) = nsoilay
c
      filen = fdir(1:lf)//'tsoi.nc'
      call readvar(filen,'tsoi','soillayer',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, tsoi'
         stop 1
      end if
      do 15 nlevel = 1,nsoilay
         call arr2vec (cdummy((nlevel-1)*nlonsub*nlatsub + 1),
     >    tsoi(1,nlevel))
 15   continue
c
      filen = fdir(1:lf)//'wisoi.nc'
      call readvar(filen,'wisoi','soillayer',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, wisoi'
         stop 1
      end if
      do 20 nlevel = 1,nsoilay
         call arr2vec (cdummy((nlevel-1)*nlonsub*nlatsub + 1),
     >    wisoi(1,nlevel))
 20   continue
c
      filen = fdir(1:lf)//'wsoi.nc'
      call readvar(filen,'wsoi','soillayer',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, wsoi'
         stop 1
      end if
      do 25 nlevel = 1,nsoilay
         call arr2vec (cdummy((nlevel-1)*nlonsub*nlatsub + 1),
     >    wsoi(1,nlevel))
 25   continue
c
c npft variables
c
      icount(3) = npft
c
      filen = fdir(1:lf)//'cbiol.nc'
      call readvar(filen,'cbiol','pft',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, cbiol'
         stop 1
      end if
      do 30 nlevel = 1,npft
         call arr2vec (cdummy((nlevel-1)*nlonsub*nlatsub + 1),
     >    cbiol(1,nlevel))
 30   continue
c
      filen = fdir(1:lf)//'cbiow.nc'
      call readvar(filen,'cbiow','pft',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, cbiow'
         stop 1
      end if
      do 35 nlevel = 1,npft
         call arr2vec (cdummy((nlevel-1)*nlonsub*nlatsub + 1),
     >    cbiow(1,nlevel))
 35   continue
c
      filen = fdir(1:lf)//'cbior.nc'
      call readvar(filen,'cbior','pft',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, cbior'
         stop 1
      end if
      do 40 nlevel = 1,npft
         call arr2vec (cdummy((nlevel-1)*nlonsub*nlatsub + 1),
     >    cbior(1,nlevel))
 40   continue
c
c single level variables
c
      icount(3) = 1
c
      filen = fdir(1:lf)//'sapfrac.nc'
      call readvar(filen,'sapfrac','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, sapfrac'
         stop 1
      end if
      call arr2vec (cdummy, sapfrac)
c
      filen = fdir(1:lf)//'clitlm.nc'
      call readvar(filen,'clitlm','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, clitlm'
         stop 1
      end if
      call arr2vec (cdummy, clitlm)
c
      filen = fdir(1:lf)//'clitls.nc'
      call readvar(filen,'clitls','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, clitls'
         stop 1
      end if
      call arr2vec (cdummy, clitls)
c
      filen = fdir(1:lf)//'clitll.nc'
      call readvar(filen,'clitll','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, clitll'
         stop 1
      end if
      call arr2vec (cdummy, clitll)
c
      filen = fdir(1:lf)//'clitrm.nc'
      call readvar(filen,'clitrm','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, clitrm'
         stop 1
      end if
      call arr2vec (cdummy, clitrm)
c
      filen = fdir(1:lf)//'clitrs.nc'
      call readvar(filen,'clitrs','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, clitrs'
         stop 1
      end if
      call arr2vec (cdummy, clitrs)
c
      filen = fdir(1:lf)//'clitrl.nc'
      call readvar(filen,'clitrl','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, clitrl'
         stop 1
      end if
      call arr2vec (cdummy, clitrl)
c
      filen = fdir(1:lf)//'clitwm.nc'
      call readvar(filen,'clitwm','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, clitwm'
         stop 1
      end if
      call arr2vec (cdummy, clitwm)
c
      filen = fdir(1:lf)//'clitws.nc'
      call readvar(filen,'clitws','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, clitws'
         stop 1
      end if
      call arr2vec (cdummy, clitws)
c
      filen = fdir(1:lf)//'clitwl.nc'
      call readvar(filen,'clitwl','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, clitwl'
         stop 1
      end if
      call arr2vec (cdummy, clitwl)
c
      filen = fdir(1:lf)//'falll.nc'
      call readvar(filen,'falll','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, falll'
         stop 1
      end if
      call arr2vec (cdummy, falll)
c
      filen = fdir(1:lf)//'fallr.nc'
      call readvar(filen,'fallr','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, fallr'
         stop 1
      end if
      call arr2vec (cdummy, fallr)
c
      filen = fdir(1:lf)//'fallw.nc'
      call readvar(filen,'fallw','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, fallw'
         stop 1
      end if
      call arr2vec (cdummy, fallw)
c
      filen = fdir(1:lf)//'totcmic.nc'
      call readvar(filen,'totcmic','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, totcmic'
         stop 1
      end if
      call arr2vec (cdummy, totcmic)
c
      filen = fdir(1:lf)//'csoislop.nc'
      call readvar(filen,'csoislop','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, csoislop'
         stop 1
      end if
      call arr2vec (cdummy, csoislop)
c
      filen = fdir(1:lf)//'csoislon.nc'
      call readvar(filen,'csoislon','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, csoislon'
         stop 1
      end if
      call arr2vec (cdummy, csoislon)
c
      filen = fdir(1:lf)//'csoipas.nc'
      call readvar(filen,'csoipas','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, csoipas'
         stop 1
      end if
      call arr2vec (cdummy, csoipas)
c
      filen = fdir(1:lf)//'gdd0.nc'
      call readvar(filen,'gdd0','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, gdd0'
         stop 1
      end if
      call arr2vec (cdummy, gdd0)
c
      filen = fdir(1:lf)//'gdd5.nc'
      call readvar(filen,'gdd5','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, gdd5'
         stop 1
      end if
      call arr2vec (cdummy, gdd5)
c
      filen = fdir(1:lf)//'tc.nc'
      call readvar(filen,'tc','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, tc'
         stop 1
      end if
      call arr2vec (cdummy, tc)
c
      filen = fdir(1:lf)//'tw.nc'
      call readvar(filen,'tw','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, tw'
         stop 1
      end if
      call arr2vec (cdummy, tw)
c
      filen = fdir(1:lf)//'wipud.nc'
      call readvar(filen,'wipud','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, wipud'
         stop 1
      end if
      call arr2vec (cdummy, wipud)
c
      filen = fdir(1:lf)//'wpud.nc'
      call readvar(filen,'wpud','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, wpud'
         stop 1
      end if
      call arr2vec (cdummy, wpud)
c
      filen = fdir(1:lf)//'agddu.nc'
      call readvar(filen,'agddu','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, agddu'
         stop 1
      end if
      call arr2vec (cdummy, agddu)
c
      filen = fdir(1:lf)//'agddl.nc'
      call readvar(filen,'agddl','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, agddl'
         stop 1
      end if
      call arr2vec (cdummy, agddl)
c
      filen = fdir(1:lf)//'tempu.nc'
      call readvar(filen,'tempu','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, tempu'
         stop 1
      end if
      call arr2vec (cdummy, tempu)
c
      filen = fdir(1:lf)//'templ.nc'
      call readvar(filen,'templ','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, templ'
         stop 1
      end if
      call arr2vec (cdummy, templ)
c
      filen = fdir(1:lf)//'a10td.nc'
      call readvar(filen,'a10td','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, a10td'
         stop 1
      end if
      call arr2vec (cdummy, a10td)
c
      filen = fdir(1:lf)//'a10ancub.nc'
      call readvar(filen,'a10ancub','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, a10ancub'
         stop 1
      end if
      call arr2vec (cdummy, a10ancub)
ccc
cc      filen = fdir(1:lf)//'a10ancuc.nc'
cc      call readvar(filen,'a10ancuc','',istart,icount,
cc     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
cc      if (istat .ne. 0) then
cc         write(*,*) 'ERROR in restart, a10ancuc'
cc         stop 1
cc      end if
cc      call arr2vec (cdummy, a10ancuc)
c
      filen = fdir(1:lf)//'a10ancls.nc'
      call readvar(filen,'a10ancls','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, a10ancls'
         stop 1
      end if
      call arr2vec (cdummy, a10ancls)
c
      filen = fdir(1:lf)//'a10ancl4.nc'
      call readvar(filen,'a10ancl4','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, a10ancl4'
         stop 1
      end if
      call arr2vec (cdummy, a10ancl4)
c
      filen = fdir(1:lf)//'a10ancl3.nc'
      call readvar(filen,'a10ancl3','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, a10ancl3'
         stop 1
      end if
      call arr2vec (cdummy, a10ancl3)
c
      filen = fdir(1:lf)//'a10scalparamu.nc'
      call readvar(filen,'a10scalparamu','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, a10scalparamu'
         stop 1
      end if
      call arr2vec (cdummy, a10scalparamu)
c
      filen = fdir(1:lf)//'a10scalparaml.nc'
      call readvar(filen,'a10scalparaml','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, a10scalparaml'
         stop 1
      end if
      call arr2vec (cdummy, a10scalparaml)
c
      filen = fdir(1:lf)//'a10daylightu.nc'
      call readvar(filen,'a10daylightu','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, a10daylightu'
         stop 1
      end if
      call arr2vec (cdummy, a10daylightu)
c
      filen = fdir(1:lf)//'a10daylightl.nc'
      call readvar(filen,'a10daylightl','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, a10daylightl'
         stop 1
      end if
      call arr2vec (cdummy, a10daylightl)
c
      filen = fdir(1:lf)//'dropu.nc'
      call readvar(filen,'dropu','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, dropu'
         stop 1
      end if
      call arr2vec (cdummy, dropu)
c
      filen = fdir(1:lf)//'dropls.nc'
      call readvar(filen,'dropls','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, dropls'
         stop 1
      end if
      call arr2vec (cdummy, dropls)
c
      filen = fdir(1:lf)//'dropl4.nc'
      call readvar(filen,'dropl4','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, dropl4'
         stop 1
      end if
      call arr2vec (cdummy, dropl4)
c
      filen = fdir(1:lf)//'dropl3.nc'
      call readvar(filen,'dropl3','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, dropl3'
         stop 1
      end if
      call arr2vec (cdummy, dropl3)
c
c calculate tcmin
c
      do i = 1, npoi
         tcmin(i) = tc(i) + deltat(i)
      enddo
c
      call existence
c
      return
      end
c
c ---------------------------------------------------------------------
      subroutine coldstart
c ---------------------------------------------------------------------
c  
      include 'implicit.h'
c
      include 'compar.h'
      include 'comsoi.h'
      include 'comsno.h'
      include 'argvs.h'
c
c initialize some model variables for cold start conditions
c
      call const (fi, npoi, 0.0)
c
      call const (hsno, npoi*nsnolay, 0.0)
      call const (tsno, npoi*nsnolay, 273.16)
c
      call const (tsoi, npoi*nsoilay, 278.16)
      call const (wsoi, npoi*nsoilay, 0.50)
      call const (wisoi, npoi*nsoilay, 0.00)
c
c return to main program
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine rdanom(imonth,iyear,iyranom,nanom,iy2,istyr,iwest,jnorth)
c ---------------------------------------------------------------------
c
c reads in anomalies for imonth+1
c
c this subroutine reads the monthly anomaly values of:
c
c temp - mean temperature
c trng - temperature range (not read in yet, saved for future use)
c prec - precipitation rate
c cld  - cloudiness
c rh   - relative humidity
c wspd - wind speed (not read in yet, saved for future use)
c wetd - wet days per month (not read in yet, saved for future use)
c
c and adds them to the climatological value for the month
c
c  
      include 'implicit.h'
c
      include 'compar.h'
      include 'combcs.h'
      include 'comveg.h'
      include 'comwork.h'
      include 'argvs.h'
c
c Arguments
c
      integer imonth, ! month
     >        iyear,  ! year
     >        iyranom,! year to start reading anoms
     >        iy2,    ! final year of the run
     >        istyr,  ! 1st year in data files
     >        iwest,  ! 1st lon index for subset
     >        jnorth, ! 1st lat index for subset
     >        nanom   ! # of years in the anomaly files
c
c Local variables
c
      real anom(npoi)
      equivalence(anom(1),cdummy(1))
c
      integer imon,   ! month + 1
     >        iyr,    ! number of years after begin of data file(istyr)
     >        istat,  ! error flag for netcdf
     >        i,      ! loop indice on land points
     >        jyear   ! iyr divided by nanom (for looping through anomaly files)
c
      integer istart(4), icount(4) ! for reading rdanom vars
      character*80 filen
c
c ---------------------------------------------------------------------
c
      data istart / 1,1,1,1 /, icount / nlon,nlat,1,1 /
      istart(1) = iwest
      istart(2) = jnorth
      icount(1) = nlonsub
      icount(2) = nlatsub
c
c determine which (if any) month to read
c If prior to November of the year before an anomaly year, then return.
c Else read anomalies for the following month.
c
      if (iyear .lt. iyranom-1) then
         return
      else if (iyear .eq. iyranom-1 .and. imonth .lt. 11) then
         return
c
c If timestep equals December of final year in anomaly file and also equals
c the final year of the run, then set January anomalies to zero and return.
c If not the final year of the run, then loop back to start of anomaly file
c (see jyear).
c
      else if (iyear .eq. istyr+nanom-1 .and. imonth .eq. 12) then
        if (iyear .eq. iy2) then
c
          print *, 'WARNING: last month of run; no anomalies for January'
          print *, 'Using climatologies for month year ='
          print *, imonth+1,iyear+1
          do 4 i = 1, npoi
            xint(i,1) = clmt(i,1)
c           xintrng(i,1) = clmtrng(i,1)
            xinprec(i,1) = clmprec(i,1)
            xincld(i,1) = clmcld(i,1)
            xinq(i,1) = clmq(i,1)
c           xinwind(i,1) = clmw(i,1)
c           xinwet(i,1) = clmwet(i,1)
 4        continue
          return
        end if
      end if
c
      iyr = iyear-istyr
      imon = imonth + 1
      if (imon .eq. 13) then
         imon = 1
         iyr = iyr + 1
      end if
      jyear = iyr/nanom
      istart(4) = (iyr - nanom*jyear)*12 + imon
c
      if (iyr.gt.0 .and. (iyr - nanom*jyear).eq.0) then
        print *, 'WARNING: Attempted to read past last month in anomaly file'
        print *, 'Looping back to the beginning of the file'
      end if
c
      if (istart(4) .gt. 0) then
        print *, 'rdanom reading month year step ='
        print *, imon,iyr+istyr-nanom*jyear,istart(4)
      else
        print *, 'WARNING, anomalies begin in year ',istyr
        print *, 'Not reading in anomalies for month year ='
        print *, imon,iyr+istyr
        return
      end if
cc
cc     dummy variable example, 4-d, whose 3rd dim (level) = 1
cc
c      aname = 'dummyv'
c      filen = 'input/anom/dummyv.danom.nc'
c      call readvar(filen,aname,'level',istart,icount,
c     > work,cdummy(1),cdummy(nlonsub+1),cdummy(2*nlonsub+1),
c     > cdummy(3*nlonsub+1),istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in rdanom, dummyv'
c         stop 1
c      end if
c      call arr2vec( work, anom )
c      do 10 i = 1, npoi
c         xindummyv(i,imon) = max (clmtrng(i,imon) + anom(i), 0.1)
c 10   continue
c
c     mean temperature
c
      aname = 'temp'
c     filen = 'input/anom/temp.danom.nc'
      filen = temp_danom_
      call readvar(filen,aname,'level',istart,icount,
     > work,cdummy(1),cdummy(nlonsub+1),cdummy(2*nlonsub+1),
     > cdummy(3*nlonsub+1),istat)
      if (istat .ne. 0) then
        write(*,*) 'ERROR in rdanom, temp'
        stop 1
      end if
      call arr2vec( work, anom )
      do 5 i = 1, npoi
         xint(i,imon) = clmt(i,imon) + anom(i)
 5    continue
cc
cc     temperature range
cc
c      aname = 'trange'
c      filen = 'input/anom/trange.danom.nc'
       filen = trange_danom_
c      call readvar(filen,aname,'level',istart,icount,
c     > work,cdummy(1),cdummy(nlonsub+1),cdummy(2*nlonsub+1),
c     > cdummy(3*nlonsub+1),istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in rdanom, trange'
c         stop 1
c      end if
c      call arr2vec( work, anom )
c      do 10 i = 1, npoi
c         xintrng(i,imon) = max (clmtrng(i,imon) + anom(i), 0.1)
c 10   continue
c
c     precipitation rate
c
      aname = 'prec'
c     filen = 'input/anom/prec.danom.nc'
      filen = prec_danom_
      call readvar(filen,aname,'level',istart,icount,
     > work,cdummy(1),cdummy(nlonsub+1),cdummy(2*nlonsub+1),
     > cdummy(3*nlonsub+1),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in rdanom, prec'
         stop 1
      end if
      call arr2vec( work, anom )
      do 15 i = 1, npoi
         xinprec(i,imon) = clmprec(i,imon) + anom(i)
 15   continue
c
c     cloudiness
c
      aname = 'cld'
c     filen = 'input/anom/cld.danom.nc'
      filen = cld_danom_
      call readvar(filen,aname,'level',istart,icount,
     > work,cdummy(1),cdummy(nlonsub+1),cdummy(2*nlonsub+1),
     > cdummy(3*nlonsub+1),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in rdanom, cld'
         stop 1
      end if
      call arr2vec( work, anom )
      do 20 i = 1, npoi
         xincld(i,imon) = clmcld(i,imon) + anom(i)
 20   continue
c
c     relative humidity
c
      aname = 'rh'
c     filen = 'input/anom/rh.danom.nc'
      filen = rh_danom_
      call readvar(filen,aname,'level',istart,icount,
     > work,cdummy(1),cdummy(nlonsub+1),cdummy(2*nlonsub+1),
     > cdummy(3*nlonsub+1),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in rdanom, rh'
         stop 1
      end if
      call arr2vec( work, anom )
      do 25 i = 1, npoi
         xinq(i,imon) = clmq(i,imon) + anom(i)
 25   continue
c
c     wind speed
c
c     aname = 'wspd'
c     filen = 'input/anom/wspd.danom.nc'
      filen = wspd_danom_
c     call readvar(filen,aname,'level',istart,icount,
c    > work,cdummy(1),cdummy(nlonsub+1),cdummy(2*nlonsub+1),
c    > cdummy(3*nlonsub+1),istat)
c     if (istat .ne. 0) then
c        write(*,*) 'ERROR in rdanom, wspd'
c        stop 1
c     end if
c     call arr2vec( work, anom )
c     do 30 i = 1, npoi
c        xinwind(i,imon) = clmw(i,imon) + anom(i)
c30   continue
c
c     wet days
c
c     aname = 'wetd'
c     filen = 'input/anom/wetd.danom.nc'
      filen = wetd_danom_
c     call readvar(filen,aname,'level',istart,icount,
c    > work,cdummy(1),cdummy(nlonsub+1),cdummy(2*nlonsub+1),
c    > cdummy(3*nlonsub+1),istat)
c     if (istat .ne. 0) then
c        write(*,*) 'ERROR in rdanom, wetd'
c        stop 1
c     end if
c     call arr2vec( work, anom )
c     do 35 i = 1, npoi
c        xinwet(i,imon) = clmwet(i,imon) + anom(i)
c35   continue
c
c
      return
      end
c
c
c
c ---------------------------------------------------------------------
      subroutine inird(file,istyr)
c ---------------------------------------------------------------------
c
c This subroutine obtains the year+1 of the year in the units attribute
c of a file.  Low-level netcdf commands are used.
c
      include 'implicit.h'
      include 'argvs.h'
c
      include 'netcdf.inc'
c
c Arguments
c
      character*(*) file
      integer istyr       
c
c Local Variables
c
      integer idies, istat, idtime, lf1
c
* Remove this block for Absoft f77, D.Pol 15-Mar-2002
* unsure if needed for SGI version
*      integer NF_OPEN,        ! netcdf function
*     >        NF_NOWRITE,     ! '
*     >        NF_NOERR,       ! '
*     >        NF_INQ_VARID,   ! '
*     >        NF_GET_ATT_TEXT ! '
c              
      character*80 units
c ---------------------------------------------------------------------
c
c open file
c
      istat = NF_OPEN(file,NF_NOWRITE,idies)
      if (istat .ne. NF_NOERR) then
         print *, 'Error in inird while trying to open file'
         print *, file
         print *, NF_STRERROR(istat)
         istyr = -1
         return
      end if
c
c get units attribute for time
c
      istat = NF_INQ_VARID(idies,'time',idtime)
      if (istat .ne. NF_NOERR) then
         print *, 'Error in inird while trying to get time id'
         print *, NF_STRERROR(istat)
         istyr = -1
         return
      end if
      units = ' '
      istat = NF_GET_ATT_TEXT(idies,idtime,'units',units)
      if (istat .ne. NF_NOERR) then
         print *, 'Error in inird while trying to get time units'
         print *, NF_STRERROR(istat)
         istyr = -1
         return
      end if
c
c put character units year into integer variable, add 1
c
      lf1 = index(units,'since') + 6
      read(units(lf1:lf1+3),'(i4)') istyr
c      read(units(12:15),'(i4)') istyr
      istyr = istyr + 1
      return
      end
c
c ---------------------------------------------------------------------
      subroutine rdday(jday, imonth, iyear, istyr, iwest, jnorth)
c ---------------------------------------------------------------------
c
c This subroutine reads in daily fields..
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'combcs.h'
      include 'comwork.h'
      include 'argvs.h'
c
c Arguments
c
      integer jday,   ! day of the year
     >        imonth, ! month
     >        iyear,  ! year
     >        istyr,  ! 1st year in data files
     >        iwest,  ! 1st lon index for subset
     >        jnorth  ! 1st lat index for subset
c
c Local variables
c
      integer i,      ! loop indice on years after istyr
     >        istat   ! error flag for netcdf
c
      character*80 filen
      integer istart(4), icount(4)
c
c ---------------------------------------------------------------------
c
      data istart / 1,1,1,1 /, icount / nlon,nlat,1,1 /
      istart(1) = iwest
      istart(2) = jnorth
      icount(1) = nlonsub
      icount(2) = nlatsub
c
      if (iyear .lt. istyr) then
        print *, 'daily data begins in year ', istyr
        print *, 'not reading in daily data'
        return
      end if
c
c count how many days since beginning of daily data
c daily data begin on Jan 1, istyr
c
      if (iyear .eq. istyr) then
        istart(4) = jday
      else
        istart(4) = 0
        do 10 i = istyr, iyear-1
          istart(4) = istart(4) + 365
          if (mod(i,4).eq.0) then
            if (mod(i,100).ne.0) then
              istart(4) = istart(4) + 1
            else if (mod(i/100,4).eq.0) then
              istart(4) = istart(4) + 1
            end if
          end if
 10     continue
        istart(4) = istart(4) + jday
      end if
c
      if (istart(4) .gt. 0) then
        print *, 'rdday reading day month year step ='
        print *, jday, imonth, iyear+istyr,istart(4)
      else
        print *, 'WARNING, anomalies begin in year ',istyr
        print *, 'Not reading in anomalies for day month year ='
        print *, jday, imonth, iyear+istyr
        return
      end if
c
c read daily precip
c
      aname = 'prec'
c     filen = 'input/anom/daily/prec.fanom.nc'
      filen = prec_fanom_
      call readvar(filen,aname,'level',istart,icount,
     > work,cdummy(1),cdummy(nlonsub+1),cdummy(2*nlonsub+1),
     > cdummy(3*nlonsub+1),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in rdday, prec'
         if (iyear .gt. 1997) then
           print *, 'Attempted to read past last day in file?'
         end if
         return
      end if
      call arr2vec( work, xinprecd(1) )
c
c read daily temp
c
      aname = 'temp'
c     filen = 'input/anom/daily/temp.danom.nc'
      filen = temp_danom_
      call readvar(filen,aname,'level',istart,icount,
     > work,cdummy(1),cdummy(nlonsub+1),cdummy(2*nlonsub+1),
     > cdummy(3*nlonsub+1),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in rdday, temp'
         if (iyear .gt. 1997) then
           print *, 'Attempted to read past last day in file?'
         end if
         return
      end if
      call arr2vec( work, xintd(1) )
c
c read daily trange
c
      aname = 'trange'
c     filen = 'input/anom/daily/trange.fanomc.nc'
      filen = trange_fanomc_
      call readvar(filen,aname,'level',istart,icount,
     > work,cdummy(1),cdummy(nlonsub+1),cdummy(2*nlonsub+1),
     > cdummy(3*nlonsub+1),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in rdday, trange'
         if (iyear .gt. 1997) then
           print *, 'Attempted to read past last day in file?'
         end if
         return
      end if
      call arr2vec( work, xintrngd(1) )
c
c read daily cloudiness
c
      aname = 'cld'
c     filen = 'input/anom/daily/cld.danom.nc'
      filen = cld_danom_
      call readvar(filen,aname,'level',istart,icount,
     > work,cdummy(1),cdummy(nlonsub+1),cdummy(2*nlonsub+1),
     > cdummy(3*nlonsub+1),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in rdday, cld'
         if (iyear .gt. 1997) then
           print *, 'Attempted to read past last day in file?'
         end if
         return
      end if
      call arr2vec( work, xincldd(1) )
c
c read daily windspeed
c
      aname = 'wspd'
c     filen = 'input/anom/daily/wspd.fanomc.nc'
      filen = wspd_fanomc_
      call readvar(filen,aname,'level',istart,icount,
     > work,cdummy(1),cdummy(nlonsub+1),cdummy(2*nlonsub+1),
     > cdummy(3*nlonsub+1),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in rdday, wspd'
         if (iyear .gt. 1997) then
           print *, 'Attempted to read past last day in file?'
         end if
         return
      end if
      call arr2vec( work, xinwindd(1) )
c
c read daily humidity
c
      aname = 'sphum'
c     filen = 'input/anom/daily/sphum.fanom.nc'
      filen = sphum_fanom_
      call readvar(filen,aname,'level',istart,icount,
     > work,cdummy(1),cdummy(nlonsub+1),cdummy(2*nlonsub+1),
     > cdummy(3*nlonsub+1),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in rdday, sphum'
         if (iyear .gt. 1997) then
           print *, 'Attempted to read past last day in file?'
         end if
         return
      end if
      call arr2vec( work, xinqd(1) )
c
      return
      end
      
      
c ---------------------------------------------------------------------
      subroutine simdate(cdate)
c ---------------------------------------------------------------------
c D.P.  3/15/02 : This subroutine replaces the intrinsic date() with fdate 
c this Y2K-compliant version. It calls another built-in c subroutine fdate() 
c that is Y2K-compliant. We also change cdate to 11 characters from 10.
c the format is now:    15-Mar-2002
c
      include 'implicit.h'
      include 'argvs.h'
      character*11 cdate
      character*24 longdate
c
      call FDATE(longdate)
      cdate=longdate(9:10)//'-'//longdate(5:7)//'-'//longdate(21:24)
c
      return
      end
c
c
c
c -------------------------------------------------------------------------
c HOW TO READ/WRITE NETCDF FILES IN IBIS
c -------------------------------------------------------------------------
c Reading/writing files in ibis is done through subroutines in
c ies-io.f.  The most important concept to understand is the
c relationship between the locations of points in an n-dimensional
c array and the values of istart and icount.  In FORTRAN, values in an
c array are stored so that the first dimension varies the fastest.  In
c C, the last dimension varies the fastest.  If you were to use ncdump
c on a netcdf file whose variable 'mydata' has 4 dimensions (latitude,
c longitude, level, time), the variable would be shown as follows:
c mydata(time, level, latitude, longitude)
c because ncdump was written in C and reflects C's ordering of dimensions.
c In this example time varies the the slowest and longitude the fastest.
c If you were to define this variable in a FORTRAN program so that
c time again varies the slowest and longitude the fastest, it would
c look like this:
c     real mydata(nlons, nlats, nlevels, ntimes)
c where nlons, nlats, nlevels, and ntimes are integer parameters of
c some specified size.  Looping through the data in the order of
c storage in memory would look like this:
c     do l = 1, ntimes
c       do k = 1, nlevels
c         do j = 1, nlats
c           do i = 1, nlons
c             mydata(i,j,k,l) = ...
c           enddo
c         enddo
c       enddo
c     enddo
c Since ies-io.f is FORTRAN code, keep in mind the FORTRAN
c representation will be flipped in comparison to what you see from
c ncdump.
c The netcdf interface reads and writes data using two integer vectors,
c istart, and icount.  Istart refers to the starting point for reading/writing
c along each dimension.  If you have a 4-d variable as above and want
c to write starting at the first point, istart would be a vector of
c length 4 and have the values (1,1,1,1). Icount refers to how far
c along each dimension data will be read/written.  If you wish to
c write the entire array in the sample FORTRAN code above, icount
c would be a vector of length 4 and have the values
c (nlons,nlats,nlevels,ntimes).
c Things get a little more complicated when you want to read/write
c only a portion of the variable.  In ibis we typically read/write a
c single time step, but possibly more than one level/pft and either an
c entire or partial lat/lon grid. Below are examples of istart and
c icount for various situations of reading/writing.
c 1) an entire lat/lon grid for only one (6th) time step and one (2nd) pft:
c istart is (1,1,2,6), icount is (nlons,nlats,1,1)
c
c 2) entire lat/lon arrays for all 9 pfts at one (6th) time step:
c istart is (1,1,1,6), icount is (nlons,nlats,npfts,1)
c
c 3) a single lat/lon point (at the ilon-th longitude and the ilat-th
c latitude) of a 3-d variable at all 100 time steps:
c istart is (ilon,ilat,1), icount is (1,1,100)
c Note that if istart and icount have been declared length 4, the 4th
c value in each is ignored when referencing a 3-d variable.
c
c 4) a subsection of a lat/lon grid, 20 points wide in longitude, 15
c points high in latitude, starting at point (ilon,ilat), at one (18th)
c level and 12 times, beginning at time step itime:
c istart is (ilon,ilat,18,itime) icount is (20,15,1,12)
c
c HOW TO ADD NEW CODE TO READ A FILE:
c To read a file, use subroutine readvar in ies-io.f. This subroutine
c assumes that the variable being read has dimensions called longitude,
c latitude, possibly time, and possibly another dimension, which is
c named in the call.  Only the bare essentials are returned.
c
c General call:
c     call readvar(filen,varname,name3d,istart,icount,values,
c    > alons,alats,vals3d,times,ierror)
c
c INPUT
c     filen - character*(*) - file name from which to read data
c     varname - character*(*) - name of variable from which to read
c     name3d - character*(*) - name of 3rd, nontime dimension variable.
c      Ignored if varname variable is only 3-d.
c     istart - integer(*) - starting points along each dimension of
c      the variable, for example the 1st point a 4-d variable is (1,1,1,1) 
c      and the 3rd level and 2nd time step of a 4d variable is (1,1,3,2).
c     icount - integer(*) - number of points along each dimension to read,
c      for example, to read in a single lat/lon grid from a 4-d variable,
c      icount would be (nlon,nlat,1,1)
c OUTPUT
c     values - real(*) - returned real values of the designated hyperslab
c     alons - real(*) - longitude vector for those points read in
c     alats - real(*) - latitude vector for those points read in
c     vals3d - real(*) - vector of values of the 3rd dimension (unchanged if
c      variable is 2- or 3-d, or if 3rd dimension has character values).
c     times - real(*) - time vector vector for those points read in (unchanged
c      if variable is 2-d).
c     ierror - integer - error code = 0 if no error, < 0 if error occurred
c
c Example: Read in an entire lat/lon grid for one (9th) pft (dimension 
c name is 'pft') and one (24th) time
c     parameter (nlons = 360, nlats=180)
c     real x(nlons,nlats), alons(nlons), alats(nlats), xjunk, time
c     integer istart(4), icount(4), ierr
c     istart(1) = 1
c     istart(2) = 1
c     istart(3) = 9
c     istart(4) = 24
c     icount(1) = nlons
c     icount(2) = nlats
c     icount(3) = 1
c     icount(4) = 1
c     call readvar('myfile.nc','myvar','pft',istart,icount,x,alons,alats,
c    > xjunk,time,ierr)
c     if (ierr .lt. 0) then
c       print *, 'Error occurred in readvar'
c     end if
c Note that above, I used a junk variable for the values of the 3rd
c dimension, because in ibis, the pft dimension is a character
c dimension.  If the 3rd dimension type is character, readvar does not
c return the value.
c
c Ibis-specific example: read in a variable which has 3rd dimension
c 'layer', using cdummy to store full array before extracting land
c points into the variable used directly in ibis.  The variable work is
c used to store returned info that won't be used later.
c Previously declared
c      istart(1) = 1        !begin at 1st longitude point
c      istart(2) = 1        !begin at 1st latitude point
c      istart(3) = 1        !begin at 1st layer
c      istart(4) = 1        !begin at 1st time
c      icount(1) = nlon     !read nlon points along longitude dimension
c      icount(2) = nlat     !read nlat points along latitude dimension
c Read in data in file input/soita.nc to ibis variable tex
c      icount(3) = nsoilay  !read nsoilay points along layer dimension
c      icount(4) = 1        !read one time step
c      filen = 'input/soita.nc'
c      aname = 'tex'
c      call readvar(filen,aname,'layer',istart,icount,
c     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
c      if (istat.lt.0) then
c         write(*,9000)
c         print *, 'while reading soita'
c         stop 1
c      end if
c For each lat/lon grid, strip off ocean points. tex only has land points.
c      do 12 j = 1, nsoilay
c        call arr2vec (cdummy((j-1)*nlonsub*nlatsub + 1), tex(1,j))
c 12   continue
c
c HOW TO ADD NEW CODE TO WRITE TO A FILE:
c Writing to a file involves several steps.  First, if the file does
c not yet exist, you must create the file using inifile or inifilec.
c Then, you must create one or more variables using inivar.  Then after
c all initializing is done, declare the end of the initialization phase
c using endini.
c Finally, once the file and variable exist, you can write data into
c the variable using writevar.
c Creating and writing to files in seperate steps may seem
c complicated at first, but this 4-step method allows much flexibility.
c Not only do you have the choice of creating a variable without a 3rd
c dimension, with a real 3rd dimension, or with a character 3rd
c dimension, you may also have more than one variable within the same
c file.
c Note that when using ies-io.f routines if you have more
c than one variable in the file, and both have 4 dimensions, that 3rd
c dimension (level, pft, etc.) must be shared by both variables.  For
c example, you may have one variable whose 3rd dimension is pft and
c another variable which does not have a 3rd non-time dimension in the
c same file, but you may not have in the same file one variable whose
c 3rd dimension is pft and another variable whose 3rd dimension is level.
c
c 1) Initialize a file
c If you wish a file to contain only latitude, longitude, and time
c dimensions, or want it to also have a 3rd real dimension, use
c inifile.  If you want to the file to have a 3rd character dimension,
c use inifilec.
c General call for inifile:
c      call inifile(idies,filen,title,source,history,nlon,alons,
c     > nlat,alats,name3rd,long3rd,units3rd,n3rd,vals3rd,pos3rd,tunits,
c     > calendar,ierror)
c
c INPUT
c     filen - character*(*) - name for new file
c     title - character*(*) - title for file
c     source - character*(*) - source for file
c     history - character*(*) - date of creation, and possibly author
c     nlon - integer - number of point along longitude direction
c     alons - real(nlon) - values of longitudes
c     nlat - integer - number of point along latitude direction
c     alats - real(nlat) - values of latitudes
c     name3rd - character*(*) - name of 3rd dimension - use '' if you
c      don't want this dimension ('' is the empty string)
c     long3rd - character*(*) - long name for 3rd dimension variable
c      (ignored if nam3rd='')
c     n3rd - integer - size of 3rd dimension (ignored if name3rd='')
c     vals3rd - real(n3rd) - values along 3rd dimension (ignored if 
c      name3rd='')
c     pos3rd - character*(*) - if the 3rd dimension exists and refers to
c      a measured level or height, use 'up' if values increase with distance
c      above earth (as in height), or 'down' if values decrease with distance 
c      (as in pressure) (ignored if name3rd=''). If 3rd dimension does not
c      refer to a height or level, use ''.
c     tunits - character*(*) - units for time, must be in the form of days
c      since yyyy-mm-dd tt:tt:tt, where yyyy is a year or 0000, mm is a
c      month, dd is a day, and tt:tt:tt is (optional) time or 00:00:00.
c     calendar - character*(*) - type of calendar.  Choose from 'noleap',
c      'gregorian','n ka BP', etc.  Use iescal (in ies.f) if orbital 
c      parameters differ from modern.
c OUTPUT
c     idies - integer - id number of new file for later use
c     ierror - integer - error code, 0 = no error, < 0 = an error occured
c
c Example: Create a file which will hold fractional snow
c cover (only 3-d) and snow thickness for each layer (4-d).
c previously defined: 
c
c     parameter (nlons = 360, nlats = 180, nsnolay = 6)
c     real lonscale(nlons), latscale(nlats), slayers(nsnolay)
c     real snowc(nlons,nlats), snowh(nlons,nlats,nsnolay)
c     integer istat
c     character*80 cdate, tunits
c     alats = ...   ! create values for latitude somehow
c     alons = ...   ! create values for longitude somehow
c     slayers = ... ! create values for snow layers somehow
c     cdate = 'created on 6/29/97'
c     tunits = 'days since 1969-12-31'
c Now initialize a file with a 3rd real variable
c     filen = 'snow.nc'
c     call inifile(idies,filen,
c    > 'file for snow variables',
c    > 'C Molling, program.f v1.01',cdate,nlons,lonscale,nlats,latscale,
c    > 'snowlayer','snow layers top to bottom','',nsnolay,slayers,
c    > 'down',tunits,'gregorian',istat)
c     if (istat .lt. 0)then
c       print *, 'Error in inifile'
c     end if
c Note above, the empty string ('') is used because the 3rd dimension,
c snowlayer, does not have any units.  The 'positive' attribute for
c the 3rd dimension is 'down' because snowlayer 1 is at the top and the last
c layer is at the bottom (distance above the center of the earth is
c decreasing as snowlayer increases).  The calendar is gregorian
c because there are leap years included.  If all years are 365 days, you
c should use 'noleap'.  Units of time are days dince a date of the
c form yyyy-mm-dd.  Time units should use this form to be compatible with
c GrADS.  Other units should be compatible with Udunits (see
c www.unidata.ucar.edu).  The returned variable idies will be used in
c the subroutine that initializes a variable.
c
c General call for inifilec
c     call inifilec(idies,filen,title,source,history,nlon,alons,
c    > nlat,alats,name3rd,long3rd,units3rd,n3rd,len3rd,chars3rd,
c    > tunits,calendar,ierror)
c
c INPUT
c     filen - character*(*) - name for new file
c     title - character*(*) - title for file
c     source - character*(*) - source for file
c     history - character*(*) - date of creation, and possibly author
c     nlon - integer - number of point along longitude direction
c     alons - real(nlon) - values of longitudes
c     nlat - integer - number of point along latitude direction
c     alats - real(nlat) - values of latitudes
c     name3rd - character*(*) - name of 3rd dimension
c     long3rd - charqcter*(*) - long name for 3rd dimension variable
c     n3rd - integer - size of 3rd dimension
c     len3rd - integer length of chracter strings in vals3rd
c     chars3rd - character*len3rd(n3rd) - values along 3rd dimension
c     tunits - character*(*) - units for time, must be in the form of days
c      since yyyy-mm-dd tt:tt:tt, where yyyy is a year or 0000, mm is a
c      month, dd is a day, and tt:tt:tt is (optional) time or 00:00:00
c     calendar - charater*(*) - type of calendar.  Choose from 'noleap',
c      'gregorian','n ka BP', etc.  Use iescal if orbital parameters
c      differ from modern.
c OUTPUT
c     idies - integer - id number of new file for later use
c     ierror - integer - error code, 0 = no error, < 0 = an error occured
c
c Example: Create a file that has a character 3rd dimension.
c Defined previously - most variables as in previous example, plus...
c      parameter (npft = 9)
c      character*80 pftdef(npft)
c      pftdef(1) = 'boreal evergreens'
c      pftdef(2) = ...
c      etc...
c Initialize file with inifilec
c      filen = 'exist.nc'
c      call inifilec(idies,filen,
c     > 'annual existence of each plant functional type',
c     > 'ibis wyearly',cdate,nlon,lonscale,nlat,latscale,
c     > 'pft','plant fuctional type','',npft,80,pftdef,
c     > tunits,'gregorian',istat)
c The '80' above refers to the length of the character strings in pftdef.
c         dimnames(3) = 'pft'
c         call inivar(idies,'exist','existence for each pft',
c     >    '',4,dimnames,OCEAN,istat)
c End initialization phase
c         call endini(idies,istat)
c
c 2) Initialize a variable and end initialization phase
c After you initialize a file, you need to initialize a variable.
c Initializing reserves space for data to be written into the file, and
c saves information about the data (such as units, descriptive name, and
c a missing value).  You can use inivar to initialize any variable of 1
c or more dimensions, provided that those dimensions already exist in
c the file.  You may initialize more than one variable in a single
c file, as stated above.
c If you have a special value that denotes 'missing', such as using
c 1.e+36 to denote ocean grid cells, use this value for valmissing.  If
c you use 0. for valmissing, it is ignored.  Pick a value for valmissing
c which is well outside the valid range for the data.  For example,
c -99. would be a fine missing value for pressure, but not a good value
c for topography.
c The character array dimnames is used to store the
c names of the dimensions upon which the new variable will depend.
c For example, a 3-d variable would only depend on longitude, latitude,
c and time.  A 4-d variable would depend on those and also pft, or
c snowlater, or level, or some other dimension.  The dimension names
c must be in the same order of varying: fastest to slowest.  See the
c discussions for istart and icount above.  For example, the 4-d
c variable to go in the file snow.nc (above inifile example)  would have
c dimnames as follows
c      dimnames(1) = 'longitude'
c      dimnames(2) = 'latitude'
c      dimnames(3) = 'snowlayer'
c      dimnames(4) = 'time'
c
c General call for inivar:
c      call inivar(idies,varname,longname,units,ndims,dimnames,
c     > valmissing,ierror)
c
c INPUT
c     idies - integer - id number of a new file from a previous call
c      to inifile, inifilec, or iesopen
c     varname - charcter*(*) - name for new variable
c     longname - character*(*) - descriptive name for variable
c     units - character*(*) - units for variable, compatible with udunits
c     ndims - integer - number of dimensions for this variable
c     dimnames - character*(*)(ndims) - name of each dimension, in order
c     valmissing - real - missing value, ignored if valmissing=0.
c OUTPUT
c     ierror - integer - error code, 0 = no error, < 0 = error occurred
c
c Example: Define a 4-d and a 3-d variable to go into the file snow.nc.
c defined previously
c        character*80 dimnames(4)
c        real OCEAN
c        dimnames(1) = 'longitude'
c        dimnames(2) = 'latitude'
c        OCEAN = 9.e+20
c define 4-d variable
c        dimnames(3) = 'snowlayer'
c        dimnames(4) = 'time'
c        call inivar(idies,'hsno','snow layer thickness',
c     >   'meters',4,dimnames,OCEAN,istat)
c define 3-d variable in same file
c        dimnames(3) = 'time'
c        call inivar(idies,'snowf','snow cover fraction',
c     >   'fraction',3,dimnames,OCEAN,istat)
c end initialization phase
c         call endini(idies,istat)
c Notice that you need to change the value of dimnames(3) for the 3-d
c variable.  The value in dimnames(4) gets ignored.
c
c 3) End the initialization phase for the file
c When you are done initializing the file and the variables in the file, you
c must call endini.  This subroutine merely closes the file.  But by closing
c the file, this is a signal to the computer to write out all changes to this
c file, synchronizing the instructions in the program with what is written on
c the hard disk.  Since most computers buffer, that is save up, data until a
c there is a large chunk to write, it is essential that all buffered
c information for a netcdf file be written to disk before attempting to write
c data to the file.  The subroutine endini accomplishes this.  The subroutine
c endini should be called after initializing the last variable in the netcdf
c file (inivar) and before calling writevar.  See the example above.
c
c 4) Write to a variable
c After you have completed the initialization phase by initializing the file
c and initializing one or more variables in the file, you can write data 
c into the space reserved for each variable with subroutine writevar.
c You need to supply istart and icount, just as when you read a
c variable.  It is perfectly legal to write values in only a portion of the
c space reserved for the variable.  If nothing is written to a
c particular grid point, it has a special fill value supplied by netcdf
c when the variable was initialized. Notice that for ease of use, I
c have designed writevar to also write in the values of the time steps
c you are writing as well as the weighting for that time step.
c Time weighting refers to the number of days in that time sample.
c For example, some monthly means will be a mean value over 31 days,
c while others will be over only 30, or 28, or 29 days.  I assumed that
c most persons will calculate data and write it out one time step at a
c time, so it is nice to write out the time values and weight as you
c go along.  Another time-related item written is the character label
c for the time step.  Since after a few years, it's hard to keep track
c of what month is represented by a large number of days, I invented the
c date variable.  Date is a 10-character long string in which you can
c put a time label.  For example the time value 4745 may not be
c informative, but the data label 'JAN1913   ' is.  I usually use 9
c characters plus a null character at the end (char(0)).  Use strings
c like 'DJF001005 ' to denote time averages (Winter years 1 through 5).
c
c General call
c     call writevar(filen,varname,istart,icount,values,
c    > times,tweights,dates,ierror)
c
c INPUT
c     filen - character*(*) - file name to which to write data
c     varname - character*(*) - name of variable to which to write
c     istart - integer(*) - starting points along each dimension of
c      the variable, for example the 1st point a 4-d variable is (1,1,1,1) 
c      and the 3rd level and 2nd time step of a 4d variable is (1,1,3,2).
c     icount - integer(*) - number of points along each dimension to write,
c      for example, to write a single lat/lon grid to a 4-d variable,
c      icount would be (nlon,nlat,1,1).
c     values - real(*) - real values of the designated hyperslab
c     times - real(*) - time vector vector for those points written (ignored
c      if variable is 2-d).
c     tweights - real(*) - number of days in each time step in times (ignored
c      if variable is 2-d).
c     dates - character*10(*) - character labels for each time step in times,
c      last character should be a null (dates ignored if variable 2-d).
c OUTPUT
c     ierror - integer - error code = 0 if no error, < 0 if error occurred
c
c Example: write data to the 4-d variable initialized above
c previously defined
c     parameter (ndim3=nlons*nlats*nsnolay)
c     real cdummy(ndim3), hsno(npts,nsnolay)
c     real time, timewght
c     integer istart(4), icount(4)
c     character*11 cdate
c     time = 730.
c     timewght = 365.
c     istart(1) = 1
c     istart(2) = 1
c     istart(3) = 1
c     istart(4) = 2
c     icount(1) = nlons
c     icount(2) = nlats
c     icount(3) = nsnolay
c     icount(4) = 1
c     cdate = 'ANN1980  '//char(0)
c put land-only data in hsno into land and sea lat/lon grid in workspace
c     do 20 k = 1, nsnolay
c        call vec2arr (hsno(1,k), cdummy((k-1)*nlonsub*nlatsub + 1))
c 20  continue
c     call writevar('snow.nc','hsno',istart,icount,cdummy,time,
c    > timewght,cdate,istat)
c     if (istat .ne. 0) then
c        write(*,*) 'ERROR writing hsno'
c        stop 1
c     end if
c
c Some other thoughts:
c It's a good idea to end all strings stored in netcdf file with a
c null character.  Some less-robust C programs may fail when
c encountering non-null-ended strings, because in C, all strings end
c with a null character.
c Udunits provides several ways to express multiplication, powers,
c etc.  I like to use ^ to denote exponents (as in m^2 instead of just
c m2) because programs like Matlab see ^ as a LaTeX command to write the
c next character as a superscript.  Likewise I try to avoid using _,
c because this is the LaTeX command to make the next character a
c subscript.  So instead of using deg_C, I use degC, to prevent the C
c from becoming a subscript.
c The format of the netcdf files written by ies-io.f subroutines is
c meant to be compatible with the COARDS and CSM conventions (see
c www.cgd.ucar.edu:80/csm/experiments/output/format.html).
c Theoretically, you should be able to use NCO (see
c www.cgd.ucar.edu:80/cms/nco/index.html) and NCL (see
c ngwww.ucar.edu/ngdoc/ng4.1alpha/ug/ncl/ncloview.html) on the files.
c A few other packages that work with this format are GrADS
c (grads.iges.org/grads/head.html) and Ncview
c (meteora.ucsd.edu/~pierce/ncview_home_page.html).
c The routines in ies-io.f are just the beginning.  You can put many
c more things in a netcdf file by using routines in ies.f or by using
c the low-level netcdf commands directly.  Use one of the pre-existing
c subroutines in ies-io.f as a template and go on from there.  You are
c free to change or distribute my code as long as you 1) acknowlege me,
c Christine C. Molling (cmolling@facstaff.wisc.edu) as the author of
c the original code and 2) do not sell it.
c Of course if there is anything wrong with the code, or you somehow
c encounter damage to your reputation/wallet because of its use, you
c are forbidden to sue me or anyone else.
c Good Luck!
