c
c ------
c comsoi
c ------
c
      real 
     >  wpudmax,                 ! normalization constant for puddles (kg m-2)
     >  zwpmax,                  ! assumed maximum fraction of soil surface 
*                                ! covered by puddles (dimensionless)
     >  bperm                    ! lower b.c. for soil profile drainage 
*                                ! (0.0 = impermeable; 1.0 = fully permeable)

c
      common /comsoi1/ wpudmax, zwpmax, bperm
c
      real 
     >  wpud(npoi),              ! liquid content of puddles per soil area (kg m-2)
     >  wipud(npoi),             ! ice content of puddles per soil area (kg m-2)
     >  z0soi(npoi),             ! roughness length of soil surface (m)
     >  albsav(npoi),            ! saturated soil surface albedo (visible waveband)
     >  albsan(npoi),            ! saturated soil surface albedo (near-ir waveband)
     >  stresstl(npoi),          ! sum of stressl over all 6 soil layers (dimensionless)
     >  stresstu(npoi),          ! sum of stressu over all 6 soil layers (dimensionless)
     >  heati(npoi),             ! net heat flux into snow surface (W m-2)
     >  heatg(npoi),             ! net heat flux into soil surface (W m-2)
     >  hvasug(npoi),            ! latent heat of vap/subl, for soil surface (J kg-1)
     >  hvasui(npoi),            ! latent heat of vap/subl, for snow surface (J kg-1)
     >  tg(npoi),                ! soil skin temperature (K)
     >  ti(npoi)                 ! snow skin temperature (K)
c
      common /comsoi2/ wpud, wipud, z0soi, albsav, albsan, stresstl, stresstu, heati,
     >     heatg, hvasug, hvasui, tg, ti
c
      real
     >  hsoi(nsoilay+1)          ! soil layer thickness (m)
c
      common /comsoi4/ hsoi
c
      real 
     >  tsoi(npoi,nsoilay),      ! soil temperature for each layer (K)
     >  wsoi(npoi,nsoilay),      ! fraction of soil pore space containing liquid water
     >  wisoi(npoi,nsoilay),     ! fraction of soil pore space containing ice
     >  consoi(npoi,nsoilay),    ! thermal conductivity of each soil layer (W m-1 K-1)
     >  csoi(npoi,nsoilay),      ! specific heat of soil, no pore spaces (J kg-1 deg-1)
     >  hydraul(npoi,nsoilay),   ! saturated hydraulic conductivity (m/s)
     >  suction(npoi,nsoilay),   ! saturated matric potential (m-h2o)
     >  bex(npoi,nsoilay),       ! exponent "b" in soil water potential
     >  sfield(npoi,nsoilay),    ! field capacity soil moisture value (fraction of pore space)
     >  swilt(npoi,nsoilay),     ! wilting soil moisture value (fraction of pore space)
     >  rhosoi(npoi,nsoilay),    ! soil density (without pores, not bulk) (kg m-3)
     >  poros(npoi,nsoilay),     ! porosity (mass of h2o per unit vol at sat / rhow)
     >  porosflo(npoi,nsoilay),  ! porosity after reduction by ice content
     >  sandpct(npoi,nsoilay),   ! percent sand of soil from input file
     >  claypct(npoi,nsoilay),   ! percent clay of soil from input file
     >  sand(npoi,nsoilay),      ! percent sand of soil
     >  clay(npoi,nsoilay),      ! percent clay of soil
     >  stressl(npoi,nsoilay),   ! soil moisture stress factor for the lower canopy (dimensionless)
     >  stressu(npoi,nsoilay),   ! soil moisture stress factor for the upper canopy (dimensionless)
     >  upsoiu(npoi,nsoilay),    ! soil water uptake from transpiration (kg_h2o m-2 s-1)
     >  upsoil(npoi,nsoilay)     ! soil water uptake from transpiration (kg_h2o m-2 s-1)
c
      common /comsoi5/ tsoi, wsoi, wisoi, consoi, csoi, hydraul, suction, bex,
     >     sfield, swilt, rhosoi, poros, porosflo, sandpct, claypct, 
     >     sand, clay, stressl, stressu, upsoiu, upsoil
c
      real 
     >  hflo(npoi,nsoilay+1)     ! downward heat transport through soil layers (W m-2)
c
      common /comsoi6/ hflo
c
      integer 
     >  ibex(npoi,nsoilay)       ! nint(bex), used for cpu speed
c
      common /comsoi7/ ibex
c
      real 
     >  qglif(npoi,4)            ! 1: fraction of soil evap (fvapg) from soil liquid
     >                           ! 2: fraction of soil evap (fvapg) from soil ice
     >                           ! 3: fraction of soil evap (fvapg) from puddle liquid
     >                           ! 4: fraction of soil evap (fvapg) from puddle ice
c
      common /comsoi8/ qglif
c
