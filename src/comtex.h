c ------
c comtex
c ------
c
c ----------------------------------------
c Soil texture-related parameters for ibis
c ----------------------------------------

      integer*2 ndat        ! number of soil types
      parameter (ndat = 11) ! excludes organics for now.


      real*4 
     >       texdat(3, ndat),  ! sand/silt/clay fractions
     >       porosdat(ndat),   ! porosity volume fraction
     >       sfielddat(ndat),  ! field capacity volume fraction
     >       swiltdat(ndat),   ! wilting point volume fraction
     >       bexdat(ndat),     ! Campbell moisture-release b exponent
     >       suctiondat(ndat), ! Air entry potential (m-H20)
     >       hydrauldat(ndat)  ! saturated hydraulic conductivity (m s-1)

      common /comtex1/
     >       texdat,           ! sand/silt/clay fractions
     >       porosdat,         ! porosity volume fraction
     >       sfielddat,        ! field capacity volume fraction
     >       swiltdat,         ! wilting point volume fraction
     >       bexdat,           ! Campbell moisture-release b exponent
     >       suctiondat,       ! Air entry potential (m-H20)
     >       hydrauldat        ! saturated hydraulic conductivity (m s-1)
