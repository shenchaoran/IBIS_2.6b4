c
c ------
c comsno
c ------
c
      real 
     >  z0sno,               ! roughness length of snow surface (m)
     >  rhos,                ! density of snow (kg m-3)
     >  consno,              ! thermal conductivity of snow (W m-1 K-1)
     >  hsnotop,             ! thickness of top snow layer (m)
     >  hsnomin,             ! minimum total thickness of snow (m)
     >  fimin,               ! minimum fractional snow cover
     >  fimax                ! maximum fractional snow cover
c
c     common /comsno1/ z0sno, rhos, consno, hsnotop, hsnomin, fimin, fimax
c
      real 
     >  fi(npoi)             ! fractional snow cover
c
c     common /comsno2/ fi
c
      real 
     >  tsno(npoi,nsnolay),  ! temperature of snow layers (K)
     >  hsno(npoi,nsnolay)   ! thickness of snow layers (m)
c
c     common /comsno3/ tsno, hsno
c
