c
c #    #    ##     #####  #    #
c ##  ##   #  #      #    #    #
c # ## #  #    #     #    ######
c #    #  ######     #    #    #
c #    #  #    #     #    #    #
c #    #  #    #     #    #    #
c
c
c --------------------------------------------------------------------
      real function ran2 (idum)
c --------------------------------------------------------------------
c
      include 'implicit.h'
c
      integer idum,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
c
      real am, eps, rnmx
c
      parameter (im1=2147483563, 
     >           im2=2147483399,
     >           am=1./im1,
     >           imm1=im1-1,
     >           ia1=40014,
     >           ia2=40692,
     >           iq1=53668,
     >           iq2=52774, 
     >           ir1=12211,
     >           ir2=3791,
     >           ntab=32,
     >           ndiv=1+imm1/ntab,
     >           eps=1.0e-7,
     >           rnmx=1.-eps)
c
      integer idum2,j,k,iv(ntab),iy
c
      save iv,iy,idum2
c
      data idum2/123456789/, iv/ntab*0/, iy/0/
c
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 10 j=ntab+8,1,-1
          k=idum/iq1
          idum=ia1*(idum-k*iq1)-k*ir1
          if (idum.lt.0) idum=idum+im1
          if (j.le.ntab) iv(j)=idum
 10     continue
        iy=iv(1)
      endif
c
      k=idum/iq1
      idum=ia1*(idum-k*iq1)-k*ir1
      if (idum.lt.0) idum=idum+im1
c
      k=idum2/iq2
      idum2=ia2*(idum2-k*iq2)-k*ir2
      if (idum2.lt.0) idum2=idum2+im2
c
      j=1+iy/ndiv
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+imm1
c
      ran2=min(am*iy,rnmx)
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine linsolve (arr, rhs, vec, mplate, nd)
c ---------------------------------------------------------------------
c
c solves multiple linear systems of equations, vectorizing
c over the number of systems. basic gaussian elimination is 
c used, with no pivoting (relies on all diagonal elements
c being and staying significantly non-zero)
c
c a template array mplate is used to detect when an operation 
c is not necessary (element already zero or would add zeros),
c assuming that every system has the same pattern of zero
c elements
c
c this template is first copied to mplatex since it 
c must be updated during the procedure in case an original-zero
c pattern location becomes non-zero
c
c the first subscript in arr, rhs, vec is over the multiple
c systems, and the others are the usual row, column subscripts
c
      include 'implicit.h'
c
      include 'compar.h'
c
c Arguments (input-output)
c
      integer nd                  ! number of equations (supplied)
c
      integer  mplate(nd,nd)      ! pattern of zero elements of arr (supplied)
c
      real arr(npoi,nd,nd),       ! equation coefficients (supplied, overwritten)
     >     rhs(npoi,nd),          ! equation right-hand sides (supplied, overwritten) 
     >     vec(npoi,nd)           ! solution (returned)
c 
c local variables
c
      integer ndx,                ! Max number of equations
     >        j, i, id, m         ! loop indices
c
      parameter (ndx=9)
c
      integer mplatex(ndx,ndx)
c
      real f(npoi)
c
      if (nd.gt.ndx) then
         write(*,900) nd, ndx
  900    format(/' *** fatal error ***'/
     >          /' number of linsolve eqns',i4,' exceeds limit',i4)
         call endrun
      endif
c
c copy the zero template so it can be changed below
c
      do 6 j=1,nd
        do 5 i=1,nd
          mplatex(i,j) = mplate(i,j)
    5   continue
    6 continue
c
c zero all array elements below the diagonal, proceeding from
c the first row to the last. note that mplatex is set non-zero
c for changed (i,j) locations, in loop 20
c
      do 10 id=1, nd-1
         do 12 i=id+1,nd
c
            if (mplatex(i,id).ne.0) then
               do 14 m=1,npoi
                  f(m) = arr(m,i,id) / arr(m,id,id)
   14          continue
c
               do 20 j=id,nd
                  if (mplatex(id,j).ne.0) then
                     do 22 m=1,npoi
                        arr(m,i,j) = arr(m,i,j) - f(m)*arr(m,id,j)
   22                continue
                     mplatex(i,j) = 1
                  endif
   20          continue
c
               do 30 m=1,npoi
                  rhs(m,i) = rhs(m,i) - f(m)*rhs(m,id)
   30          continue
            endif
c
   12    continue
   10 continue
c
c all array elements below the diagonal are zero, so can
c immediately solve the equations in reverse order
c
      do 50 id=nd,1,-1
c
         call const (f, npoi, 0.0)
         if (id.lt.nd) then
            do 52 j=id+1,nd
               if (mplatex(id,j).ne.0) then
                  do 54 m=1,npoi
                     f(m) = f(m) + arr(m,id,j)*vec(m,j)
   54             continue
               endif
   52       continue
         endif
c
         do 56 m=1,npoi
            vec(m,id) = (rhs(m,id) - f(m)) / arr(m,id,id)
   56    continue
c
   50 continue
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine tridia (ns, nd, ne, a, b, c, y, x, alpha, gamma)
c ---------------------------------------------------------------------
c
      include 'implicit.h'
c
      include 'compar.h'
c
c     purpose:
c     to compute the solution of many tridiagonal linear systems.
c
c      arguments:
c
c      ns ..... the number of systems to be solved.
c
c      nd ..... first dimension of arrays (ge ns).
c
c      ne ..... the number of unknowns in each system.
c               this must be > 2. second dimension of arrays.
c
c      a ...... the subdiagonals of the matrices are stored
c               in locations a(j,2) through a(j,ne).
c
c      b ...... the main diagonals of the matrices are stored
c               in locations b(j,1) through b(j,ne).
c
c      c ...... the super-diagonals of the matrices are stored in
c               locations c(j,1) through c(j,ne-1).
c
c      y ...... the right hand side of the equations is stored in
c               y(j,1) through y(j,ne).
c
c      x ...... the solutions of the systems are returned in
c               locations x(j,1) through x(j,ne).
c
c      alpha .. work array dimensioned alpha(nd,ne)
c
c      gamma .. work array dimensioned gamma(nd,ne)
c
c       history:  based on a streamlined version of the old ncar
c                 ulib subr trdi used in the phoenix climate
c                 model of schneider and thompson (j.g.r., 1981).
c                 revised by starley thompson to solve multiple
c                 systems and vectorize well on the cray-1.
c                 later revised to include a parameter statement
c                 to define loop limits and thus enable cray short
c                 vector loops.
c
c       algorithm:  lu decomposition followed by solution.
c                   note: this subr executes satisfactorily
c                   if the input matrix is diagonally dominant
c                   and non-singular.  the diagonal elements are
c                   used to pivot, and no tests are made to determine
c                   singularity. if a singular or numerically singular
c                   matrix is used as input a divide by zero or
c                   floating point overflow will result.
c
c       last revision date:      4 february 1988
c
c
c Arguments
c
      integer ns,     ! number of systems to be solved.
     >        nd,     ! first dimension of arrays (ge ns)
     >        ne      ! number of unknowns in each system. (>2)
      
      real 
     >  a(nd,ne),     ! subdiagonals of matrices stored in a(j,2)...a(j,ne).
     >  b(nd,ne),     ! main diagonals of matrices stored in b(j,1)...b(j,ne).
     >  c(nd,ne),     ! super-diagonals of matrices stored in c(j,1)...c(j,ne-1).
     >  y(nd,ne),     ! right hand side of equations stored in y(j,1)...y(j,ne).
     >  x(nd,ne),     ! solutions of the systems returned in x(j,1)...x(j,ne).
     >  alpha(nd,ne), ! work array 
     >  gamma(nd,ne)  ! work array
c
c local variables
c
      integer nm1,    !
     >  j, i, ib      ! loop indices

c
      nm1 = ne-1
c
c obtain the lu decompositions
c
      do 10 j=1,ns
         alpha(j,1) = 1./b(j,1)
         gamma(j,1) = c(j,1)*alpha(j,1)
   10 continue
      do 11 i=2,nm1
         do 12 j=1,ns
            alpha(j,i) = 1./(b(j,i)-a(j,i)*gamma(j,i-1))
            gamma(j,i) = c(j,i)*alpha(j,i)
   12    continue
   11 continue
c
c solve
c
      do 20 j=1,ns
         x(j,1) = y(j,1)*alpha(j,1)
   20 continue
      do 21 i=2,nm1
         do 22 j=1,ns
            x(j,i) = (y(j,i)-a(j,i)*x(j,i-1))*alpha(j,i)
   22    continue
   21 continue
      do 23 j=1,ns
         x(j,ne) = (y(j,ne)-a(j,ne)*x(j,nm1))/
     >             (b(j,ne)-a(j,ne)*gamma(j,nm1))
   23 continue
      do 24 i=1,nm1
         ib = ne-i
         do 25 j=1,ns
            x(j,ib) = x(j,ib)-gamma(j,ib)*x(j,ib+1)
   25    continue
   24 continue
c
      return
      end
c
